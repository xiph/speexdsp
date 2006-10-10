/* Copyright (C) 2003 Epic Games 
   Written by Jean-Marc Valin

   File: preprocess.c
   Preprocessor with denoising based on the algorithm by Ephraim and Malah

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:

   1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

   3. The name of the author may not be used to endorse or promote products
   derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
   IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
   OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
   INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
   HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
   STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include "speex/speex_preprocess.h"
#include "misc.h"
#include "fftwrap.h"
#include "filterbank.h"

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))

#ifndef M_PI
#define M_PI 3.14159263
#endif

#define SQRT_M_PI_2 0.88623
#define LOUDNESS_EXP 2.5

#define NB_BANDS 8

#define SPEEX_PROB_START_DEFAULT    0.35f
#define SPEEX_PROB_CONTINUE_DEFAULT 0.20f

#define ZMIN .1
#define ZMAX .316
#define ZMIN_1 10
#define LOG_MIN_MAX_1 0.86859


/** Speex pre-processor state. */
struct SpeexPreprocessState_ {
   int    frame_size;        /**< Number of samples processed each time */
   int    ps_size;           /**< Number of points in the power spectrum */
   int    sampling_rate;     /**< Sampling rate of the input/output */
   
   /* parameters */
   int    denoise_enabled;
   int    agc_enabled;
   float  agc_level;
   int    vad_enabled;
   int    dereverb_enabled;
   float  reverb_decay;
   float  reverb_level;
   float  speech_prob_start;
   float  speech_prob_continue;
   
   FilterBank *bank;
   
   float *frame;             /**< Processing frame (2*ps_size) */
   float *ft;                /**< Processing frame in freq domain (2*ps_size) */
   float *ps;                /**< Current power spectrum */
   float *gain2;             /**< Adjusted gains */
   float *window;            /**< Analysis/Synthesis window */
   float *noise;             /**< Noise estimate */
   float *reverb_estimate;   /**< Estimate of reverb energy */
   float *old_ps;            /**< Power spectrum for last frame */
   float *gain;              /**< Ephraim Malah gain */
   float *prior;             /**< A-priori SNR */
   float *post;              /**< A-posteriori SNR */

   float *S;                 /**< Smoothed power spectrum */
   float *Smin;              /**< See Cohen paper */
   float *Stmp;              /**< See Cohen paper */
   float *update_prob;       /**< Propability of speech presence for noise update */

   float *zeta;              /**< Smoothed a priori SNR */
   float  Zpeak;
   float  Zlast;

   float *loudness_weight;   /**< Perceptual loudness curve */

   spx_int32_t *echo_noise;

   float *noise_bands;
   float *noise_bands2;
   int    noise_bandsN;
   float *speech_bands;
   float *speech_bands2;
   int    speech_bandsN;

   float *inbuf;             /**< Input buffer (overlapped analysis) */
   float *outbuf;            /**< Output buffer (for overlap and add) */

   float  speech_prob;
   int    last_speech;
   float  loudness;          /**< loudness estimate */
   float  loudness2;         /**< loudness estimate */
   int    nb_adapt;          /**< Number of frames used for adaptation so far */
   int    nb_loudness_adapt; /**< Number of frames used for loudness adaptation so far */
   int    consec_noise;      /**< Number of consecutive noise frames */
   int    nb_preprocess;     /**< Number of frames processed so far */
   void *fft_lookup;   /**< Lookup table for the FFT */

};


static void conj_window(float *w, int len)
{
   int i;
   for (i=0;i<len;i++)
   {
      float x=4*((float)i)/len;
      int inv=0;
      if (x<1)
      {
      } else if (x<2)
      {
         x=2-x;
         inv=1;
      } else if (x<3)
      {
         x=x-2;
         inv=1;
      } else {
         x=4-x;
      }
      x*=1.9979;
      w[i]=(.5-.5*cos(x))*(.5-.5*cos(x));
      if (inv)
         w[i]=1-w[i];
      w[i]=sqrt(w[i]);
   }
}

/* This function approximates the gain function 
   y = gamma(1.25)^2 * M(-.25;1;-x) / sqrt(x)  
   which multiplied by xi/(1+xi) is the optimal gain
   in the loudness domain ( sqrt[amplitude] )
*/
static inline float hypergeom_gain(float x)
{
   int ind;
   float integer, frac;
   static const float table[21] = {
      0.82157f, 1.02017f, 1.20461f, 1.37534f, 1.53363f, 1.68092f, 1.81865f,
      1.94811f, 2.07038f, 2.18638f, 2.29688f, 2.40255f, 2.50391f, 2.60144f,
      2.69551f, 2.78647f, 2.87458f, 2.96015f, 3.04333f, 3.12431f, 3.20326f};
      
   integer = floor(2*x);
   ind = (int)integer;
   if (ind<0)
      return 1;
   if (ind>19)
      return 1+.1296/x;
   frac = 2*x-integer;
   return ((1-frac)*table[ind] + frac*table[ind+1])/sqrt(x+.0001f);
}

static inline float qcurve(float x)
{
   return 1.f/(1.f+.1f/(x*x));
}

SpeexPreprocessState *speex_preprocess_state_init(int frame_size, int sampling_rate)
{
   int i;
   int N, N3, N4;

   SpeexPreprocessState *st = (SpeexPreprocessState *)speex_alloc(sizeof(SpeexPreprocessState));
   st->frame_size = frame_size;

   /* Round ps_size down to the nearest power of two */
#if 0
   i=1;
   st->ps_size = st->frame_size;
   while(1)
   {
      if (st->ps_size & ~i)
      {
         st->ps_size &= ~i;
         i<<=1;
      } else {
         break;
      }
   }
   
   
   if (st->ps_size < 3*st->frame_size/4)
      st->ps_size = st->ps_size * 3 / 2;
#else
   st->ps_size = st->frame_size;
#endif

   N = st->ps_size;
   N3 = 2*N - st->frame_size;
   N4 = st->frame_size - N3;
   
   st->sampling_rate = sampling_rate;
   st->denoise_enabled = 1;
   st->agc_enabled = 0;
   st->agc_level = 8000;
   st->vad_enabled = 0;
   st->dereverb_enabled = 0;
   st->reverb_decay = .5;
   st->reverb_level = .2;

   st->speech_prob_start = SPEEX_PROB_START_DEFAULT;
   st->speech_prob_continue = SPEEX_PROB_CONTINUE_DEFAULT;

   st->bank = filterbank_new(24, 4000, 8000, N, 1);
         
   st->frame = (float*)speex_alloc(2*N*sizeof(float));
   st->ft = (float*)speex_alloc(2*N*sizeof(float));
   st->ps = (float*)speex_alloc(N*sizeof(float));
   st->gain2 = (float*)speex_alloc(N*sizeof(float));
   st->window = (float*)speex_alloc(2*N*sizeof(float));
   st->noise = (float*)speex_alloc(N*sizeof(float));
   st->reverb_estimate = (float*)speex_alloc(N*sizeof(float));
   st->old_ps = (float*)speex_alloc(N*sizeof(float));
   st->gain = (float*)speex_alloc(N*sizeof(float));
   st->prior = (float*)speex_alloc(N*sizeof(float));
   st->post = (float*)speex_alloc(N*sizeof(float));
   st->loudness_weight = (float*)speex_alloc(N*sizeof(float));
   st->inbuf = (float*)speex_alloc(N3*sizeof(float));
   st->outbuf = (float*)speex_alloc(N3*sizeof(float));
   st->echo_noise = (spx_int32_t*)speex_alloc(N*sizeof(float));

   st->S = (float*)speex_alloc(N*sizeof(float));
   st->Smin = (float*)speex_alloc(N*sizeof(float));
   st->Stmp = (float*)speex_alloc(N*sizeof(float));
   st->update_prob = (float*)speex_alloc(N*sizeof(float));

   st->zeta = (float*)speex_alloc(N*sizeof(float));
   st->Zpeak = 0;
   st->Zlast = 0;

   st->noise_bands = (float*)speex_alloc(NB_BANDS*sizeof(float));
   st->noise_bands2 = (float*)speex_alloc(NB_BANDS*sizeof(float));
   st->speech_bands = (float*)speex_alloc(NB_BANDS*sizeof(float));
   st->speech_bands2 = (float*)speex_alloc(NB_BANDS*sizeof(float));
   st->noise_bandsN = st->speech_bandsN = 1;

   conj_window(st->window, 2*N3);
   for (i=2*N3;i<2*st->ps_size;i++)
      st->window[i]=1;
   
   if (N4>0)
   {
      for (i=N3-1;i>=0;i--)
      {
         st->window[i+N3+N4]=st->window[i+N3];
         st->window[i+N3]=1;
      }
   }
   for (i=0;i<N;i++)
   {
      st->noise[i]=1e4;
      st->reverb_estimate[i]=0.;
      st->old_ps[i]=1e4;
      st->gain[i]=1;
      st->post[i]=1;
      st->prior[i]=1;
   }

   for (i=0;i<N3;i++)
   {
      st->inbuf[i]=0;
      st->outbuf[i]=0;
   }

   for (i=0;i<N;i++)
   {
      float ff=((float)i)*.5*sampling_rate/((float)N);
      st->loudness_weight[i] = .35f-.35f*ff/16000.f+.73f*exp(-.5f*(ff-3800)*(ff-3800)/9e5f);
      if (st->loudness_weight[i]<.01f)
         st->loudness_weight[i]=.01f;
      st->loudness_weight[i] *= st->loudness_weight[i];
   }

   st->speech_prob = 0;
   st->last_speech = 1000;
   st->loudness = pow(6000,LOUDNESS_EXP);
   st->loudness2 = 6000;
   st->nb_loudness_adapt = 0;

   st->fft_lookup = spx_fft_init(2*N);

   st->nb_adapt=0;
   st->consec_noise=0;
   st->nb_preprocess=0;
   return st;
}

void speex_preprocess_state_destroy(SpeexPreprocessState *st)
{
   speex_free(st->frame);
   speex_free(st->ft);
   speex_free(st->ps);
   speex_free(st->gain2);
   speex_free(st->window);
   speex_free(st->noise);
   speex_free(st->reverb_estimate);
   speex_free(st->old_ps);
   speex_free(st->gain);
   speex_free(st->prior);
   speex_free(st->post);
   speex_free(st->loudness_weight);
   speex_free(st->echo_noise);

   speex_free(st->S);
   speex_free(st->Smin);
   speex_free(st->Stmp);
   speex_free(st->update_prob);
   speex_free(st->zeta);

   speex_free(st->noise_bands);
   speex_free(st->noise_bands2);
   speex_free(st->speech_bands);
   speex_free(st->speech_bands2);

   speex_free(st->inbuf);
   speex_free(st->outbuf);

   spx_fft_destroy(st->fft_lookup);

   speex_free(st);
}

static void speex_compute_agc(SpeexPreprocessState *st)
{
   int i;
   int N = st->ps_size;
   float scale=.5f/N;
   float agc_gain;
   int freq_start, freq_end;
   float active_bands = 0;

   freq_start = (int)(300.0f*2*N/st->sampling_rate);
   freq_end   = (int)(2000.0f*2*N/st->sampling_rate);
   for (i=freq_start;i<freq_end;i++)
   {
      if (st->S[i] > 20.f*st->Smin[i]+1000.f)
         active_bands+=1;
   }
   active_bands /= (freq_end-freq_start+1);

   if (active_bands > .2f)
   {
      float loudness=0.f;
      float rate, rate2=.2f;
      st->nb_loudness_adapt++;
      rate=2.0f/(1+st->nb_loudness_adapt);
      if (rate < .05f)
         rate = .05f;
      if (rate < .1f && pow(loudness, LOUDNESS_EXP) > st->loudness)
         rate = .1f;
      if (rate < .2f && pow(loudness, LOUDNESS_EXP) > 3.f*st->loudness)
         rate = .2f;
      if (rate < .4f && pow(loudness, LOUDNESS_EXP) > 10.f*st->loudness)
         rate = .4f;

      for (i=2;i<N;i++)
      {
         loudness += scale*st->ps[i] * st->gain2[i] * st->gain2[i] * st->loudness_weight[i];
      }
      loudness=sqrt(loudness);
      /*if (loudness < 2*pow(st->loudness, 1.0/LOUDNESS_EXP) &&
        loudness*2 > pow(st->loudness, 1.0/LOUDNESS_EXP))*/
      st->loudness = (1-rate)*st->loudness + (rate)*pow(loudness, LOUDNESS_EXP);
      
      st->loudness2 = (1-rate2)*st->loudness2 + rate2*pow(st->loudness, 1.0f/LOUDNESS_EXP);

      loudness = pow(st->loudness, 1.0f/LOUDNESS_EXP);

      /*fprintf (stderr, "%f %f %f\n", loudness, st->loudness2, rate);*/
   }
   
   agc_gain = st->agc_level/st->loudness2;
   /*fprintf (stderr, "%f %f %f %f\n", active_bands, st->loudness, st->loudness2, agc_gain);*/
   if (agc_gain>200)
      agc_gain = 200;

   for (i=0;i<N;i++)
      st->gain2[i] *= agc_gain;
   
}

static void preprocess_analysis(SpeexPreprocessState *st, spx_int16_t *x)
{
   int i;
   int N = st->ps_size;
   int N3 = 2*N - st->frame_size;
   int N4 = st->frame_size - N3;
   float *ps=st->ps;

   /* 'Build' input frame */
   for (i=0;i<N3;i++)
      st->frame[i]=st->inbuf[i];
   for (i=0;i<st->frame_size;i++)
      st->frame[N3+i]=x[i];
   
   /* Update inbuf */
   for (i=0;i<N3;i++)
      st->inbuf[i]=x[N4+i];

   /* Windowing */
   for (i=0;i<2*N;i++)
      st->frame[i] *= st->window[i];

   /* Perform FFT */
   spx_fft_float(st->fft_lookup, st->frame, st->ft);
         
   /* Power spectrum */
   ps[0]=1;
   for (i=1;i<N;i++)
      ps[i]=1+st->ft[2*i-1]*st->ft[2*i-1] + st->ft[2*i]*st->ft[2*i];

}

static void update_noise_prob(SpeexPreprocessState *st)
{
   int i;
   int N = st->ps_size;

   for (i=1;i<N-1;i++)
      st->S[i] = 100.f+ .8f*st->S[i] + .05f*st->ps[i-1]+.1f*st->ps[i]+.05f*st->ps[i+1];
   
   if (st->nb_preprocess<1)
   {
      for (i=1;i<N-1;i++)
         st->Smin[i] = st->Stmp[i] = st->S[i]+100.f;
   }

   if (st->nb_preprocess%200==0)
   {
      for (i=1;i<N-1;i++)
      {
         st->Smin[i] = min(st->Stmp[i], st->S[i]);
         st->Stmp[i] = st->S[i];
      }
   } else {
      for (i=1;i<N-1;i++)
      {
         st->Smin[i] = min(st->Smin[i], st->S[i]);
         st->Stmp[i] = min(st->Stmp[i], st->S[i]);      
      }
   }
   for (i=1;i<N-1;i++)
   {
      st->update_prob[i] *= .2f;
      if (st->S[i] > 2.5*st->Smin[i])
         st->update_prob[i] += .8f;
      /*fprintf (stderr, "%f ", st->S[i]/st->Smin[i]);*/
      /*fprintf (stderr, "%f ", st->update_prob[i]);*/
   }

}

#define NOISE_OVERCOMPENS 1.4

int speex_preprocess(SpeexPreprocessState *st, spx_int16_t *x, spx_int32_t *echo)
{
   int i;
   int N = st->ps_size;
   int N3 = 2*N - st->frame_size;
   int N4 = st->frame_size - N3;
   /*float scale=.5f/N;*/
   float *ps=st->ps;
   float Zframe=0, Pframe;
   float beta, beta_1;
   
   st->nb_adapt++;
   beta =1.0f/st->nb_adapt;
   if (beta < .05f)
      beta=.05f;
   beta_1 = 1.0f-beta;
   
   /* Deal with residual echo if provided */
   if (echo)
      for (i=0;i<N;i++)
         st->echo_noise[i] = (.3f*st->echo_noise[i] + echo[i]);
   else
      for (i=0;i<N;i++)
         st->echo_noise[i] = 0;

   preprocess_analysis(st, x);

   update_noise_prob(st);

   st->nb_preprocess++;

   /* Noise estimation always updated for the 10 first frames */
   if (st->nb_adapt<10)
   {
      for (i=1;i<N-1;i++)
         st->update_prob[i] = 0;
   }

   for (i=1;i<N-1;i++)
   {
      if (st->update_prob[i]<.5f/* || st->ps[i] < st->noise[i]*/)
         st->noise[i] = beta_1*st->noise[i] + beta*st->ps[i];
   }

   /* Special case for first frame */
   if (st->nb_adapt==1)
      for (i=0;i<N;i++)
         st->old_ps[i] = ps[i];

   /* Compute a posteriori SNR */
   for (i=1;i<N;i++)
   {
      float gamma;
      float tot_noise = 1.f+ NOISE_OVERCOMPENS*st->noise[i] + st->echo_noise[i] + st->reverb_estimate[i];
      st->post[i] = ps[i]/tot_noise - 1.f;
      if (st->post[i]>100.f)
         st->post[i]=100.f;
      gamma = .15+.85*st->prior[i]*st->prior[i]/((1+st->prior[i])*(1+st->prior[i]));
      /* A priori SNR update */
      st->prior[i] = gamma*max(0.0f,st->post[i]) +
            (1.f-gamma)* (.8*st->gain[i]*st->gain[i]*st->old_ps[i]/tot_noise + .2*st->prior[i]);
      if (st->prior[i]>100.f)
         st->prior[i]=100.f;
   }

#if 0
   for (i=0;i<N;i++)
   {
      fprintf (stderr, "%f ", st->prior[i]);
   }
   fprintf (stderr, "\n");
#endif


   for (i=1;i<N;i++)
   {
      st->zeta[i] = .7f*st->zeta[i] + .3f*st->prior[i];
   }

   {
      int freq_start = (int)(300.0f*2.f*N/st->sampling_rate);
      int freq_end   = (int)(2000.0f*2.f*N/st->sampling_rate);
      for (i=freq_start;i<freq_end;i++)
      {
         Zframe += st->zeta[i];         
      }
      Zframe /= (freq_end-freq_start);
   }
   st->Zlast = Zframe;

   Pframe = qcurve(Zframe);

   /*fprintf (stderr, "%f\n", Pframe);*/
   /* Compute gain according to the Ephraim-Malah algorithm */
   for (i=1;i<N;i++)
   {
      float MM;
      float theta;
      float prior_ratio;
      float p, q;
      float zeta1;
      float P1;

      prior_ratio = st->prior[i]/(1.0001f+st->prior[i]);
      theta = (1.f+st->post[i])*prior_ratio;

      if (i==1 || i==N-1)
         zeta1 = st->zeta[i];
      else
         zeta1 = .25f*st->zeta[i-1] + .5f*st->zeta[i] + .25f*st->zeta[i+1];
      P1 = qcurve (zeta1);
      
      /* FIXME: add global prob (P2) */
      q = 1-Pframe*P1;
      q = 1-P1;
      if (q>.95f)
         q=.95f;
      p=1.f/(1.f + (q/(1.f-q))*(1.f+st->prior[i])*exp(-theta));
      /*p=1;*/

      /* Optimal estimator for loudness domain */
      MM = hypergeom_gain(theta);

      st->gain[i] = prior_ratio * MM;
      /*Put some (very arbitraty) limit on the gain*/
      if (st->gain[i]>2.f)
      {
         st->gain[i]=2.f;
      }
      
      st->reverb_estimate[i] = st->reverb_decay*st->reverb_estimate[i] + st->reverb_decay*st->reverb_level*st->gain[i]*st->gain[i]*st->ps[i];
      if (st->denoise_enabled)
      {
         /*st->gain2[i] = p*p*st->gain[i];*/
         st->gain2[i]=(p*sqrt(st->gain[i])+.2*(1-p)) * (p*sqrt(st->gain[i])+.2*(1-p));
         /*st->gain2[i] = pow(st->gain[i], p) * pow(.1f,1.f-p);*/
      } else {
         st->gain2[i]=1.f;
      }
   }
   
   st->gain2[0]=st->gain[0]=0.f;
   st->gain2[N-1]=st->gain[N-1]=0.f;
   
   if (1) {
      float m[24];
      filterbank_compute_bank(st->bank, st->gain2, m);
      filterbank_compute_psd(st->bank,m, st->gain2);
   }
   /*
   for (i=30;i<N-2;i++)
   {
      st->gain[i] = st->gain2[i]*st->gain2[i] + (1-st->gain2[i])*.333*(.6*st->gain2[i-1]+st->gain2[i]+.6*st->gain2[i+1]+.4*st->gain2[i-2]+.4*st->gain2[i+2]);
   }
   for (i=30;i<N-2;i++)
      st->gain2[i] = st->gain[i];
   */
   if (st->agc_enabled)
      speex_compute_agc(st);

#if 0
   if (!is_speech)
   {
      for (i=0;i<N;i++)
         st->gain2[i] = 0;
   }
#if 0
 else {
      for (i=0;i<N;i++)
         st->gain2[i] = 1;
   }
#endif
#endif

   /* Apply computed gain */
   for (i=1;i<N;i++)
   {
      st->ft[2*i-1] *= st->gain2[i];
      st->ft[2*i] *= st->gain2[i];
   }

   /* Get rid of the DC and very low frequencies */
   st->ft[0]=0;
   st->ft[1]=0;
   st->ft[2]=0;
   /* Nyquist frequency is mostly useless too */
   st->ft[2*N-1]=0;

   /* Inverse FFT with 1/N scaling */
   spx_ifft_float(st->fft_lookup, st->ft, st->frame);
         
   /*for (i=0;i<2*N;i++)
    st->frame[i] *= scale;*/

   {
      float max_sample=0;
      for (i=0;i<2*N;i++)
         if (fabs(st->frame[i])>max_sample)
            max_sample = fabs(st->frame[i]);
      if (max_sample>28000.f)
      {
         float damp = 28000.f/max_sample;
         for (i=0;i<2*N;i++)
            st->frame[i] *= damp;
      }
   }

   for (i=0;i<2*N;i++)
      st->frame[i] *= st->window[i];

   /* Perform overlap and add */
   for (i=0;i<N3;i++)
      x[i] = st->outbuf[i] + st->frame[i];
   for (i=0;i<N4;i++)
      x[N3+i] = st->frame[N3+i];
   
   /* Update outbuf */
   for (i=0;i<N3;i++)
      st->outbuf[i] = st->frame[st->frame_size+i];

   /* Save old power spectrum */
   for (i=1;i<N;i++)
      st->old_ps[i] = ps[i];

   return 1;
}

void speex_preprocess_estimate_update(SpeexPreprocessState *st, spx_int16_t *x, spx_int32_t *echo)
{
   int i;
   int N = st->ps_size;
   int N3 = 2*N - st->frame_size;

   float *ps=st->ps;

   preprocess_analysis(st, x);

   update_noise_prob(st);

   st->nb_preprocess++;
   
   for (i=1;i<N-1;i++)
   {
      if (st->update_prob[i]<.5f || st->ps[i] < st->noise[i])
      {
         if (echo)
            st->noise[i] = .95f*st->noise[i] + .1f*max(1.0f,st->ps[i]-1.0*echo[i]);
         else
            st->noise[i] = .95f*st->noise[i] + .1f*st->ps[i];
      }
   }

   for (i=0;i<N3;i++)
      st->outbuf[i] = x[st->frame_size-N3+i]*st->window[st->frame_size+i];

   /* Save old power spectrum */
   for (i=1;i<N;i++)
      st->old_ps[i] = ps[i];

   for (i=1;i<N;i++)
      st->reverb_estimate[i] *= st->reverb_decay;
}


int speex_preprocess_ctl(SpeexPreprocessState *state, int request, void *ptr)
{
   int i;
   SpeexPreprocessState *st;
   st=(SpeexPreprocessState*)state;
   switch(request)
   {
   case SPEEX_PREPROCESS_SET_DENOISE:
      st->denoise_enabled = (*(int*)ptr);
      break;
   case SPEEX_PREPROCESS_GET_DENOISE:
      (*(int*)ptr) = st->denoise_enabled;
      break;

   case SPEEX_PREPROCESS_SET_AGC:
      st->agc_enabled = (*(int*)ptr);
      break;
   case SPEEX_PREPROCESS_GET_AGC:
      (*(int*)ptr) = st->agc_enabled;
      break;

   case SPEEX_PREPROCESS_SET_AGC_LEVEL:
      st->agc_level = (*(float*)ptr);
      if (st->agc_level<1)
         st->agc_level=1;
      if (st->agc_level>32768)
         st->agc_level=32768;
      break;
   case SPEEX_PREPROCESS_GET_AGC_LEVEL:
      (*(float*)ptr) = st->agc_level;
      break;

   case SPEEX_PREPROCESS_SET_VAD:
      speex_warning("The VAD has been removed pending a complete rewrite");
      st->vad_enabled = (*(int*)ptr);
      break;
   case SPEEX_PREPROCESS_GET_VAD:
      (*(int*)ptr) = st->vad_enabled;
      break;
   
   case SPEEX_PREPROCESS_SET_DEREVERB:
      st->dereverb_enabled = (*(int*)ptr);
      for (i=0;i<st->ps_size;i++)
         st->reverb_estimate[i]=0;
      break;
   case SPEEX_PREPROCESS_GET_DEREVERB:
      (*(int*)ptr) = st->dereverb_enabled;
      break;

   case SPEEX_PREPROCESS_SET_DEREVERB_LEVEL:
      st->reverb_level = (*(float*)ptr);
      break;
   case SPEEX_PREPROCESS_GET_DEREVERB_LEVEL:
      (*(float*)ptr) = st->reverb_level;
      break;
   
   case SPEEX_PREPROCESS_SET_DEREVERB_DECAY:
      st->reverb_decay = (*(float*)ptr);
      break;
   case SPEEX_PREPROCESS_GET_DEREVERB_DECAY:
      (*(float*)ptr) = st->reverb_decay;
      break;

   case SPEEX_PREPROCESS_SET_PROB_START:
      st->speech_prob_start = (*(int*)ptr) / 100.0;
      if ( st->speech_prob_start > 1 || st->speech_prob_start < 0 )
         st->speech_prob_start = SPEEX_PROB_START_DEFAULT;
      break;
   case SPEEX_PREPROCESS_GET_PROB_START:
      (*(int*)ptr) = st->speech_prob_start * 100;
      break;

   case SPEEX_PREPROCESS_SET_PROB_CONTINUE:
      st->speech_prob_continue = (*(int*)ptr) / 100.0;
      if ( st->speech_prob_continue > 1 || st->speech_prob_continue < 0 )
         st->speech_prob_continue = SPEEX_PROB_CONTINUE_DEFAULT;
      break;
   case SPEEX_PREPROCESS_GET_PROB_CONTINUE:
      (*(int*)ptr) = st->speech_prob_continue * 100;
      break;

      default:
      speex_warning_int("Unknown speex_preprocess_ctl request: ", request);
      return -1;
   }
   return 0;
}
