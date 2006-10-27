/* Copyright (C) 2003 Epic Games (written by Jean-Marc Valin)
   Copyright (C) 2004-2006 Epic Games 
   
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


/*
   Recommended papers:
   
   Y. Ephraim and D. Malah, "Speech enhancement using minimum mean-square error
   short-time spectral amplitude estimator". IEEE Transactions on Acoustics, 
   Speech and Signal Processing, vol. ASSP-32, no. 6, pp. 1109-1121, 1984.
   
   Y. Ephraim and D. Malah, "Speech enhancement using minimum mean-square error
   log-spectral amplitude estimator". IEEE Transactions on Acoustics, Speech and 
   Signal Processing, vol. ASSP-33, no. 2, pp. 443-445, 1985.
   
   I. Cohen and B. Berdugo, "Speech enhancement for non-stationary noise environments".
   Signal Processing, vol. 81, no. 2, pp. 2403-2418, 2001.

   Stefan Gustafsson, Rainer Martin, Peter Jax, and Peter Vary. "A psychoacoustic 
   approach to combined acoustic echo cancellation and noise reduction". IEEE 
   Transactions on Speech and Audio Processing, 2002.
   
   J.-M. Valin, J. Rouat, and F. Michaud, "Microphone array post-filter for separation
   of simultaneous non-stationary sources". In Proceedings IEEE International 
   Conference on Acoustics, Speech, and Signal Processing, 2004.
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include "speex/speex_preprocess.h"
#include "speex/speex_echo.h"
#include "misc.h"
#include "fftwrap.h"
#include "filterbank.h"

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))

#ifndef M_PI
#define M_PI 3.14159263
#endif

#define LOUDNESS_EXP 2.5

#define NB_BANDS 24

#define SPEEX_PROB_START_DEFAULT       0.35f
#define SPEEX_PROB_CONTINUE_DEFAULT    0.20f
#define NOISE_SUPPRESS_DEFAULT       -25
#define ECHO_SUPPRESS_DEFAULT        -45
#define ECHO_SUPPRESS_ACTIVE_DEFAULT -15

#ifndef NULL
#define NULL 0
#endif

/** Speex pre-processor state. */
struct SpeexPreprocessState_ {
   /* Basic info */
   int    frame_size;        /**< Number of samples processed each time */
   int    ps_size;           /**< Number of points in the power spectrum */
   int    sampling_rate;     /**< Sampling rate of the input/output */
   int    nbands;
   FilterBank *bank;
   
   /* Parameters */
   int    denoise_enabled;
   int    agc_enabled;
   float  agc_level;
   int    vad_enabled;
   int    dereverb_enabled;
   float  reverb_decay;
   float  reverb_level;
   float  speech_prob_start;
   float  speech_prob_continue;
   int    noise_suppress;
   int    echo_suppress;
   int    echo_suppress_active;
   SpeexEchoState *echo_state;
   
   /* DSP-related arrays */
   spx_word16_t *frame;      /**< Processing frame (2*ps_size) */
   spx_word16_t *ft;         /**< Processing frame in freq domain (2*ps_size) */
   spx_word32_t *ps;         /**< Current power spectrum */
   float *gain2;             /**< Adjusted gains */
   float *gain_floor;        /**< Minimum gain allowed */
   spx_word16_t *window;     /**< Analysis/Synthesis window */
   spx_word32_t *noise;      /**< Noise estimate */
   spx_word32_t *reverb_estimate; /**< Estimate of reverb energy */
   spx_word32_t *old_ps;     /**< Power spectrum for last frame */
   float *gain;              /**< Ephraim Malah gain */
   float *prior;             /**< A-priori SNR */
   float *post;              /**< A-posteriori SNR */

   spx_word32_t *S;          /**< Smoothed power spectrum */
   spx_word32_t *Smin;       /**< See Cohen paper */
   spx_word32_t *Stmp;       /**< See Cohen paper */
   float *update_prob;       /**< Propability of speech presence for noise update */

   float *zeta;              /**< Smoothed a priori SNR */

   float *loudness_weight;   /**< Perceptual loudness curve */

   spx_word32_t *echo_noise;
   spx_word32_t *residual_echo;

   /* Misc */
   spx_word16_t *inbuf;      /**< Input buffer (overlapped analysis) */
   spx_word16_t *outbuf;     /**< Output buffer (for overlap and add) */

   int    was_speech;
   float  loudness;          /**< loudness estimate */
   float  loudness2;         /**< loudness estimate */
   int    nb_adapt;          /**< Number of frames used for adaptation so far */
   int    nb_loudness_adapt; /**< Number of frames used for loudness adaptation so far */
   int    min_count;         /**< Number of frames processed so far */
   void  *fft_lookup;        /**< Lookup table for the FFT */

};


static void conj_window(spx_word16_t *w, int len)
{
   int i;
   for (i=0;i<len;i++)
   {
      float tmp;
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
      tmp=(.5-.5*cos(x))*(.5-.5*cos(x));
      if (inv)
         tmp=1-tmp;
      w[i]=QCONST16(.999,15)*sqrt(tmp);
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
   return 1.f/(1.f+.15f/(x));
}

SpeexPreprocessState *speex_preprocess_state_init(int frame_size, int sampling_rate)
{
   int i;
   int N, N3, N4, M;

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
   st->reverb_decay = .0;
   st->reverb_level = .0;
   st->noise_suppress = NOISE_SUPPRESS_DEFAULT;
   st->echo_suppress = ECHO_SUPPRESS_DEFAULT;
   st->echo_suppress_active = ECHO_SUPPRESS_ACTIVE_DEFAULT;

   st->speech_prob_start = SPEEX_PROB_START_DEFAULT;
   st->speech_prob_continue = SPEEX_PROB_CONTINUE_DEFAULT;

   st->echo_state = NULL;
   
   st->nbands = NB_BANDS;
   M = st->nbands;
   st->bank = filterbank_new(M, sampling_rate, N, 1);
   
   st->frame = (spx_word16_t*)speex_alloc(2*N*sizeof(float));
   st->window = (spx_word16_t*)speex_alloc(2*N*sizeof(float));
   st->ft = (spx_word16_t*)speex_alloc(2*N*sizeof(float));
   
   st->ps = (spx_word32_t*)speex_alloc((N+M)*sizeof(float));
   st->noise = (spx_word32_t*)speex_alloc((N+M)*sizeof(float));
   st->echo_noise = (spx_word32_t*)speex_alloc((N+M)*sizeof(float));
   st->residual_echo = (spx_word32_t*)speex_alloc((N+M)*sizeof(float));
   st->reverb_estimate = (spx_word32_t*)speex_alloc((N+M)*sizeof(float));
   st->old_ps = (spx_word32_t*)speex_alloc((N+M)*sizeof(float));
   st->prior = (float*)speex_alloc((N+M)*sizeof(float));
   st->post = (float*)speex_alloc((N+M)*sizeof(float));
   st->gain = (float*)speex_alloc((N+M)*sizeof(float));
   st->gain2 = (float*)speex_alloc((N+M)*sizeof(float));
   st->gain_floor = (float*)speex_alloc((N+M)*sizeof(float));
   st->zeta = (float*)speex_alloc((N+M)*sizeof(float));
   
   st->S = (spx_word32_t*)speex_alloc(N*sizeof(float));
   st->Smin = (spx_word32_t*)speex_alloc(N*sizeof(float));
   st->Stmp = (spx_word32_t*)speex_alloc(N*sizeof(float));
   st->update_prob = (float*)speex_alloc(N*sizeof(float));
   
   st->loudness_weight = (float*)speex_alloc(N*sizeof(float));
   st->inbuf = (spx_word16_t*)speex_alloc(N3*sizeof(float));
   st->outbuf = (spx_word16_t*)speex_alloc(N3*sizeof(float));

   conj_window(st->window, 2*N3);
   for (i=2*N3;i<2*st->ps_size;i++)
      st->window[i]=QCONST16(.999,15);
   
   if (N4>0)
   {
      for (i=N3-1;i>=0;i--)
      {
         st->window[i+N3+N4]=st->window[i+N3];
         st->window[i+N3]=1;
      }
   }
   for (i=0;i<N+M;i++)
   {
      st->noise[i]=1.;
      st->reverb_estimate[i]=0.;
      st->old_ps[i]=1.;
      st->gain[i]=1;
      st->post[i]=1;
      st->prior[i]=1;
   }

   for (i=0;i<N;i++)
      st->update_prob[i] = 1;
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

   st->was_speech = 0;
   st->loudness = pow(6000,LOUDNESS_EXP);
   st->loudness2 = 6000;
   st->nb_loudness_adapt = 0;

   st->fft_lookup = spx_fft_init(2*N);

   st->nb_adapt=0;
   st->min_count=0;
   return st;
}

void speex_preprocess_state_destroy(SpeexPreprocessState *st)
{
   speex_free(st->frame);
   speex_free(st->ft);
   speex_free(st->ps);
   speex_free(st->gain2);
   speex_free(st->gain_floor);
   speex_free(st->window);
   speex_free(st->noise);
   speex_free(st->reverb_estimate);
   speex_free(st->old_ps);
   speex_free(st->gain);
   speex_free(st->prior);
   speex_free(st->post);
   speex_free(st->loudness_weight);
   speex_free(st->echo_noise);
   speex_free(st->residual_echo);

   speex_free(st->S);
   speex_free(st->Smin);
   speex_free(st->Stmp);
   speex_free(st->update_prob);
   speex_free(st->zeta);

   speex_free(st->inbuf);
   speex_free(st->outbuf);

   spx_fft_destroy(st->fft_lookup);
   filterbank_destroy(st->bank);
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
   spx_word32_t *ps=st->ps;

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
      st->frame[i] = MULT16_16_Q15(st->window[i], st->frame[i]);

   /* Perform FFT */
   spx_fft(st->fft_lookup, st->frame, st->ft);
         
   /* Power spectrum */
   ps[0]=1;
   for (i=1;i<N;i++)
      ps[i]=1+MULT16_16(st->ft[2*i-1],st->ft[2*i-1]) + MULT16_16(st->ft[2*i],st->ft[2*i]);

   filterbank_compute_bank32(st->bank, ps, ps+N);
}

static void update_noise_prob(SpeexPreprocessState *st)
{
   int i;
   int min_range;
   int N = st->ps_size;

   for (i=1;i<N-1;i++)
      st->S[i] = 100.f+ .8f*st->S[i] + .05f*st->ps[i-1]+.1f*st->ps[i]+.05f*st->ps[i+1];
   st->S[0] = 100.f+ .8f*st->S[0] + .2*st->ps[0];
   st->S[N-1] = 100.f+ .8f*st->S[N-1] + .2*st->ps[N-1];
   
   if (st->nb_adapt==1)
   {
      for (i=0;i<N;i++)
         st->Smin[i] = st->Stmp[i] = 0.f;
   }

   if (st->nb_adapt < 100)
      min_range = 15;
   else if (st->nb_adapt < 1000)
      min_range = 50;
   else if (st->nb_adapt < 10000)
      min_range = 150;
   else
      min_range = 300;
   if (st->min_count > min_range)
   {
      st->min_count = 0;
      for (i=0;i<N;i++)
      {
         st->Smin[i] = min(st->Stmp[i], st->S[i]);
         st->Stmp[i] = st->S[i];
      }
   } else {
      for (i=0;i<N;i++)
      {
         st->Smin[i] = min(st->Smin[i], st->S[i]);
         st->Stmp[i] = min(st->Stmp[i], st->S[i]);      
      }
   }
   for (i=0;i<N;i++)
   {
      st->update_prob[i] *= .2f;
      if (st->S[i] > 2.5*st->Smin[i])
         st->update_prob[i] += .8f;
      /*fprintf (stderr, "%f ", st->S[i]/st->Smin[i]);*/
      /*fprintf (stderr, "%f ", st->update_prob[i]);*/
   }

}

#define NOISE_OVERCOMPENS 1.

void speex_echo_get_residual(SpeexEchoState *st, spx_word32_t *Yout, int len);

int speex_preprocess(SpeexPreprocessState *st, spx_int16_t *x, spx_int32_t *echo)
{
   return speex_preprocess_run(st, x);
}

int speex_preprocess_run(SpeexPreprocessState *st, spx_int16_t *x)
{
   int i;
   int M;
   int N = st->ps_size;
   int N3 = 2*N - st->frame_size;
   int N4 = st->frame_size - N3;
   spx_word32_t *ps=st->ps;
   float Zframe=0, Pframe;
   float beta, beta_1;
   float echo_floor;
   float noise_floor;
   
   st->nb_adapt++;
   st->min_count++;
   
   beta =1.0f/st->nb_adapt;
   if (beta < .03f)
      beta=.03f;
   beta_1 = 1.0f-beta;
   M = st->nbands;
   /* Deal with residual echo if provided */
   if (st->echo_state)
   {
      speex_echo_get_residual(st->echo_state, st->residual_echo, N);
      for (i=0;i<N;i++)
         st->echo_noise[i] = MAX32(.6f*st->echo_noise[i], st->residual_echo[i]);
      filterbank_compute_bank32(st->bank, st->echo_noise, st->echo_noise+N);
   } else {
      for (i=0;i<N+M;i++)
         st->echo_noise[i] = 0;
   }
   preprocess_analysis(st, x);

   update_noise_prob(st);

   /* Noise estimation always updated for the 10 first frames */
   /*if (st->nb_adapt<10)
   {
      for (i=1;i<N-1;i++)
         st->update_prob[i] = 0;
   }
   */
   
   /* Update the noise estimate for the frequencies where it can be */
   for (i=0;i<N;i++)
   {
      if (st->update_prob[i]<.5f || st->ps[i] < st->noise[i])
         st->noise[i] = beta_1*st->noise[i] + beta*NOISE_OVERCOMPENS*st->ps[i];
   }
   filterbank_compute_bank32(st->bank, st->noise, st->noise+N);

   /* Special case for first frame */
   if (st->nb_adapt==1)
      for (i=0;i<N+M;i++)
         st->old_ps[i] = ps[i];

   /* Compute a posteriori SNR */
   for (i=0;i<N+M;i++)
   {
      float gamma = .1;
      spx_word32_t tot_noise = 1+st->noise[i] + st->echo_noise[i] + st->reverb_estimate[i];
      st->post[i] = 1.f*ps[i]/tot_noise - 1.f;
      if (st->post[i]>100.f)
         st->post[i]=100.f;
      /*gamma = .15+.85*st->prior[i]*st->prior[i]/((1+st->prior[i])*(1+st->prior[i]));*/
      gamma = .1+.9*(st->old_ps[i]/(1.f+st->old_ps[i]+tot_noise))*(st->old_ps[i]/(1.f+st->old_ps[i]+tot_noise));
      /* A priori SNR update */
      st->prior[i] = gamma*max(0.0f,st->post[i]) + (1.f-gamma)*st->old_ps[i]/tot_noise;
      if (st->prior[i]>100.f)
         st->prior[i]=100.f;
   }

   /*print_vec(st->post, N+M, "");*/

   /* Recursive average of the a priori SNR. A bit smoothed for the psd components */
   st->zeta[0] = .7f*st->zeta[0] + .3f*st->prior[0];
   for (i=1;i<N-1;i++)
      st->zeta[i] = .7f*st->zeta[i] + .3f*(.5f*st->prior[i]+.25f*st->prior[i+1]+.25f*st->prior[i-1]);
   for (i=N-1;i<N+M;i++)
      st->zeta[i] = .7f*st->zeta[i] + .3f*st->prior[i];

   /* Speech probability of presence for the entire frame is based on the average filterbank a priori SNR */
   Zframe = 0;
   for (i=N;i<N+M;i++)
      Zframe += st->zeta[i];
   Zframe /= st->nbands;
   Pframe = .1+.9*qcurve(Zframe);
   
   noise_floor = exp(.2302585f*st->noise_suppress);
   echo_floor = exp(.2302585f* (st->echo_suppress*(1-Pframe) + st->echo_suppress_active*Pframe));
   /*print_vec(&Pframe, 1, "");*/
   /*for (i=N;i<N+M;i++)
      st->gain2[i] = qcurve (st->zeta[i]);
   filterbank_compute_psd(st->bank,st->gain2+N, st->gain2);*/
   
   for (i=N;i<N+M;i++)
   {
      float theta, MM;
      float prior_ratio;
      float q;
      float P1;

      /* Compute the gain floor based on different floors for the background noise and residual echo */
      st->gain_floor[i] = sqrt((noise_floor*st->noise[i] + echo_floor*st->echo_noise[i])/(1+st->noise[i] + st->echo_noise[i]));
      prior_ratio = st->prior[i]/(1.f+st->prior[i]);
      theta = (1.f+st->post[i])*prior_ratio;

      MM = hypergeom_gain(theta);
      st->gain[i] = prior_ratio * MM;
      /* Bound on the gain */
      if (st->gain[i]>1.f)
         st->gain[i]=1.f;
      
      /* Save old Bark power spectrum */
      st->old_ps[i] = .2*st->old_ps[i] + .8*st->gain[i]*st->gain[i]*ps[i];

      P1 = .2+.8*qcurve (st->zeta[i]);
      q = 1-Pframe*P1;
      st->gain2[i]=1.f/(1.f + (q/(1.f-q))*(1.f+st->prior[i])*exp(-theta));
   }
   filterbank_compute_psd(st->bank,st->gain2+N, st->gain2);
   filterbank_compute_psd(st->bank,st->gain+N, st->gain);
   
   /* Use 1 for linear gain resolution (best) or 0 for Bark gain resolution (faster) */
   if (1)
   {
      filterbank_compute_psd(st->bank,st->gain_floor+N, st->gain_floor);
   
      /* Compute gain according to the Ephraim-Malah algorithm */
      for (i=0;i<N;i++)
      {
         float MM;
         float theta;
         float prior_ratio;
         float p;
         float g;
         
         /* Wiener filter gain */
         prior_ratio = st->prior[i]/(1.f+st->prior[i]);
         theta = (1.f+st->post[i])*prior_ratio;
         p = st->gain2[i];
         
         /* Optimal estimator for loudness domain */
         MM = hypergeom_gain(theta);
         g = prior_ratio * MM;
         
         /* Constrain the gain to be close to the Bark scale gain */
         if (g > 3*st->gain[i])
            g = 3*st->gain[i];
         st->gain[i] = g;
         
         /* Bound on the gain */
         if (st->gain[i]>1.f)
            st->gain[i]=1.f;
         
         /* Save old power spectrum */
         st->old_ps[i] = .2*st->old_ps[i] + .8*st->gain[i]*st->gain[i]*ps[i];
         
         /* Apply gain floor */
         if (st->gain[i] < st->gain_floor[i])
            st->gain[i] = st->gain_floor[i];
         
         /* Exponential decay model for reverberation (unused) */
         /*st->reverb_estimate[i] = st->reverb_decay*st->reverb_estimate[i] + st->reverb_decay*st->reverb_level*st->gain[i]*st->gain[i]*st->ps[i];*/
         
         /* Take into account speech probability of presence (loudness domain MMSE estimator) */
         st->gain2[i]=(p*sqrt(st->gain[i])+sqrt(st->gain_floor[i])*(1-p)) * (p*sqrt(st->gain[i])+sqrt(st->gain_floor[i])*(1-p));
         
         /* Use this if you want a log-domain MMSE estimator instead */
         /*st->gain2[i] = pow(st->gain[i], p) * pow(st->gain_floor[i],1.f-p);*/
         
      }
   } else {
      for (i=N;i<N+M;i++)
      {
         float p = st->gain2[i];
         if (st->gain[i] < st->gain_floor[i])
            st->gain[i] = st->gain_floor[i];
         st->gain2[i]=(p*sqrt(st->gain[i])+sqrt(st->gain_floor[i])*(1-p)) * (p*sqrt(st->gain[i])+sqrt(st->gain_floor[i])*(1-p));
      }
      filterbank_compute_psd(st->bank,st->gain2+N, st->gain2);
      
   }
   
   if (!st->denoise_enabled)
   {
      for (i=0;i<N+M;i++)
         st->gain2[i]=1.f;
   }
   
   if (st->agc_enabled)
      speex_compute_agc(st);

   /* Apply computed gain */
   for (i=1;i<N;i++)
   {
      st->ft[2*i-1] *= st->gain2[i];
      st->ft[2*i] *= st->gain2[i];
   }
   st->ft[0] *= st->gain2[0];
   st->ft[2*N-1] *= st->gain2[N-1];

   /* Inverse FFT with 1/N scaling */
   spx_ifft(st->fft_lookup, st->ft, st->frame);

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
      st->frame[i] = MULT16_16_Q15(st->window[i], st->frame[i]);

   /* Perform overlap and add */
   for (i=0;i<N3;i++)
      x[i] = st->outbuf[i] + st->frame[i];
   for (i=0;i<N4;i++)
      x[N3+i] = st->frame[N3+i];
   
   /* Update outbuf */
   for (i=0;i<N3;i++)
      st->outbuf[i] = st->frame[st->frame_size+i];

   if (st->vad_enabled)
   {
      if (Pframe > st->speech_prob_start || (st->was_speech && Pframe > st->speech_prob_continue))
      {
         st->was_speech=1;
         return 1;
      } else
      {
         st->was_speech=0;
         return 0;
      }
   } else {
      return 1;
   }
}

void speex_preprocess_estimate_update(SpeexPreprocessState *st, spx_int16_t *x, spx_int32_t *echo)
{
   int i;
   int N = st->ps_size;
   int N3 = 2*N - st->frame_size;
   int M;
   spx_word32_t *ps=st->ps;

   M = st->nbands;
   st->min_count++;
   
   preprocess_analysis(st, x);

   update_noise_prob(st);
   
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
      st->outbuf[i] = MULT16_16_Q15(x[st->frame_size-N3+i],st->window[st->frame_size+i]);

   /* Save old power spectrum */
   for (i=0;i<N+M;i++)
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
      st->vad_enabled = (*(spx_int32_t*)ptr);
      break;
   case SPEEX_PREPROCESS_GET_VAD:
      (*(spx_int32_t*)ptr) = st->vad_enabled;
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

   case SPEEX_PREPROCESS_SET_NOISE_SUPPRESS:
      st->noise_suppress = -ABS(*(spx_int32_t*)ptr);
      break;
   case SPEEX_PREPROCESS_GET_NOISE_SUPPRESS:
      (*(spx_int32_t*)ptr) = st->noise_suppress;
      break;
   case SPEEX_PREPROCESS_SET_ECHO_SUPPRESS:
      st->echo_suppress = -ABS(*(spx_int32_t*)ptr);
      break;
   case SPEEX_PREPROCESS_GET_ECHO_SUPPRESS:
      (*(spx_int32_t*)ptr) = st->echo_suppress;
      break;
   case SPEEX_PREPROCESS_SET_ECHO_SUPPRESS_ACTIVE:
      st->echo_suppress_active = -ABS(*(spx_int32_t*)ptr);
      break;
   case SPEEX_PREPROCESS_GET_ECHO_SUPPRESS_ACTIVE:
      (*(spx_int32_t*)ptr) = st->echo_suppress_active;
      break;
   case SPEEX_PREPROCESS_SET_ECHO_STATE:
      st->echo_state = (SpeexEchoState*)ptr;
      break;
   case SPEEX_PREPROCESS_GET_ECHO_STATE:
      ptr = (void*)st->echo_state;
      break;

   default:
      speex_warning_int("Unknown speex_preprocess_ctl request: ", request);
      return -1;
   }
   return 0;
}
