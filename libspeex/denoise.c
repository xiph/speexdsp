/* Copyright (C) 2003 Epic Games 
   Written by Jean-Marc Valin

   File: denoise.c


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

#include <math.h>
#include "speex_denoise.h"
#include <stdio.h>
#include "misc.h"

#define STABILITY_TIME 20
#define NB_LAST_PS 10

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))

#ifndef M_PI
#define M_PI 3.14159263
#endif

#define SQRT_M_PI_2 0.88623

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

DenoiseState *denoise_state_init(int frame_size)
{
   int i;
   int N, N3, N4;

   DenoiseState *st = (DenoiseState *)speex_alloc(sizeof(DenoiseState));
   st->frame_size = frame_size;

   /* Round ps_size down to the nearest power of two */
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

   N = st->ps_size;
   N3 = 2*N - st->frame_size;
   N4 = st->frame_size - N3;

   st->frame = (float*)speex_alloc(2*N*sizeof(float));
   st->ps = (float*)speex_alloc(N*sizeof(float));
   st->gain2 = (float*)speex_alloc(N*sizeof(float));
   st->window = (float*)speex_alloc(2*N*sizeof(float));
   st->noise = (float*)speex_alloc(N*sizeof(float));
   st->old_ps = (float*)speex_alloc(N*sizeof(float));
   st->gain = (float*)speex_alloc(N*sizeof(float));
   st->prior = (float*)speex_alloc(N*sizeof(float));
   st->post = (float*)speex_alloc(N*sizeof(float));
   st->min_ps = (float*)speex_alloc(N*sizeof(float));
   st->last_energy = (float*)speex_alloc(STABILITY_TIME*sizeof(float));
   st->last_ps = (float*)speex_alloc(NB_LAST_PS*N*sizeof(float));
   st->loudness_weight = (float*)speex_alloc(N*sizeof(float));
   st->inbuf = (float*)speex_alloc(N3*sizeof(float));
   st->outbuf = (float*)speex_alloc(N3*sizeof(float));

   conj_window(st->window, 2*N3);
   for (i=N3-1;i>=0;i--)
   {
      st->window[i+N3+N4]=st->window[i+N3];
      st->window[i+N3]=1;
   }

   for (i=0;i<N;i++)
   {
      st->noise[i]=1e4;
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
      float ff=((float)i)*128.0/4000.0;
      st->loudness_weight[i] = .35-.35*ff/16000+.73*exp(-.5*(ff-3800)*(ff-3800)/9e5);
      st->loudness_weight[i] *= st->loudness_weight[i];
   }

   drft_init(&st->fft_lookup,2*N);


   st->nb_adapt=0;
   st->consec_noise=0;
   st->nb_denoise=0;
   st->nb_min_estimate=0;
   st->last_update=0;
   st->last_id=0;
   return st;
}

void denoise_state_destroy(DenoiseState *st)
{
   speex_free(st->frame);
   speex_free(st->ps);
   speex_free(st->gain2);
   speex_free(st->window);
   speex_free(st->noise);
   speex_free(st->old_ps);
   speex_free(st->gain);
   speex_free(st->prior);
   speex_free(st->post);
   speex_free(st->min_ps);
   speex_free(st->last_energy);
   speex_free(st->last_ps);
   speex_free(st->loudness_weight);

   speex_free(st->inbuf);
   speex_free(st->outbuf);

   drft_clear(&st->fft_lookup);
   
   speex_free(st);
}

static void update_noise(DenoiseState *st, float *ps)
{
   int i;
   float beta;
   st->nb_adapt++;
   beta=1.0/st->nb_adapt;
   if (beta < .05)
      beta=.05;
   
   for (i=0;i<st->ps_size;i++)
      st->noise[i] = (1-beta)*st->noise[i] + beta*ps[i];   
}

int denoise(DenoiseState *st, float *x)
{
   int i;
   float mean_post=0;
   float mean_prior=0;
   float energy;
   int N = st->ps_size;
   int N3 = 2*N - st->frame_size;
   int N4 = st->frame_size - N3;
   float scale=.5/N;
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
   drft_forward(&st->fft_lookup, st->frame);

   /************************************************************** 
    *  Denoise in spectral domain using Ephraim-Malah algorithm  *
    **************************************************************/

   /* Power spectrum */
   ps[0]=1;
   for (i=1;i<N;i++)
      ps[i]=1+st->frame[2*i-1]*st->frame[2*i-1] + st->frame[2*i]*st->frame[2*i];

   energy=0;
   for (i=1;i<N;i++)
      energy += log(100+ps[i]);
   energy /= 160;
   st->last_energy[st->nb_denoise%STABILITY_TIME]=energy;

   if (st->nb_denoise>=STABILITY_TIME)
   {
      float E=0, E2=0;
      float std;
      for (i=0;i<STABILITY_TIME;i++)
      {
         E+=st->last_energy[i];
         E2+=st->last_energy[i]*st->last_energy[i];
      }
      E2=E2/STABILITY_TIME;
      E=E/STABILITY_TIME;
      std = sqrt(E2-E*E);
      if (std<.15 && st->last_update>20)
      {
         update_noise(st, &st->last_ps[st->last_id*N]);
      }
      /*fprintf (stderr, "%f\n", std);*/
   }

   st->nb_denoise++;
#if 0
   if (st->nb_min_estimate<50)
   {
      float ener=0;
      for (i=1;i<N;i++)
         ener += ps[i];
      /*fprintf (stderr, "%f\n", ener);*/
      if (ener < st->min_ener || st->nb_min_estimate==0)
      {
         st->min_ener = ener;
         for (i=1;i<N;i++)
            st->min_ps[i] = ps[i];
      }
      st->nb_min_estimate++;
   } else {
      float noise_ener=0;
      st->nb_min_estimate=0;
      for (i=1;i<N;i++)
         noise_ener += st->noise[i];
      /*fprintf (stderr, "%f %f\n", noise_ener, st->min_ener);*/
      if (0&&(st->last_update>50 && st->min_ener > 3*noise_ener) || st->last_update>50)
      {
         for (i=1;i<N;i++)
         {
            if (st->noise[i] < st->min_ps[i])
               st->noise[i] = st->min_ps[i];
         }
         /*fprintf (stderr, "tata %d\n",st->last_update);*/
         st->last_update=0;
      } else {
         /*fprintf (stderr, "+");*/
      }
   }
#endif

   /* Noise estimation always updated for the 20 first times */
   if (st->nb_adapt<20)
   {
      update_noise(st, ps);
      st->last_update=0;
   }

   /* Compute a posteriori SNR */
   for (i=1;i<N;i++)
   {
      st->post[i] = ps[i]/(1+st->noise[i]) - 1;
      if (st->post[i]>100)
         st->post[i]=100;
      if (st->post[i]<0)
        st->post[i]=0;
      mean_post+=st->post[i];
   }
   mean_post /= N;
   if (mean_post<0)
      mean_post=0;

   /* Special case for first frame */
   if (st->nb_adapt==1)
      for (i=1;i<N;i++)
         st->old_ps[i] = ps[i];

   /* Compute a priori SNR */
   {
      /* A priori update rate */
      float gamma;
      float min_gamma=0.05;
      gamma = 1.0/st->nb_denoise;

      /*Make update rate smaller when there's no speech*/
      if (mean_post<3)
         min_gamma *= (mean_post+.1);
      else
         min_gamma *= 3.1;

      if (gamma<min_gamma)
         gamma=min_gamma;

      for (i=1;i<N;i++)
      {
         
         /* A priori SNR update */
         st->prior[i] = gamma*max(0.0,st->post[i]) +
         (1-gamma)*st->gain[i]*st->gain[i]*st->old_ps[i]/st->noise[i];
         
         if (st->prior[i]>100)
            st->prior[i]=100;
         
         mean_prior+=st->prior[i];
      }
   }
   mean_prior /= N;

#if 0
   for (i=0;i<N;i++)
   {
      fprintf (stderr, "%f ", st->prior[i]);
   }
   fprintf (stderr, "\n");
#endif
   /*fprintf (stderr, "%f %f\n", mean_prior,mean_post);*/

   /* If SNR is low (both a priori and a posteriori), update the noise estimate*/
   if (mean_prior<.23 && mean_post < .5 && st->nb_adapt>=20)
   {
      st->consec_noise++;
   } else {
      st->consec_noise=0;
   }

   if (st->consec_noise>=3)
   {
      update_noise(st, st->old_ps);
      st->last_update=0;
   } else {
      st->last_update++;
   }

   /* Compute gain according to the Ephraim-Malah algorithm */
   for (i=1;i<N;i++)
   {
      float MM;
      float theta;
      float prior_ratio;

      prior_ratio = st->prior[i]/(1.0001+st->prior[i]);
      theta = (1+st->post[i])*prior_ratio;

      /* Approximation of:
         exp(-theta/2)*((1+theta)*I0(theta/2) + theta.*I1(theta/2))
         because I don't feel like computing Bessel functions
      */
      /*MM = -.22+1.155*sqrt(theta+1.1);*/
      MM=-.22+1.163*sqrt(theta+1.1)-.0015*theta;

      st->gain[i] = SQRT_M_PI_2*sqrt(prior_ratio/(1.0001+st->post[i]))*MM;
      if (st->gain[i]>1)
      {
         st->gain[i]=1;
      }
      /*st->gain[i] = prior_ratio;*/
   }
   st->gain[0]=0;
   st->gain[N-1]=0;

   for (i=1;i<N-1;i++)
   {
      st->gain2[i]=st->gain[i];
      if (st->gain2[i]<.1)
         st->gain2[i]=.1;
   }
   st->gain2[N-1]=0;

   {
      float loudness=0;
      for (i=2;i<N;i++)
      {
         loudness += scale*st->ps[i] * st->gain2[i] * st->gain2[i] * st->loudness_weight[i];
      }
      loudness=sqrt(loudness);
      fprintf (stderr, "%f\n", loudness);
   }

   /* Apply computed gain */
   for (i=1;i<N;i++)
   {
      st->frame[2*i-1] *= st->gain2[i];
      st->frame[2*i] *= st->gain2[i];
   }
   /* Get rid of the DC and very low frequencies */
   st->frame[0]=0;
   st->frame[1]=0;
   st->frame[2]=0;
   /* Nyquist frequency is mostly useless too */
   st->frame[2*N-1]=0;

   /* Inverse FFT with 1/N scaling */
   drft_backward(&st->fft_lookup, st->frame);

   for (i=0;i<2*N;i++)
      st->frame[i] *= scale*st->window[i];

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

   for (i=1;i<N;i++)
      st->last_ps[st->last_id*N+i] = ps[i];
   st->last_id++;
   if (st->last_id>=NB_LAST_PS)
      st->last_id=0;

   return 1;
}
