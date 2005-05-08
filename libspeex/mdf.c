/* Copyright (C) Jean-Marc Valin

   File: speex_echo.c
   Echo cancelling based on the MDF algorithm described in:

   J. S. Soo, K. K. Pang Multidelay block frequency adaptive filter, 
   IEEE Trans. Acoust. Speech Signal Process., Vol. ASSP-38, No. 2, 
   February 1990.

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

#include "misc.h"
#include <speex/speex_echo.h>
#include "smallft.h"
#include <math.h>
/*#include <stdio.h>*/

#define BETA .65
//#define BETA 0

#define min(a,b) ((a)<(b) ? (a) : (b))
#define max(a,b) ((a)>(b) ? (a) : (b))

static inline float inner_prod(float *x, float *y, int N)
{
   int i;
   float ret=0;
   for (i=0;i<N;i++)
      ret += x[i]*y[i];
   return ret;
}

static inline void power_spectrum(float *X, float *ps, int N)
{
   int i, j;
   ps[0]=X[0]*X[0];
   for (i=1,j=1;i<N-1;i+=2,j++)
   {
      ps[j] =  X[i]*X[i] + X[i+1]*X[i+1];
   }
   ps[j]=X[i]*X[i];
}

static inline void spectral_mul_accum(float *X, float *Y, float *acc, int N)
{
   int i;
   acc[0] += X[0]*Y[0];
   for (i=1;i<N-1;i+=2)
   {
      acc[i] += (X[i]*Y[i] - X[i+1]*Y[i+1]);
      acc[i+1] += (X[i+1]*Y[i] + X[i]*Y[i+1]);
   }
   acc[i] += X[i]*Y[i];
}

static inline void spectral_mul_conj(float *X, float *Y, float *prod, int N)
{
   int i;
   prod[0] = X[0]*Y[0];
   for (i=1;i<N-1;i+=2)
   {
      prod[i] = (X[i]*Y[i] + X[i+1]*Y[i+1]);
      prod[i+1] = (-X[i+1]*Y[i] + X[i]*Y[i+1]);
   }
   prod[i] = X[i]*Y[i];
}


static inline void weighted_spectral_mul_conj(float *w, float *X, float *Y, float *prod, int N)
{
   int i, j;
   prod[0] = w[0]*X[0]*Y[0];
   for (i=1,j=1;i<N-1;i+=2,j++)
   {
      prod[i] = w[j]*(X[i]*Y[i] + X[i+1]*Y[i+1]);
      prod[i+1] = w[j]*(-X[i+1]*Y[i] + X[i]*Y[i+1]);
   }
   prod[i] = w[j]*X[i]*Y[i];
}


/** Creates a new echo canceller state */
SpeexEchoState *speex_echo_state_init(int frame_size, int filter_length)
{
   int i,N,M;
   SpeexEchoState *st = (SpeexEchoState *)speex_alloc(sizeof(SpeexEchoState));

   st->frame_size = frame_size;
   st->window_size = 2*frame_size;
   N = st->window_size;
   M = st->M = (filter_length+st->frame_size-1)/frame_size;
   st->cancel_count=0;
   st->adapt_rate = .01f;

   st->fft_lookup = (struct drft_lookup*)speex_alloc(sizeof(struct drft_lookup));
   spx_drft_init(st->fft_lookup, N);
   
   st->x = (float*)speex_alloc(N*sizeof(float));
   st->d = (float*)speex_alloc(N*sizeof(float));
   st->y = (float*)speex_alloc(N*sizeof(float));
   st->Yps = (float*)speex_alloc(N*sizeof(float));
   st->last_y = (float*)speex_alloc(N*sizeof(float));
   st->Yf = (float*)speex_alloc((st->frame_size+1)*sizeof(float));
   st->Rf = (float*)speex_alloc((st->frame_size+1)*sizeof(float));
   st->Xf = (float*)speex_alloc((st->frame_size+1)*sizeof(float));
   st->fratio = (float*)speex_alloc((st->frame_size+1)*sizeof(float));

   st->X = (float*)speex_alloc(M*N*sizeof(float));
   st->D = (float*)speex_alloc(N*sizeof(float));
   st->Y = (float*)speex_alloc(N*sizeof(float));
   st->E = (float*)speex_alloc(N*sizeof(float));
   st->W = (float*)speex_alloc(M*N*sizeof(float));
   st->PHI = (float*)speex_alloc(N*sizeof(float));
   st->power = (float*)speex_alloc((frame_size+1)*sizeof(float));
   st->power_1 = (float*)speex_alloc((frame_size+1)*sizeof(float));
   st->grad = (float*)speex_alloc(N*M*sizeof(float));
   
   for (i=0;i<N*M;i++)
   {
      st->W[i] = 0;
   }
   
   st->adapted = 0;
   return st;
}

void speex_echo_reset(SpeexEchoState *st)
{
   int i, M, N;
   st->cancel_count=0;
   st->adapt_rate = .01f;
   N = st->window_size;
   M = st->M;
   for (i=0;i<N*M;i++)
   {
      st->W[i] = 0;
      st->X[i] = 0;
   }
}

/** Destroys an echo canceller state */
void speex_echo_state_destroy(SpeexEchoState *st)
{
   spx_drft_clear(st->fft_lookup);
   speex_free(st->fft_lookup);
   speex_free(st->x);
   speex_free(st->d);
   speex_free(st->y);
   speex_free(st->last_y);
   speex_free(st->Yps);
   speex_free(st->Yf);
   speex_free(st->Rf);
   speex_free(st->Xf);
   speex_free(st->fratio);

   speex_free(st->X);
   speex_free(st->D);
   speex_free(st->Y);
   speex_free(st->E);
   speex_free(st->W);
   speex_free(st->PHI);
   speex_free(st->power);
   speex_free(st->power_1);
   speex_free(st->grad);

   speex_free(st);
}

/** Performs echo cancellation on a frame */
void speex_echo_cancel(SpeexEchoState *st, short *ref, short *echo, short *out, float *Yout)
{
   int i,j;
   int N,M;
   float scale;
   float ESR;
   float SER;
   float Sry=0,Srr=0,Syy=0,Sey=0,See=0,Sxx=0;

   N = st->window_size;
   M = st->M;
   scale = 1.0f/N;
   st->cancel_count++;

   /* Copy input data to buffer */
   for (i=0;i<st->frame_size;i++)
   {
      st->x[i] = st->x[i+st->frame_size];
      st->x[i+st->frame_size] = echo[i];

      st->d[i] = st->d[i+st->frame_size];
      st->d[i+st->frame_size] = ref[i];
   }

   /* Shift memory: this could be optimized eventually*/
   for (i=0;i<N*(M-1);i++)
      st->X[i]=st->X[i+N];

   /* Copy new echo frame */
   for (i=0;i<N;i++)
      st->X[(M-1)*N+i]=st->x[i];

   /* Convert x (echo input) to frequency domain */
   spx_drft_forward(st->fft_lookup, &st->X[(M-1)*N]);

   /* Compute filter response Y */
   for (i=0;i<N;i++)
      st->Y[i] = 0;
   for (j=0;j<M;j++)
      spectral_mul_accum(&st->X[j*N], &st->W[j*N], st->Y, N);
   
   /* Convert Y (filter response) to time domain */
   for (i=0;i<N;i++)
      st->y[i] = st->Y[i];
   spx_drft_backward(st->fft_lookup, st->y);
   for (i=0;i<N;i++)
      st->y[i] *= scale;

   /* Transform d (reference signal) to frequency domain */
   for (i=0;i<N;i++)
      st->D[i]=st->d[i];
   spx_drft_forward(st->fft_lookup, st->D);

   /* Compute error signal (signal with echo removed) */ 
   for (i=0;i<st->frame_size;i++)
   {
      float tmp_out;
      tmp_out = (float)ref[i] - st->y[i+st->frame_size];
      
      st->E[i] = 0;
      st->E[i+st->frame_size] = tmp_out;
      
      /* Saturation */
      if (tmp_out>32767)
         tmp_out = 32767;
      else if (tmp_out<-32768)
         tmp_out = -32768;
      out[i] = tmp_out;  
   }
   
   /* Compute power spectrum of output (D-Y) and filter response (Y) */
   for (i=0;i<N;i++)
      st->D[i] -= st->Y[i];
   power_spectrum(st->D, st->Rf, N);
   power_spectrum(st->Y, st->Yf, N);
   
   /* Compute frequency-domain adaptation mask */
   for (j=0;j<=st->frame_size;j++)
   {
      float r;
      r = .1*st->Yf[j] / (1+st->Rf[j]);
      if (r>1)
         r = 1;
      st->fratio[j] = r;
      //printf ("%f ", r);
   }
   //printf ("\n");

   /* Compute a bunch of correlations */
   Sry = inner_prod(st->y+st->frame_size, st->d+st->frame_size, st->frame_size);
   Sey = inner_prod(st->y+st->frame_size, st->E+st->frame_size, st->frame_size);
   See = inner_prod(st->E+st->frame_size, st->E+st->frame_size, st->frame_size);
   Syy = inner_prod(st->y+st->frame_size, st->y+st->frame_size, st->frame_size);
   Srr = inner_prod(st->d+st->frame_size, st->d+st->frame_size, st->frame_size);
   Sxx = inner_prod(st->x+st->frame_size, st->x+st->frame_size, st->frame_size);
      
   SER = Srr / (1+Sxx);
   ESR = .1*Syy / (1+See);
   if (ESR>1)
      ESR = 1;
   
   if (Sey/(1+Syy) < -.09 && ESR > .3)
   {
      for (i=0;i<M*N;i++)
         st->W[i] *= max(0.8,1+Sey/(1+Syy)-.05);
   }
   if (Sey/(1+Syy) > .2 && (ESR > .3 || SER < 1))
   {
      for (i=0;i<M*N;i++)
         st->W[i] *= 1.05;
   }
   
   if (ESR>.6 && st->cancel_count > 20)
   //if (st->cancel_count > 40)
   {
      if (!st->adapted)
         fprintf(stderr, "Adapted at %d\n", st->cancel_count); 
      st->adapted = 1;
   }
   //printf ("%f %f %f %f %f %f %f %f\n", Srr, Syy, Sxx, See, ESR, SER, Sry, Sey);
   for (i=0;i<=st->frame_size;i++)
   {
      st->fratio[i]  = (.2*ESR+.8*min(.005+ESR,st->fratio[i]));
      printf ("%f ", st->fratio[i]);
   }
   printf ("\n");
   
   
   if (st->adapted)
   {
      st->adapt_rate = .95f/(2+M);
   } else {
      if (SER<.1)
         st->adapt_rate =.8/(2+M);
      else if (SER<1)
         st->adapt_rate =.4/(2+M);
      else if (SER<10)
         st->adapt_rate =.2/(2+M);
      else if (SER<30)
         st->adapt_rate =.08/(2+M);
      else
         st->adapt_rate = 0;
   }
   

   /* Compute input power in each frequency bin */
   {
#if 0
      float s;
      float tmp, tmp2;
      int m;
      if (st->cancel_count<M)
         s = 1.0f/st->cancel_count;
      else
         s = 1.0f/M;
      
      for (i=1,j=1;i<N-1;i+=2,j++)
      {
#if 0
         tmp=0;
         for (m=0;m<M;m++)
         {
            float E = st->X[m*N+i]*st->X[m*N+i] + st->X[m*N+i+1]*st->X[m*N+i+1];
            tmp += E;
            if (st->power[j] < .02*E)
               st->power[j] = .02*E;
            
         }
         tmp *= s;
         if (st->cancel_count<M)
            st->power[j] = tmp;
         else
            st->power[j] = BETA*st->power[j] + (1-BETA)*tmp;
#else
         float E;
         float ss = 1.0f/st->cancel_count;
         if (ss < .3/M)
            ss=.3/M;
         E = (st->X[(M-1)*N+i]*st->X[(M-1)*N+i] + st->X[(M-1)*N+i+1]*st->X[(M-1)*N+i+1]);
         st->power[j] = (1-ss)*st->power[j] + ss*E;
#endif
      }
      tmp=tmp2=0;
      for (m=0;m<M;m++)
      {
         tmp += st->X[m*N]*st->X[m*N];
         tmp2 += st->X[(m+1)*N-1]*st->X[(m+1)*N-1];
      }
      tmp *= s;
      tmp2 *= s;
      if (st->cancel_count<M)
      {
         st->power[0] = tmp;
         st->power[st->frame_size] = tmp2;
      } else {
         st->power[0] = BETA*st->power[0] + (1-BETA)*tmp;
         st->power[st->frame_size] = BETA*st->power[st->frame_size] + (1-BETA)*tmp2;
      }
#else
      
      float ss = 1.0f/st->cancel_count;
      if (ss < .3/M)
         ss=.3/M;
      power_spectrum(&st->X[(M-1)*N], st->Xf, N);
      for (j=0;j<=st->frame_size;j++)
         st->power[j] = (1-ss)*st->power[j] + ss*st->Xf[j];
#endif      
      if (st->adapted)
      {
         for (i=0;i<=st->frame_size;i++)
         {
            st->power_1[i] = st->fratio[i] /(1+st->power[i]);
         }
      } else {
         for (i=0;i<=st->frame_size;i++)
            st->power_1[i] = 1.0f/(1e5f+st->power[i]);
      }
   }

   
   /* Convert error to frequency domain */
   spx_drft_forward(st->fft_lookup, st->E);

   /* Compute weight gradient */
   for (j=0;j<M;j++)
   {
      weighted_spectral_mul_conj(st->power_1, &st->X[j*N], st->E, st->PHI, N);

      for (i=0;i<N;i++)
      {
         st->grad[j*N+i] = st->PHI[i];
      }
   }
   
   /* Update weights */
   for (i=0;i<M*N;i++)
      st->W[i] += st->adapt_rate*st->grad[i];

   /* AUMDF weight constraint */
   for (j=0;j<M;j++)
   {
      /* Remove the "if" to make this an MDF filter */
      if (st->cancel_count%M == j)
      {
         spx_drft_backward(st->fft_lookup, &st->W[j*N]);
         for (i=0;i<N;i++)
            st->W[j*N+i]*=scale;
         for (i=st->frame_size;i<N;i++)
         {
            st->W[j*N+i]=0;
         }
         spx_drft_forward(st->fft_lookup, &st->W[j*N]);
      }
   }


   /* Compute spectrum of estimated echo for use in an echo post-filter (if necessary)*/
   if (Yout)
   {
      if (st->adapted)
      {
         for (i=0;i<st->frame_size;i++)
            st->last_y[i] = st->last_y[st->frame_size+i];
         for (i=0;i<st->frame_size;i++)
            st->last_y[st->frame_size+i] = st->y[st->frame_size+i];
      } else {
         for (i=0;i<N;i++)
            st->last_y[i] = st->x[i];
      }
      for (i=0;i<N;i++)
         st->Yps[i] = (.5-.5*cos(2*M_PI*i/N))*st->last_y[i];
      spx_drft_forward(st->fft_lookup, st->Yps);
#if 0
      for (i=1;i<st->frame_size;i++)
         st->Yps[i] = .1*st->Yps[2*i-1]*st->Yps[2*i-1] + st->Yps[2*i]*st->Yps[2*i];
      st->Yps[0] = .1*st->Yps[0]*st->Yps[0];
      st->Yps[st->frame_size] = .1*st->Yps[N-1]*st->Yps[N-1];
#else
      power_spectrum(st->Yps, st->Yps, N);
#endif      
      for (i=0;i<=st->frame_size;i++)
         Yout[i] = .3*st->Yps[i];
   }

}

