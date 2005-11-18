/* Copyright (C) 2003-2005 Jean-Marc Valin

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
#include "speex/speex_echo.h"
#include "smallft.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#undef BETA
#define BETA .65

#define min(a,b) ((a)<(b) ? (a) : (b))
#define max(a,b) ((a)>(b) ? (a) : (b))

/** Compute inner product of two real vectors */
static inline float inner_prod(float *x, float *y, int N)
{
   int i;
   float ret=0;
   for (i=0;i<N;i++)
      ret += x[i]*y[i];
   return ret;
}

/** Compute power spectrum of a half-complex (packed) vector */
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

/** Compute cross-power spectrum of a half-complex (packed) vectors and add to acc */
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

/** Compute weighted cross-power spectrum of a half-complex (packed) vector with conjugate */
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
   st->sum_adapt = 0;
         
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
   st->Yh = (float*)speex_alloc((st->frame_size+1)*sizeof(float));
   st->Eh = (float*)speex_alloc((st->frame_size+1)*sizeof(float));

   st->X = (float*)speex_alloc(M*N*sizeof(float));
   st->Y = (float*)speex_alloc(N*sizeof(float));
   st->E = (float*)speex_alloc(N*sizeof(float));
   st->W = (float*)speex_alloc(M*N*sizeof(float));
   st->PHI = (float*)speex_alloc(M*N*sizeof(float));
   st->power = (float*)speex_alloc((frame_size+1)*sizeof(float));
   st->power_1 = (float*)speex_alloc((frame_size+1)*sizeof(float));
   
   for (i=0;i<N*M;i++)
   {
      st->W[i] = st->PHI[i] = 0;
   }

   st->adapted = 0;
   st->Pey = st->Pyy = 0;
   return st;
}

/** Resets echo canceller state */
void speex_echo_state_reset(SpeexEchoState *st)
{
   int i, M, N;
   st->cancel_count=0;
   N = st->window_size;
   M = st->M;
   for (i=0;i<N*M;i++)
   {
      st->W[i] = 0;
      st->X[i] = 0;
   }
   for (i=0;i<=st->frame_size;i++)
      st->power[i] = 0;
   
   st->adapted = 0;
   st->sum_adapt = 0;
   st->Pey = st->Pyy = 0;

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
   speex_free(st->Yh);
   speex_free(st->Eh);

   speex_free(st->X);
   speex_free(st->Y);
   speex_free(st->E);
   speex_free(st->W);
   speex_free(st->PHI);
   speex_free(st->power);
   speex_free(st->power_1);

   speex_free(st);
}


/** Performs echo cancellation on a frame */
void speex_echo_cancel(SpeexEchoState *st, short *ref, short *echo, short *out, float *Yout)
{
   int i,j;
   int N,M;
   float scale;
   float Syy=0,See=0;
   float leak_estimate;
   float ss;
   float adapt_rate;

   N = st->window_size;
   M = st->M;
   scale = 1.0f/N;
   st->cancel_count++;
   ss = 1.0f/st->cancel_count;
   if (ss < .4/M)
      ss=.4/M;

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
   
   /* Compute a bunch of correlations */
   See = inner_prod(st->E+st->frame_size, st->E+st->frame_size, st->frame_size);
   Syy = inner_prod(st->y+st->frame_size, st->y+st->frame_size, st->frame_size);
   
   /* Convert error to frequency domain */
   spx_drft_forward(st->fft_lookup, st->E);
   for (i=0;i<st->frame_size;i++)
      st->y[i] = 0;
   for (i=0;i<N;i++)
      st->Y[i] = st->y[i];
   spx_drft_forward(st->fft_lookup, st->Y);
   
   /* Compute power spectrum of echo (X), error (E) and filter response (Y) */
   power_spectrum(st->E, st->Rf, N);
   power_spectrum(st->Y, st->Yf, N);
   power_spectrum(&st->X[(M-1)*N], st->Xf, N);
   
   /* Smooth echo energy estimate over time */
   for (j=0;j<=st->frame_size;j++)
      st->power[j] = (1-ss)*st->power[j] + ss*st->Xf[j];

   {
      float Pey = 0, Pyy=0;
      float alpha;
      for (j=0;j<=st->frame_size;j++)
      {
         float E, Y, Eh, Yh;
         E = (st->Rf[j]);
         Y = (st->Yf[j]);
         Eh = st->Eh[j] + E;
         Yh = st->Yh[j] + Y;
         Pey += Eh*Yh;
         Pyy += Yh*Yh;
         st->Eh[j] = .95*Eh - E;
         st->Yh[j] = .95*Yh - Y;
      }
      alpha = .02*Syy / (1+See);
      if (alpha > .02)
         alpha = .02;
      st->Pey = (1-alpha)*st->Pey + alpha*Pey;
      st->Pyy = (1-alpha)*st->Pyy + alpha*Pyy;
      if (st->Pey< .001*st->Pyy)
         st->Pey = .001*st->Pyy;
      leak_estimate = st->Pey / (1+st->Pyy);
      if (leak_estimate > 1)
         leak_estimate = 1;
      /*printf ("%f\n", leak_estimate);*/
   }
   
   if (!st->adapted)
   {
      float Sxx;
      Sxx = inner_prod(st->x+st->frame_size, st->x+st->frame_size, st->frame_size);

      /* We consider that the filter is adapted if the following is true*/
      if (st->sum_adapt > 1)
         st->adapted = 1;

      /* Temporary adaption rate if filter is not adapted correctly */
      adapt_rate = .2f * Sxx / (1e4+See);
      if (adapt_rate>.2)
         adapt_rate = .2;
      adapt_rate /= M;
      
      /* How much have we adapted so far? */
      st->sum_adapt += adapt_rate;
   }

   if (st->adapted)
   {
      adapt_rate = 1.f/M;
      for (i=0;i<=st->frame_size;i++)
      {
         float r;
         /* Compute frequency-domain adaptation mask */
         r = leak_estimate*st->Yf[i] / (1+st->Rf[i]);
         if (r>1)
            r = 1;
         st->power_1[i] = adapt_rate*r/(1.f+st->power[i]);
      }
   } else {
      for (i=0;i<=st->frame_size;i++)
         st->power_1[i] = adapt_rate/(1.f+st->power[i]);      
   }

   /* Compute weight gradient */
   for (j=0;j<M;j++)
   {
      weighted_spectral_mul_conj(st->power_1, &st->X[j*N], st->E, st->PHI+N*j, N);
   }

   /* Gradient descent */
   for (i=0;i<M*N;i++)
      st->W[i] += st->PHI[i];
   
   /* AUMDF weight constraint */
   for (j=0;j<M;j++)
   {
      /* Remove the "if" to make this an MDF filter */
      if (j==M-1 || st->cancel_count%(M-1) == j)
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
         /* If the filter is adapted, take the filtered echo */
         for (i=0;i<st->frame_size;i++)
            st->last_y[i] = st->last_y[st->frame_size+i];
         for (i=0;i<st->frame_size;i++)
            st->last_y[st->frame_size+i] = st->y[st->frame_size+i];
      } else {
         /* If filter isn't adapted yet, all we can do is take the echo signal directly */
         for (i=0;i<N;i++)
            st->last_y[i] = st->x[i];
      }
      
      /* Apply hanning window (should pre-compute it)*/
      for (i=0;i<N;i++)
         st->Yps[i] = (.5-.5*cos(2*M_PI*i/N))*st->last_y[i];
      
      /* Compute power spectrum of the echo */
      spx_drft_forward(st->fft_lookup, st->Yps);
      power_spectrum(st->Yps, st->Yps, N);
      
      /* Estimate residual echo */
      for (i=0;i<=st->frame_size;i++)
         Yout[i] = 2.f*leak_estimate*st->Yps[i];
   }

}

