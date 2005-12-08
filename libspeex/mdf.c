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

#define MDF_C

#include "misc.h"
//#include "speex/speex_echo.h"
#include "smallft.h"
#include "fftwrap.h"

#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#undef BETA
#define BETA .65

#define min(a,b) ((a)<(b) ? (a) : (b))
#define max(a,b) ((a)>(b) ? (a) : (b))

#ifdef FIXED_POINT
#define WEIGHT_SHIFT 11
#define WEIGHT_SCALING 2048
#define WEIGHT_SCALING_1 0.00048828f
//#define WEIGHT_SCALING (100*16*128.f)
//#define WEIGHT_SCALING_1 (.0625f*0.0078125f)
#else
#define WEIGHT_SCALING 1.f
#define WEIGHT_SCALING_1 1.f
#endif

/** Speex echo cancellation state. */
typedef struct {
   int frame_size;           /**< Number of samples processed each time */
   int window_size;
   int M;
   int cancel_count;
   int adapted;
   float sum_adapt;
   spx_word16_t *e;
   spx_word16_t *x;
   spx_word16_t *X;
   spx_word16_t *d;
   spx_word16_t *y;
   spx_word16_t *last_y;
   spx_word32_t *Yps;
   spx_word16_t *Y;
   spx_word16_t *E;
   spx_word16_t *PHI;
   spx_word16_t *W;
   spx_word32_t *power;
   float *power_1;
   spx_word32_t *Rf;
   spx_word32_t *Yf;
   spx_word32_t *Xf;
   spx_word32_t *Eh;
   spx_word32_t *Yh;
   float Pey;
   float Pyy;
   /*struct drft_lookup *fft_lookup;*/
   void *fft_table;


} SpeexEchoState;


/** Compute inner product of two real vectors */
static inline float inner_prod(spx_word16_t *x, spx_word16_t *y, int N)
{
   int i;
   float ret=0;
   for (i=0;i<N;i++)
      ret += (1.f*x[i])*y[i];
   return ret;
}

/** Compute power spectrum of a half-complex (packed) vector */
static inline void power_spectrum(spx_word16_t *X, spx_word32_t *ps, int N)
{
   int i, j;
   ps[0]=MULT16_16(X[0],X[0]);
   for (i=1,j=1;i<N-1;i+=2,j++)
   {
      ps[j] =  MULT16_16(X[i],X[i]) + MULT16_16(X[i+1],X[i+1]);
   }
   ps[j]=MULT16_16(X[i],X[i]);
}

/** Compute cross-power spectrum of a half-complex (packed) vectors and add to acc */
#ifdef FIXED_POINT
static inline void spectral_mul_accum(spx_word16_t *X, spx_word16_t *Y, spx_word16_t *acc, int N, int M)
{
   int i,j;
   spx_word32_t tmp1=0,tmp2=0;
   for (j=0;j<M;j++)
   {
      tmp1 = MAC16_16(tmp1, X[j*N],Y[j*N]);
   }
   acc[0] = SHR32(tmp1,WEIGHT_SHIFT);
   for (i=1;i<N-1;i+=2)
   {
      tmp1 = tmp2 = 0;
      for (j=0;j<M;j++)
      {
         tmp1 = SUB32(MAC16_16(tmp1, X[j*N+i],Y[j*N+i]), MULT16_16(X[j*N+i+1],Y[j*N+i+1]));
         tmp2 = MAC16_16(MAC16_16(tmp2, X[j*N+i+1],Y[j*N+i]), X[j*N+i], Y[j*N+i+1]);
      }
      acc[i] = SHR32(tmp1,WEIGHT_SHIFT);
      acc[i+1] = SHR32(tmp2,WEIGHT_SHIFT);
   }
   tmp1 = tmp2 = 0;
   for (j=0;j<M;j++)
   {
      tmp1 = MAC16_16(tmp1, X[(j+1)*N-1],Y[(j+1)*N-1]);
   }
   acc[N-1] = SHR32(tmp1,WEIGHT_SHIFT);
}
#else
static inline void spectral_mul_accum(spx_word16_t *X, spx_word16_t *Y, spx_word16_t *acc, int N, int M)
{
   int i,j;
   for (i=0;i<N;i++)
      acc[i] = 0;
   for (j=0;j<M;j++)
   {
      acc[0] += X[0]*Y[0];
      for (i=1;i<N-1;i+=2)
      {
         acc[i] += (X[i]*Y[i] - X[i+1]*Y[i+1]);
         acc[i+1] += (X[i+1]*Y[i] + X[i]*Y[i+1]);
      }
      acc[i] += X[i]*Y[i];
      X += N;
      Y += N;
   }
}
#endif

/** Compute weighted cross-power spectrum of a half-complex (packed) vector with conjugate */
static inline void weighted_spectral_mul_conj(float *w, spx_word16_t *X, spx_word16_t *Y, spx_word16_t *prod, int N)
{
   int i, j;
   prod[0] = w[0]*MULT16_16(X[0],Y[0]);
   for (i=1,j=1;i<N-1;i+=2,j++)
   {
      prod[i] = w[j]*MAC16_16(MULT16_16(X[i],Y[i]), X[i+1],Y[i+1]);
      prod[i+1] = w[j]*MAC16_16(MULT16_16(-X[i+1],Y[i]), X[i],Y[i+1]);
   }
   prod[i] = w[j]*MULT16_16(X[i],Y[i]);
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

   st->fft_table = spx_fft_init(N);
   
   st->e = (spx_word16_t*)speex_alloc(N*sizeof(spx_word16_t));
   st->x = (spx_word16_t*)speex_alloc(N*sizeof(spx_word16_t));
   st->d = (spx_word16_t*)speex_alloc(N*sizeof(spx_word16_t));
   st->y = (spx_word16_t*)speex_alloc(N*sizeof(spx_word16_t));
   st->Yps = (spx_word32_t*)speex_alloc(N*sizeof(spx_word32_t));
   st->last_y = (spx_word16_t*)speex_alloc(N*sizeof(spx_word16_t));
   st->Yf = (spx_word32_t*)speex_alloc((st->frame_size+1)*sizeof(spx_word32_t));
   st->Rf = (spx_word32_t*)speex_alloc((st->frame_size+1)*sizeof(spx_word32_t));
   st->Xf = (spx_word32_t*)speex_alloc((st->frame_size+1)*sizeof(spx_word32_t));
   st->Yh = (spx_word32_t*)speex_alloc((st->frame_size+1)*sizeof(spx_word32_t));
   st->Eh = (spx_word32_t*)speex_alloc((st->frame_size+1)*sizeof(spx_word32_t));

   st->X = (spx_word16_t*)speex_alloc(M*N*sizeof(spx_word16_t));
   st->Y = (spx_word16_t*)speex_alloc(N*sizeof(spx_word16_t));
   st->E = (spx_word16_t*)speex_alloc(N*sizeof(spx_word16_t));
   st->W = (spx_word16_t*)speex_alloc(M*N*sizeof(spx_word16_t));
   st->PHI = (spx_word16_t*)speex_alloc(M*N*sizeof(spx_word16_t));
   st->power = (spx_word32_t*)speex_alloc((frame_size+1)*sizeof(spx_word32_t));
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
   spx_fft_destroy(st->fft_table);

   speex_free(st->e);
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

extern int fixed_point;
/** Performs echo cancellation on a frame */
void speex_echo_cancel(SpeexEchoState *st, short *ref, short *echo, short *out, float *Yout)
{
   int i,j;
   int N,M;
   float Syy=0,See=0;
   float leak_estimate;
   float ss;
   float adapt_rate;
   
   N = st->window_size;
   M = st->M;
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

   /* Convert x (echo input) to frequency domain */
   spx_fft(st->fft_table, st->x, &st->X[(M-1)*N]);
   
   /* Compute filter response Y */
   spectral_mul_accum(st->X, st->W, st->Y, N, M);
   
   spx_ifft(st->fft_table, st->Y, st->y);

   /* Compute error signal (signal with echo removed) */ 
   for (i=0;i<st->frame_size;i++)
   {
      spx_word32_t tmp_out;
      tmp_out = SUB32(EXTEND32(ref[i]), EXTEND32(st->y[i+st->frame_size]));
      
      st->e[i] = 0;
      /* Do we need saturation? */
      st->e[i+st->frame_size] = tmp_out;
      
      /* Saturation */
      if (tmp_out>32767)
         tmp_out = 32767;
      else if (tmp_out<-32768)
         tmp_out = -32768;
      out[i] = tmp_out;
   }
   
   /* Compute a bunch of correlations */
   See = inner_prod(st->e+st->frame_size, st->e+st->frame_size, st->frame_size);
   Syy = inner_prod(st->y+st->frame_size, st->y+st->frame_size, st->frame_size);
   
   /* Convert error to frequency domain */
   spx_fft(st->fft_table, st->e, st->E);
   for (i=0;i<st->frame_size;i++)
      st->y[i] = 0;
   spx_fft(st->fft_table, st->y, st->Y);

   /* Compute power spectrum of echo (X), error (E) and filter response (Y) */
   power_spectrum(st->E, st->Rf, N);
   power_spectrum(st->Y, st->Yf, N);
   power_spectrum(&st->X[(M-1)*N], st->Xf, N);
   
   /* Smooth echo energy estimate over time */
   for (j=0;j<=st->frame_size;j++)
      st->power[j] = (1-ss)*st->power[j] + 1 + ss*st->Xf[j];

   {
      float Pey = 0, Pyy=0;
      float alpha;
      for (j=0;j<=st->frame_size;j++)
      {
         spx_word32_t E, Y, Eh, Yh;
         E = (st->Rf[j]);
         Y = (st->Yf[j]);
         Eh = st->Eh[j] + E;
         Yh = st->Yh[j] + Y;
         Pey += Eh*1.f*Yh;
         Pyy += Yh*1.f*Yh;
         st->Eh[j] = .95*Eh - E;
         st->Yh[j] = .95*Yh - Y;
      }
      alpha = .02*Syy / (1e4+See);
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
         r = leak_estimate*st->Yf[i] / (1.f+st->Rf[i]);
         if (r>.5)
            r = .5;
         st->power_1[i] = WEIGHT_SCALING*adapt_rate*r/(1.f+st->power[i]);
         /*printf ("%f ", st->power_1[i]);*/
      }
   } else {
      for (i=0;i<=st->frame_size;i++)
         st->power_1[i] = WEIGHT_SCALING*adapt_rate/(1.f+st->power[i]);
   }
   /* Compute weight gradient */
   for (j=0;j<M;j++)
   {
      weighted_spectral_mul_conj(st->power_1, &st->X[j*N], st->E, st->PHI+N*j, N);
   }

   /* Gradient descent */
   for (i=0;i<M*N;i++)
   {
      st->W[i] += st->PHI[i];
      //printf ("%f ", st->PHI[i]);
   }
   /*if (st->cancel_count==1100)
      for (i=0;i<M*N;i++)
   printf ("%f ", st->W[i]);*/
   /*printf ("\n");*/
   
   /* AUMDF weight constraint */
   for (j=0;j<M;j++)
   {
      /* Remove the "if" to make this an MDF filter */
      //if(1)
      if (j==M-1 || st->cancel_count%(M-1) == j)
      {
         float w[N];
#ifdef FIXED_POINT
         float w2[N];
         for (i=0;i<N;i++)
            w2[i] = .03125*st->W[j*N+i];
         spx_ifft_float(st->fft_table, w2, w);
         for (i=0;i<st->frame_size;i++)
         {
            w[i]=0;
         }
         for (i=st->frame_size;i<N;i++)
         {
            w[i]*=4;
         }
         spx_fft_float(st->fft_table, w, w2);
         for (i=0;i<N;i++)
         {
            w2[i]*=.25;
         }
         for (i=0;i<N;i++)
            st->W[j*N+i] -= 32*w2[i];
#else
         float w2[N];
         fixed_point = 0;
         spx_ifft_float(st->fft_table, &st->W[j*N], w);
         for (i=st->frame_size;i<N;i++)
         {
            w[i]=0;
         }
         spx_fft_float(st->fft_table, w, &st->W[j*N]);
         fixed_point=1;
#endif
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
         st->y[i] = (.5-.5*cos(2*M_PI*i/N))*st->last_y[i];
      
      /* Compute power spectrum of the echo */
      spx_fft(st->fft_table, st->y, st->Y);
      power_spectrum(st->Y, st->Yps, N);
      
      /* Estimate residual echo */
      for (i=0;i<=st->frame_size;i++)
         Yout[i] = 2.f*leak_estimate*st->Yps[i];
   }

}

