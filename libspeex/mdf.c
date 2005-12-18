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
#include "speex/speex_echo.h"
#include "smallft.h"
#include "fftwrap.h"
#include "pseudofloat.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define min(a,b) ((a)<(b) ? (a) : (b))
#define max(a,b) ((a)>(b) ? (a) : (b))

#ifdef FIXED_POINT
#define WEIGHT_SHIFT 11
#define WEIGHT_SCALING 2048
#else
#define WEIGHT_SCALING 1.f
#define WEIGHT_SHIFT 0
#endif

#ifdef FIXED_POINT
static const spx_float_t MAX_ALPHA = ((spx_float_t){16777, -21});
static const spx_float_t ALPHA0 = ((spx_float_t){26214, -19});
static const spx_float_t MIN_LEAK = ((spx_float_t){16777, -24});
static const spx_float_t SPEC_AVERAGE = ((spx_float_t){20972, -20});
static const spx_float_t SPEC_AVERAGE_1 = ((spx_float_t){32113,-15});
#else
static const spx_float_t MAX_ALPHA = .008f;
static const spx_float_t ALPHA0 = .05f;
static const spx_float_t MIN_LEAK = .001f;
static const spx_float_t SPEC_AVERAGE = .02f;
static const spx_float_t SPEC_AVERAGE_1 = .98f;
#endif

/** Speex echo cancellation state. */
struct SpeexEchoState_ {
   int frame_size;           /**< Number of samples processed each time */
   int window_size;
   int M;
   int cancel_count;
   int adapted;
   spx_word32_t sum_adapt;
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
   spx_float_t *power_1;
   spx_word32_t *Rf;
   spx_word32_t *Yf;
   spx_word32_t *Xf;
   spx_word32_t *Eh;
   spx_word32_t *Yh;
   spx_float_t Pey;
   spx_float_t Pyy;
   spx_word16_t *window;
   /*struct drft_lookup *fft_lookup;*/
   void *fft_table;
   spx_word16_t memX, memD, memE;
   spx_word16_t preemph;
};

static inline spx_word32_t inner_prod(const spx_word16_t *x, const spx_word16_t *y, int len)
{
   spx_word32_t sum=0;
   len >>= 2;
   while(len--)
   {
      spx_word32_t part=0;
      part = MAC16_16(part,*x++,*y++);
      part = MAC16_16(part,*x++,*y++);
      part = MAC16_16(part,*x++,*y++);
      part = MAC16_16(part,*x++,*y++);
      /* HINT: If you had a 40-bit accumulator, you could shift only at the end */
      sum = ADD32(sum,SHR32(part,6));
   }
   return sum;
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
   acc[0] = PSHR32(tmp1,WEIGHT_SHIFT);
   for (i=1;i<N-1;i+=2)
   {
      tmp1 = tmp2 = 0;
      for (j=0;j<M;j++)
      {
         tmp1 = SUB32(MAC16_16(tmp1, X[j*N+i],Y[j*N+i]), MULT16_16(X[j*N+i+1],Y[j*N+i+1]));
         tmp2 = MAC16_16(MAC16_16(tmp2, X[j*N+i+1],Y[j*N+i]), X[j*N+i], Y[j*N+i+1]);
      }
      acc[i] = PSHR32(tmp1,WEIGHT_SHIFT);
      acc[i+1] = PSHR32(tmp2,WEIGHT_SHIFT);
   }
   tmp1 = tmp2 = 0;
   for (j=0;j<M;j++)
   {
      tmp1 = MAC16_16(tmp1, X[(j+1)*N-1],Y[(j+1)*N-1]);
   }
   acc[N-1] = PSHR32(tmp1,WEIGHT_SHIFT);
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
static inline void weighted_spectral_mul_conj(spx_float_t *w, spx_word16_t *X, spx_word16_t *Y, spx_word16_t *prod, int N)
{
   int i, j;
   prod[0] = FLOAT_MUL32(w[0],MULT16_16(X[0],Y[0]));
   for (i=1,j=1;i<N-1;i+=2,j++)
   {
      prod[i] = FLOAT_MUL32(w[j],MAC16_16(MULT16_16(X[i],Y[i]), X[i+1],Y[i+1]));
      prod[i+1] = FLOAT_MUL32(w[j],MAC16_16(MULT16_16(-X[i+1],Y[i]), X[i],Y[i+1]));
   }
   prod[i] = FLOAT_MUL32(w[j],MULT16_16(X[i],Y[i]));
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
   st->power_1 = (spx_float_t*)speex_alloc((frame_size+1)*sizeof(spx_float_t));
   st->window = (spx_word16_t*)speex_alloc(N*sizeof(spx_word16_t));
#ifdef FIXED_POINT   
   for (i=0;i<N;i++)
      st->window[i] = 32767.*(.5-.5*cos(2*M_PI*i/N));
#else
   for (i=0;i<N;i++)
      st->window[i] = .5-.5*cos(2*M_PI*i/N);
#endif
   for (i=0;i<N*M;i++)
   {
      st->W[i] = st->PHI[i] = 0;
   }
   st->memX=st->memD=st->memE=0;
   st->preemph = QCONST16(.9,15);
   st->adapted = 0;
   st->Pey = st->Pyy = FLOAT_ONE;
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
   st->Pey = st->Pyy = FLOAT_ONE;

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
   speex_free(st->window);

   speex_free(st);
}

extern int fixed_point;
/** Performs echo cancellation on a frame */
void speex_echo_cancel(SpeexEchoState *st, short *ref, short *echo, short *out, float *Yout)
{
   int i,j;
   int N,M;
   spx_word32_t Syy,See;
   spx_word16_t leak_estimate;
   spx_word16_t ss, ss_1;
   spx_float_t Pey = FLOAT_ONE, Pyy=FLOAT_ONE;
   spx_float_t alpha;
   float RER;
   spx_word32_t tmp32;
   spx_word16_t M_1;
   
   N = st->window_size;
   M = st->M;
   st->cancel_count++;
#ifdef FIXED_POINT
   ss=DIV32_16(11469,M);
   ss_1 = SUB16(32767,ss);
   M_1 = DIV32_16(32767,M);
#else
   ss=.35/M;
   ss_1 = 1-ss;
   M_1 = 1.f/M;
#endif

   /* Copy input data to buffer */
   for (i=0;i<st->frame_size;i++)
   {
      st->x[i] = st->x[i+st->frame_size];
      st->x[i+st->frame_size] = SHL16(SUB16(echo[i], MULT16_16_P15(st->preemph, st->memX)),1);
      st->memX = echo[i];
      
      st->d[i] = st->d[i+st->frame_size];
      st->d[i+st->frame_size] = SHL16(SUB16(ref[i], MULT16_16_P15(st->preemph, st->memD)),1);
      st->memD = ref[i];
   }

   /* Shift memory: this could be optimized eventually*/
   for (i=0;i<N*(M-1);i++)
      st->X[i]=st->X[i+N];

   /* Convert x (echo input) to frequency domain */
   spx_fft(st->fft_table, st->x, &st->X[(M-1)*N]);
   
   /* Compute filter response Y */
   spectral_mul_accum(st->X, st->W, st->Y, N, M);
   
   spx_ifft(st->fft_table, st->Y, st->y);

#if 1
   spectral_mul_accum(st->X, st->PHI, st->Y, N, M);   
   spx_ifft(st->fft_table, st->Y, st->e);
#endif

   /* Compute error signal (for the output with de-emphasis) */ 
   for (i=0;i<st->frame_size;i++)
   {
      spx_word32_t tmp_out;
#if 1
      spx_word16_t y = MULT16_16_Q15(st->window[i+st->frame_size],st->e[i+st->frame_size]) + MULT16_16_Q15(st->window[i],st->y[i+st->frame_size]);
      tmp_out = SUB32(EXTEND32(st->d[i+st->frame_size]), EXTEND32(y));
#else
      tmp_out = SUB32(EXTEND32(st->d[i+st->frame_size]), EXTEND32(st->y[i+st->frame_size]));
#endif

      tmp_out = PSHR32(tmp_out,1);
      /* Saturation */
      if (tmp_out>32767)
         tmp_out = 32767;
      else if (tmp_out<-32768)
         tmp_out = -32768;
      tmp_out = ADD32(tmp_out, EXTEND32(MULT16_16_P15(st->preemph, st->memE)));
      out[i] = tmp_out;
      st->memE = tmp_out;
   }

   /* Compute error signal (filter update version) */ 
   for (i=0;i<st->frame_size;i++)
   {
      st->e[i] = 0;
      st->e[i+st->frame_size] = st->d[i+st->frame_size] - st->y[i+st->frame_size];
   }

   /* Compute a bunch of correlations */
   See = inner_prod(st->e+st->frame_size, st->e+st->frame_size, st->frame_size);
   See = ADD32(See, SHR32(10000,6));
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
      st->power[j] = MULT16_32_Q15(ss_1,st->power[j]) + 1 + MULT16_32_Q15(ss,st->Xf[j]);

   /* Compute filtered spectra and (cross-)correlations */
   for (j=st->frame_size;j>=0;j--)
   {
      spx_float_t Eh, Yh;
      Eh = PSEUDOFLOAT(st->Rf[j] - st->Eh[j]);
      Yh = PSEUDOFLOAT(st->Yf[j] - st->Yh[j]);
      Pey = FLOAT_ADD(Pey,FLOAT_MULT(Eh,Yh));
      Pyy = FLOAT_ADD(Pyy,FLOAT_MULT(Yh,Yh));
#ifdef FIXED_POINT
      st->Eh[j] = MAC16_32_Q15(MULT16_32_Q15(32113,st->Eh[j]), 655, st->Rf[j]);
      st->Yh[j] = MAC16_32_Q15(MULT16_32_Q15(32113,st->Yh[j]), 655, st->Yf[j]);
#else
      st->Eh[j] = .98*st->Eh[j] + .02*st->Rf[j];
      st->Yh[j] = .98*st->Yh[j] + .02*st->Yf[j];
#endif
   }
   
   /* Compute correlation updatete rate */
   tmp32 = MULT16_32_Q15(QCONST16(.05,15),Syy);
   if (tmp32 > MULT16_32_Q15(QCONST16(.008,15),See))
      tmp32 = MULT16_32_Q15(QCONST16(.008,15),See);
   alpha = FLOAT_DIV32(tmp32, See);
   spx_float_t alpha_1 = FLOAT_SUB(FLOAT_ONE, alpha);
   /* Update correlations (recursive average) */
   st->Pey = FLOAT_ADD(FLOAT_MULT(alpha_1,st->Pey) , FLOAT_MULT(alpha,Pey));
   st->Pyy = FLOAT_ADD(FLOAT_MULT(alpha_1,st->Pyy) , FLOAT_MULT(alpha,Pyy));
   if (FLOAT_LT(st->Pyy, FLOAT_ONE))
      st->Pyy = FLOAT_ONE;
   /* We don't really hope to get better than 33 dB (MIN_LEAK-3dB) attenuation anyway */
   if (FLOAT_LT(st->Pey, FLOAT_MULT(MIN_LEAK,st->Pyy)))
      st->Pey = FLOAT_MULT(MIN_LEAK,st->Pyy);
   if (FLOAT_GT(st->Pey, st->Pyy))
      st->Pey = st->Pyy;
   /* leak_estimate is the limear regression result */
   leak_estimate = FLOAT_EXTRACT16(FLOAT_SHL(FLOAT_DIVU(st->Pey, st->Pyy),14));
   if (leak_estimate > 16383)
      leak_estimate = 32767;
   else
      leak_estimate = SHL16(leak_estimate,1);

   /*printf ("%f\n", leak_estimate);*/
   
   RER = 3.*MULT16_32_Q15(leak_estimate,Syy) / See;
   if (RER > .5)
      RER = .5;
   
   /* We consider that the filter has had minimal adaptation if the following is true*/
   if (!st->adapted && st->sum_adapt > QCONST32(1,15))
   {
      st->adapted = 1;
   }

   if (st->adapted)
   {
      for (i=0;i<=st->frame_size;i++)
      {
         spx_word32_t r, e;
         /* Compute frequency-domain adaptation mask */
         r = MULT16_32_Q15(leak_estimate,SHL32(st->Yf[i],3));
         e = SHL32(st->Rf[i],3)+1;
#ifdef FIXED_POINT
         if (r>SHR32(e,1))
            r = SHR32(e,1);
#else
         if (r>.5*e)
            r = .5*e;
#endif
         r = MULT16_32_Q15(QCONST16(.8,15),r) + MULT16_32_Q15(QCONST16(.2,15),(spx_word32_t)(RER*e));
         /*st->power_1[i] = adapt_rate*r/(e*(1+st->power[i]));*/
         st->power_1[i] = FLOAT_SHL(FLOAT_DIV32_FLOAT(MULT16_32_Q15(M_1,r),FLOAT_MUL32U(e,st->power[i]+10)),WEIGHT_SHIFT);
      }
   } else {
      spx_word32_t Sxx;
      spx_word16_t adapt_rate=0;

      Sxx = inner_prod(st->x+st->frame_size, st->x+st->frame_size, st->frame_size);
      /* Temporary adaption rate if filter is not adapted correctly */

      tmp32 = MULT16_32_Q15(QCONST16(.15f, 15), Sxx);
#ifdef FIXED_POINT
      if (Sxx > SHR32(See,2))
         Sxx = SHR32(See,2);
#else
      if (Sxx > .25See)
         Sxx = .25*See;      
#endif
      adapt_rate = FLOAT_EXTRACT16(FLOAT_SHL(FLOAT_DIV32(MULT16_32_Q15(M_1,Sxx), See),15));
      
      for (i=0;i<=st->frame_size;i++)
         st->power_1[i] = FLOAT_SHL(FLOAT_DIV32(EXTEND32(adapt_rate),ADD32(st->power[i],10)),WEIGHT_SHIFT-15);


      /* How much have we adapted so far? */
      st->sum_adapt = ADD32(st->sum_adapt,adapt_rate);
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
      /* Old value of W in PHI */
      st->PHI[i] = st->W[i] - st->PHI[i];
   }
   
   /* AUMDF weight constraint */
   for (j=0;j<M;j++)
   {
      /* Remove the "if" to make this an MDF filter */
      if (j==M-1 || st->cancel_count%(M-1) == j)
      {
         spx_word16_t w[N];
#ifdef FIXED_POINT
         spx_word16_t w2[N];
         for (i=0;i<N;i++)
            w2[i] = PSHR16(st->W[j*N+i],5);
         spx_ifft(st->fft_table, w2, w);
         for (i=0;i<st->frame_size;i++)
         {
            w[i]=0;
         }
         for (i=st->frame_size;i<N;i++)
         {
            w[i]=SHL(w[i],2);
         }
         spx_fft(st->fft_table, w, w2);
         for (i=0;i<N;i++)
         {
            w2[i]=PSHR(w2[i],4);
         }
         for (i=0;i<N;i++)
            st->W[j*N+i] -= SHL16(w2[i],5);
#else
         spx_ifft(st->fft_table, &st->W[j*N], w);
         for (i=st->frame_size;i<N;i++)
         {
            w[i]=0;
         }
         spx_fft(st->fft_table, w, &st->W[j*N]);
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
            st->last_y[st->frame_size+i] = ref[i]-out[i];
      } else {
         /* If filter isn't adapted yet, all we can do is take the echo signal directly */
         for (i=0;i<N;i++)
            st->last_y[i] = st->x[i];
      }
      
      /* Apply hanning window (should pre-compute it)*/
      for (i=0;i<N;i++)
         st->y[i] = MULT16_16_Q15(st->window[i],st->last_y[i]);
      
      /* Compute power spectrum of the echo */
      spx_fft(st->fft_table, st->y, st->Y);
      power_spectrum(st->Y, st->Yps, N);
      
      /* Estimate residual echo */
      for (i=0;i<=st->frame_size;i++)
         Yout[i] = N*N*max(.2,3.f*leak_estimate)*st->Yps[i];
   }

}

