/* Copyright (C) Jean-Marc Valin

   File: speex_echo.c


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

#include "misc.h"
#include "speex_echo.h"
#include "smallft.h"
#include <math.h>

/** Creates a new echo canceller state */
SpeexEchoState *speex_echo_state_init(int frame_size, int filter_length)
{
   int N,M;
   SpeexEchoState *st = (SpeexEchoState *)speex_alloc(sizeof(SpeexEchoState));

   st->frame_size = frame_size;
   st->window_size = 2*frame_size;
   N = st->window_size;
   M = st->M = (filter_length+N-1)/frame_size;
   st->cancel_count=0;

   drft_init(&st->fft_lookup, N);
   
   st->x = (float*)speex_alloc(N*sizeof(float));
   st->y = (float*)speex_alloc(N*sizeof(float));

   st->X = (float*)speex_alloc(M*N*sizeof(float));
   st->Y = (float*)speex_alloc(N*sizeof(float));
   st->E = (float*)speex_alloc(N*sizeof(float));
   st->W = (float*)speex_alloc(M*N*sizeof(float));
   st->PHI = (float*)speex_alloc(N*sizeof(float));
   st->power = (float*)speex_alloc((frame_size+1)*sizeof(float));
   st->power_1 = (float*)speex_alloc((frame_size+1)*sizeof(float));
   return st;
}

/** Destroys an echo canceller state */
void speex_echo_state_destroy(SpeexEchoState *st)
{
}

/** Performs echo cancellation a frame */
void speex_echo_cancel(SpeexEchoState *st, float *ref, float *echo, float *out, float *Yout)
{
   int i,j,m;
   int N,M;
   float scale;

   N = st->window_size;
   M = st->M;
   scale = 1.0/N;
   st->cancel_count++;

   for (i=0;i<st->frame_size;i++)
   {
      st->x[i] = st->x[i+st->frame_size];
      st->x[i+st->frame_size] = echo[i];
   }

   /* Shift memory: this could be optimized eventually*/
   for (i=0;i<N*(M-1);i++)
      st->X[i]=st->X[i+N];

   for (i=0;i<N;i++)
      st->X[(M-1)*N+i]=st->x[i];

   drft_forward(&st->fft_lookup, &st->X[(M-1)*N]);

   /* Compute filter response Y */
   for (i=1;i<N-1;i+=2)
   {
      st->Y[i] = st->Y[i+1] = 0;
      for (j=0;j<M;j++)
      {
         st->Y[i] += st->X[j*N+i]*st->W[j*N+i] - st->X[j*N+i+1]*st->W[j*N+i+1];
         st->Y[i+1] += st->X[j*N+i+1]*st->W[j*N+i] + st->X[j*N+i]*st->W[j*N+i+1];
      }
   }
   st->Y[0] = st->Y[N-1] = 0;
   for (j=0;j<M;j++)
   {
      st->Y[0] += st->X[j*N]*st->W[j*N];
      st->Y[N-1] += st->X[(j+1)*N-1]*st->W[(j+1)*N-1];
   }

   if (Yout)
      for (i=0;i<N;i++)
         Yout[i] = st->Y[i];

   for (i=0;i<N;i++)
      st->y[i] = st->Y[i];
   
   /* Filter response in time domain */
   drft_backward(&st->fft_lookup, st->y);
   for (i=0;i<N;i++)
      st->y[i] *= scale;

   /* Compute error signal (echo canceller output) */
   for (i=0;i<st->frame_size;i++)
   {
      out[i] = ref[i] - st->y[i+st->frame_size];
      st->E[i] = 0;
      st->E[i+st->frame_size] = out[i];
   }

   drft_forward(&st->fft_lookup, st->E);


   /* Compute input power in each band */
   {
      float s;
      float tmp, tmp2;

      if (st->cancel_count<M)
         s = 1.0/st->cancel_count;
      else
         s = 1.0/M;
      
      for (i=1,j=1;i<N-1;i+=2,j++)
      {
         tmp=0;
         for (m=0;m<M;m++)
         {
            tmp += st->X[m*N+i]*st->X[m*N+i] + st->X[m*N+i+1]*st->X[m*N+i+1];
         }
         tmp *= s;
         if (st->cancel_count<M)
            st->power[j] = tmp;
         else
            st->power[j] = .8*st->power[j] + .2*tmp;
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
         st->power[0] = .8*st->power[0] + .2*tmp;
         st->power[st->frame_size] = .8*st->power[st->frame_size] + .2*tmp2;
      }
      
      for (i=0;i<=st->frame_size;i++)
         st->power_1[i] = 1.0/(1e-10+st->power[i]);
   }

   /* Update filter weights */
   for (j=0;j<M;j++)
   {
      for (i=1;i<N-1;i+=2)
      {
         st->PHI[i] = st->X[j*N+i]*st->E[i] + st->X[j*N+i+1]*st->E[i+1];
         st->PHI[i+1] = -st->X[j*N+i+1]*st->E[i] + st->X[j*N+i]*st->E[i+1];
      }
      st->PHI[0] = st->X[j*N]*st->E[0];
      st->PHI[N-1] = st->X[(j+1)*N-1]*st->E[N-1];

      /* Optionally perform some transform (normalization?) on PHI */
      
      for (i=1,m=1;i<N-1;i+=2,m++)
         st->W[j*N+i] += .1*st->PHI[i]*st->power_1[m];
      st->W[j*N] += .1*st->PHI[0]*st->power_1[0];
      st->W[(j+1)*N-1] += .1*st->PHI[N-1]*st->power_1[st->frame_size];
      /*
      for (i=0,j=0;i<=N;i++,j++)
         st->W[j*N+i] += .001*st->PHI[i];
      */
   }
}

