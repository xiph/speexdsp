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

#include "misc.h"
#include "speex_echo.h"
#include "smallft.h"
#include <math.h>
/*#include <stdio.h>*/

#define BETA .65

/** Creates a new echo canceller state */
SpeexEchoState *speex_echo_state_init(int frame_size, int filter_length)
{
   int i,N,M;
   SpeexEchoState *st = (SpeexEchoState *)speex_alloc(sizeof(SpeexEchoState));

   st->frame_size = frame_size;
   st->window_size = 2*frame_size;
   N = st->window_size;
   M = st->M = (filter_length+N-1)/frame_size;
   st->cancel_count=0;
   st->adapt_rate = .01f;

   st->fft_lookup = (struct drft_lookup*)speex_alloc(sizeof(struct drft_lookup));
   drft_init(st->fft_lookup, N);
   
   st->x = (float*)speex_alloc(N*sizeof(float));
   st->d = (float*)speex_alloc(N*sizeof(float));
   st->y = (float*)speex_alloc(N*sizeof(float));

   st->X = (float*)speex_alloc(M*N*sizeof(float));
   st->D = (float*)speex_alloc(N*sizeof(float));
   st->Y = (float*)speex_alloc(N*sizeof(float));
   st->E = (float*)speex_alloc(N*sizeof(float));
   st->W = (float*)speex_alloc(M*N*sizeof(float));
   st->PHI = (float*)speex_alloc(N*sizeof(float));
   st->power = (float*)speex_alloc((frame_size+1)*sizeof(float));
   st->power_1 = (float*)speex_alloc((frame_size+1)*sizeof(float));
   st->grad = (float*)speex_alloc(N*M*sizeof(float));
   st->old_grad = (float*)speex_alloc(N*M*sizeof(float));
   
   for (i=0;i<N*M;i++)
   {
      st->W[i] = 0;
   }
   return st;
}

/** Destroys an echo canceller state */
void speex_echo_state_destroy(SpeexEchoState *st)
{
   drft_clear(st->fft_lookup);
   speex_free(st->fft_lookup);
   speex_free(st->x);
   speex_free(st->d);
   speex_free(st->y);

   speex_free(st->X);
   speex_free(st->D);
   speex_free(st->Y);
   speex_free(st->E);
   speex_free(st->W);
   speex_free(st->PHI);
   speex_free(st->power);
   speex_free(st->power_1);
   speex_free(st->grad);
   speex_free(st->old_grad);

   speex_free(st);
}

/** Performs echo cancellation on a frame */
void speex_echo_cancel(SpeexEchoState *st, short *ref, short *echo, short *out, int *Yout)
{
   int i,j,m;
   int N,M;
   float scale;
   float powE=0, powY=0, powD=0;
   float spectral_dist=0;

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

   for (i=0;i<N;i++)
      st->X[(M-1)*N+i]=st->x[i];

   /* Convert x (echo input) to frequency domain */
   drft_forward(st->fft_lookup, &st->X[(M-1)*N]);


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


   /* Transform d (reference signal) to frequency domain */
   for (i=0;i<N;i++)
      st->D[i]=st->d[i];
   drft_forward(st->fft_lookup, st->D);

   /* Evaluate "spectral distance" between Y and D t odetect crosstalk */
   for (i=1;i<N-1;i+=2)
   {
      float Sdd, Syy, Sdy;
      Sdd = 1e4 + st->D[i]*st->D[i] + st->D[i+1]*st->D[i+1];
      Syy = 1e4 + st->Y[i]*st->Y[i] + st->Y[i+1]*st->Y[i+1];
      Sdy = st->Y[i]*st->D[i] + st->Y[i+1]*st->D[i+1];
      spectral_dist += Sdy/sqrt(Sdd*Syy);
   }
   spectral_dist *= 2*scale;
   /*printf ("%f\n", spectral_dist);*/

   /* Copy spectrum of Y to Yout for use in an echo post-filter */
   if (Yout)
   {
      for (i=1,j=1;i<N-1;i+=2,j++)
      {
         Yout[j] =  st->Y[i]*st->Y[i] + st->Y[i+1]*st->Y[i+1];
      }
      Yout[0] = Yout[st->frame_size] = 0;
      for (i=0;i<=st->frame_size;i++)
         Yout[i] *= .1;
   }

   for (i=0;i<N;i++)
      st->y[i] = st->Y[i];
   
   /* Convery Y (filter response) to time domain */
   drft_backward(st->fft_lookup, st->y);
   for (i=0;i<N;i++)
      st->y[i] *= scale;

   /* Compute error signal (echo canceller output) */
   for (i=0;i<st->frame_size;i++)
   {
      out[i] = ref[i] - st->y[i+st->frame_size];
      st->E[i] = 0;
      st->E[i+st->frame_size] = out[i];
   }

   /* Convert error to frequency domain */
   drft_forward(st->fft_lookup, st->E);

   for (i=0;i<st->frame_size;i++)
   {
      powD += N*((float)ref[i])*ref[i];
   }
#if 0
   for (i=1;i<N-1;i+=2)
   {
      float tmp;
      tmp = st->Y[i]*st->Y[i] + st->Y[i+1]*st->Y[i+1];
      powY += 1*tmp + 0*(st->E[i]*st->E[i] + st->E[i+1]*st->E[i+1]);
      tmp = st->E[i]*st->E[i] + st->E[i+1]*st->E[i+1] - 4*tmp;
      if (tmp<0)
         tmp=0;
      powE += tmp; 
   }
#else
   for (i=3;i<N-3;i+=2)
   {
      float tmp;
      tmp = .5*(st->Y[i]*st->Y[i] + st->Y[i+1]*st->Y[i+1]) + .25*
      (st->Y[i-2]*st->Y[i-2] + st->Y[i-1]*st->Y[i-1] 
       + st->Y[i+2]*st->Y[i+2] + st->Y[i+3]*st->Y[i+3]);
      powY += 1*tmp + 0*(st->E[i]*st->E[i] + st->E[i+1]*st->E[i+1]);
      tmp = st->E[i]*st->E[i] + st->E[i+1]*st->E[i+1] - .5*tmp;
      if (tmp<0)
         tmp=0;
      powE += tmp; 
   }
#endif

   /* Compute input power in each frequency bin */
   {
      float s;
      float tmp, tmp2;

      if (st->cancel_count<M)
         s = 1.0f/st->cancel_count;
      else
         s = 1.0f/M;
      
      for (i=1,j=1;i<N-1;i+=2,j++)
      {
         tmp=0;
         for (m=0;m<M;m++)
         {
            float E = st->X[m*N+i]*st->X[m*N+i] + st->X[m*N+i+1]*st->X[m*N+i+1];
            tmp += E;
            if (st->power[j] < .2*E)
               st->power[j] = .2*E;

         }
         tmp *= s;
         if (st->cancel_count<M)
            st->power[j] = tmp;
         else
            st->power[j] = BETA*st->power[j] + (1-BETA)*tmp;
      }
      tmp=tmp2=0;
      for (m=0;m<M;m++)
      {
         tmp += st->X[m*N]*st->X[m*N];
         tmp2 += st->X[(m+1)*N-1]*st->X[(m+1)*N-1];
         /*FIXME: Should put a bound on energy like several lines above */
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
      
      for (i=0;i<=st->frame_size;i++)
         st->power_1[i] = 1.0f/(1e6f+st->power[i]);
   }

   /* Compute weight gradient */
   for (j=0;j<M;j++)
   {
      for (i=1,m=1;i<N-1;i+=2,m++)
      {
         st->PHI[i] = st->power_1[m] 
         * (st->X[j*N+i]*st->E[i] + st->X[j*N+i+1]*st->E[i+1]);
         st->PHI[i+1] = st->power_1[m] 
         * (-st->X[j*N+i+1]*st->E[i] + st->X[j*N+i]*st->E[i+1]);
      }
      st->PHI[0] = st->power_1[0] * st->X[j*N]*st->E[0];
      st->PHI[N-1] = st->power_1[st->frame_size] * st->X[(j+1)*N-1]*st->E[N-1];
      

#if 0 /* Set to 1 to enable MDF instead of AUMDF (and comment out weight constraint below) */
      drft_backward(st->fft_lookup, st->PHI);
      for (i=0;i<N;i++)
         st->PHI[i]*=scale;
      for (i=st->frame_size;i<N;i++)
        st->PHI[i]=0;
      drft_forward(st->fft_lookup, st->PHI);
#endif
     

      for (i=0;i<N;i++)
      {
         st->old_grad[j*N+i] = st->PHI[i];
         st->grad[j*N+i] = st->PHI[i];
      }

      
   }

   /* Adjust adaptation rate */
   if (st->cancel_count>2*M)
   {
      if (st->cancel_count<8*M)
      {
         st->adapt_rate = .15f/(2+M);
      } else {
         st->adapt_rate = spectral_dist*(1.0f/(2+M));
         if (st->adapt_rate>.15f/(2+M))
            st->adapt_rate=.15f/(2+M);
         if (st->adapt_rate<0)
            st->adapt_rate=0;
      }
   } else
      st->adapt_rate = .0f;

   /* Update weights */
   for (i=0;i<M*N;i++)
      st->W[i] += st->adapt_rate*st->grad[i];

   /* AUMDF weight constraint */
   for (j=0;j<M;j++)
   {
      if (st->cancel_count%M == j)
      {
         drft_backward(st->fft_lookup, &st->W[j*N]);
         for (i=0;i<N;i++)
            st->W[j*N+i]*=scale;
         for (i=st->frame_size;i<N;i++)
         {
            st->W[j*N+i]=0;
         }
         drft_forward(st->fft_lookup, &st->W[j*N]);
      }

   }

   /*fprintf (stderr, "%f %f %f %f\n", st->adapt_rate, powE, powY, powD);*/
}

