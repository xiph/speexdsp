/* Copyright (C) 2007 Jean-Marc Valin
      
   File: resample.c
   Resample code

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
      
typedef struct {
   int in_rate;
   int out_rate;
   int num_rate;
   int den_rate;
   int last_sample;
   int samp_frac_num;
   int filt_len;
   float *mem;
} SpeexResamplerState;


SpeexResamplerState *speex_resampler_init(int in_rate, int out_rate)
{
   SpeexResamplerState *st = (SpeexResamplerState *)speex_alloc(sizeof(SpeexResamplerState));
   int fact, i;
   st->in_rate = in_rate;
   st->out_rate = out_rate;
   st->num_rate = in_rate;
   st->den_rate = out_rate;
   /* FIXME: This is terribly inefficient, but who cares (at least for now)? */
   for (fact=2;fact<=spx_sqrt(MAX32(in_rate, out_rate));fact++)
   {
      while ((st->num_rate % fact == 0) && (st->den_rate % fact == 0))
      {
         st->num_rate /= fact;
         st->den_rate /= fact;
      }
   }
   st->last_sample = 0;
   st->filt_len = 64;
   st->mem = speex_alloc((st->filt_len-1) * sizeof(float));
   for (i=0;i<st->filt_len-1;i++)
      st->mem[i] = 0;
   return st;
}

void speex_resampler_init(SpeexResamplerState *st)
{
   speex_free(st->mem);
   speex_free(st);
}

int speex_resample_float(SpeexResamplerState *st, const float *in, int len, float *out)
{
   int i=0;
   int N = st->filt_len;
   while (1)
   {
      int j;
      float sum=0;
      for (j=0;j<N;j++)
      {
         sum += in[st->last_sample-N+1+j]*sinc((j-N/2)-((float)st->samp_frac_num)/st->den_rate);
      }
      out[i++] = sum;
      
      st->last_sample += st->num_rate/st->den_rate;
      st->samp_frac_num += st->num_rate%st->den_rate;
      if (st->samp_frac_num >= st->den_rate)
      {
         st->samp_frac_num >= st->den_rate;
         st->last_sample++;
      }
      if (st->last_sample >= len)
      {
         break;
      }      
   }
   return i;
}

