/* Copyright (C) 2007 Jean-Marc Valin
      
   File: resample.c
   Resampling code
      
   The design goals of this code are:
      - Very fast algorithm
      - Low memory requirement
      - Good *perceptual* quality (and not best SNR)

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
#include <math.h>
#include <stdio.h>
            
/*#define float double*/
#define FILTER_SIZE 64
#define OVERSAMPLE 8

typedef enum {SPEEX_RESAMPLER_DIRECT=0, SPEEX_RESAMPLER_INTERPOLATE=1} SpeexSincType;

typedef struct {
   int in_rate;
   int out_rate;
   int num_rate;
   int den_rate;
   int last_sample;
   int samp_frac_num;
   int filt_len;
   float *mem;
   float *sinc_table;
   SpeexSincType type;
} SpeexResamplerState;

static float sinc(float x, int N)
{
   /*fprintf (stderr, "%f ", x);*/
   if (fabs(x)<1e-6)
      return 1;
   else if (fabs(x) > .5f*N)
      return 0;
   /*FIXME: Can it really be any slower than this? */
   return sin(M_PI*x)/(M_PI*x) * (.5+.5*cos(2*x*M_PI/N));
}

SpeexResamplerState *speex_resampler_init(int nb_channels, int in_rate, int out_rate, int in_rate_den, int out_rate_den)
{
   int fact, i;
   float cutoff;
   SpeexResamplerState *st = (SpeexResamplerState *)speex_alloc(sizeof(SpeexResamplerState));
   st->in_rate = in_rate;
   st->out_rate = out_rate;
   st->num_rate = in_rate;
   st->den_rate = out_rate;
   /* FIXME: This is terribly inefficient, but who cares (at least for now)? */
   for (fact=2;fact<=sqrt(MAX32(in_rate, out_rate));fact++)
   {
      while ((st->num_rate % fact == 0) && (st->den_rate % fact == 0))
      {
         st->num_rate /= fact;
         st->den_rate /= fact;
      }
   }
   st->last_sample = 0;
   st->filt_len = FILTER_SIZE;
   st->mem = (float*)speex_alloc((st->filt_len-1) * sizeof(float));
   for (i=0;i<st->filt_len-1;i++)
      st->mem[i] = 0;
   /* FIXME: Is there a danger of overflow? */
   if (in_rate*out_rate_den > out_rate*in_rate_den)
   {
      /* down-sampling */
      cutoff = .92f * out_rate*in_rate_den / (in_rate*out_rate_den);
   } else {
      /* up-sampling */
      cutoff = .97;
   }
   if (st->den_rate <= OVERSAMPLE)
   {
      st->sinc_table = (float *)speex_alloc(st->filt_len*st->den_rate*sizeof(float));
      for (i=0;i<st->den_rate;i++)
      {
         int j;
         for (j=0;j<st->filt_len;j++)
         {
            st->sinc_table[i*st->filt_len+j] = sinc(cutoff*((j-st->filt_len/2+1)-((float)i)/st->den_rate), st->filt_len);
         }
      }
      st->type = SPEEX_RESAMPLER_DIRECT;
      fprintf (stderr, "resampler uses direct sinc table and normalised cutoff %f\n", cutoff);
   } else {
      st->sinc_table = (float *)speex_alloc(st->filt_len*st->den_rate*sizeof(float));
      for (i=-4;i<OVERSAMPLE*st->filt_len+4;i++)
         st->sinc_table[i+4] = sinc(cutoff*(i/(float)OVERSAMPLE - st->filt_len/2), st->filt_len);
      st->type = SPEEX_RESAMPLER_INTERPOLATE;
      fprintf (stderr, "resampler uses interpolated sinc table and normalised cutoff %f\n", cutoff);
   }
   return st;
}

void speex_resampler_destroy(SpeexResamplerState *st)
{
   speex_free(st->mem);
   if (st->sinc_table)
      speex_free(st->sinc_table);
   speex_free(st);
}

void speex_resample_float(SpeexResamplerState *st, int index, const float *in, int *in_len, float *out, int *out_len)
{
   int j=0;
   int N = st->filt_len;
   int out_sample = 0;
   while (1)
   {
      int j;
      float sum=0;
      
      if (st->type == SPEEX_RESAMPLER_DIRECT)
      {
         /* We already have all the filter coefficients pre-computed in the table */
         /* Do the memory part */
         for (j=0;st->last_sample-N+1+j < 0;j++)
         {
            sum += st->mem[st->last_sample+j]*st->sinc_table[st->samp_frac_num*st->filt_len+j];
         }
         /* Do the new part */
         for (;j<N;j++)
         {
            sum += in[st->last_sample-N+1+j]*st->sinc_table[st->samp_frac_num*st->filt_len+j];
         }
      } else {
         /* We need to interpolate the sinc filter */
         float accum[4] = {0.f,0.f, 0.f, 0.f};
         float interp[4];
         float alpha = ((float)st->samp_frac_num)/st->den_rate;
         int offset = st->samp_frac_num*OVERSAMPLE/st->den_rate;
         float frac = alpha*OVERSAMPLE - offset;
         /* This code is written like this to make it easy to optimise with SIMD.
            For most DSPs, it would be best to split the loops in two because most DSPs 
            have only two accumulators */
         for (j=0;st->last_sample-N+1+j < 0;j++)
         {
            accum[0] += st->mem[st->last_sample+j]*st->sinc_table[4+(j+1)*OVERSAMPLE-offset-2];
            accum[1] += st->mem[st->last_sample+j]*st->sinc_table[4+(j+1)*OVERSAMPLE-offset-1];
            accum[2] += st->mem[st->last_sample+j]*st->sinc_table[4+(j+1)*OVERSAMPLE-offset];
            accum[3] += st->mem[st->last_sample+j]*st->sinc_table[4+(j+1)*OVERSAMPLE-offset+1];
         }
         /* Do the new part */
         for (;j<N;j++)
         {
            accum[0] += in[st->last_sample-N+1+j]*st->sinc_table[4+(j+1)*OVERSAMPLE-offset-2];
            accum[1] += in[st->last_sample-N+1+j]*st->sinc_table[4+(j+1)*OVERSAMPLE-offset-1];
            accum[2] += in[st->last_sample-N+1+j]*st->sinc_table[4+(j+1)*OVERSAMPLE-offset];
            accum[3] += in[st->last_sample-N+1+j]*st->sinc_table[4+(j+1)*OVERSAMPLE-offset+1];
         }
         /* Compute interpolation coefficients. I'm not sure whether this corresponds to cubic interpolation
            but I know it's MMSE-optimal on a sinc */
         interp[0] =  -0.16667f*frac + 0.16667f*frac*frac*frac;
         interp[1] = frac + 0.5f*frac*frac - 0.5f*frac*frac*frac;
         interp[2] = 1.f - 0.5f*frac - frac*frac + 0.5f*frac*frac*frac;
         interp[3] = -0.33333f*frac + 0.5f*frac*frac - 0.16667f*frac*frac*frac;
         /*sum = frac*accum[1] + (1-frac)*accum[2];*/
         sum = interp[0]*accum[0] + interp[1]*accum[1] + interp[2]*accum[2] + interp[3]*accum[3];
      }
      out[out_sample++] = sum;
      
      st->last_sample += st->num_rate/st->den_rate;
      st->samp_frac_num += st->num_rate%st->den_rate;
      if (st->samp_frac_num >= st->den_rate)
      {
         st->samp_frac_num -= st->den_rate;
         st->last_sample++;
      }
      if (st->last_sample >= *in_len || out_sample >= *out_len)
         break;
   }
   if (st->last_sample < *in_len)
      *in_len = st->last_sample;
   *out_len = out_sample;
   st->last_sample -= *in_len;
   
   /* FIXME: The details of this are untested */
   for (j=0;j<N-1-*in_len;j++)
      st->mem[j] = st->mem[j-*in_len];
   for (;j<N-1;j++)
      st->mem[j] = in[j+*in_len-N+1];
   
}

void speex_resample_set_rate(SpeexResamplerState *st, int in_rate, int out_rate, int in_rate_den, int out_rate_den);

void speex_resample_set_input_stride(SpeexResamplerState *st, int stride);

void speex_resample_set_output_stride(SpeexResamplerState *st, int stride);

void speex_resample_skip_zeros(SpeexResamplerState *st);


#define NN 256

int main(int argc, char **argv)
{
   int i;
   short *in;
   short *out;
   float *fin, *fout;
   SpeexResamplerState *st = speex_resampler_init(1, 8000, 13501, 1, 1);
   in = speex_alloc(NN*sizeof(short));
   out = speex_alloc(2*NN*sizeof(short));
   fin = speex_alloc(NN*sizeof(float));
   fout = speex_alloc(2*NN*sizeof(float));
   while (1)
   {
      int in_len;
      int out_len;
      fread(in, sizeof(short), NN, stdin);
      if (feof(stdin))
         break;
      for (i=0;i<NN;i++)
         fin[i]=in[i];
      in_len = NN;
      out_len = 2*NN;
      speex_resample_float(st, 0, fin, &in_len, fout, &out_len);
      for (i=0;i<out_len;i++)
         out[i]=floor(.5+fout[i]);
      fwrite(out, sizeof(short), out_len, stdout);
   }
   speex_resampler_destroy(st);
   speex_free(in);
   speex_free(out);
   speex_free(fin);
   speex_free(fout);
   return 0;
}

