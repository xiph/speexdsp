/* Copyright (C) 2007 Jean-Marc Valin
      
   File: resample.c
   Arbitrary resampling code

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
   The design goals of this code are:
      - Very fast algorithm
      - SIMD-friendly algorithm
      - Low memory requirement
      - Good *perceptual* quality (and not best SNR)

   The code is working, but it's in a very early stage, so it may have
   artifacts, noise or subliminal messages from satan. Also, the API 
   isn't stable and I can actually promise that I *will* change the API
   some time in the future.

TODO list:
      - Variable length filter (depending on frequency/ratio and quality)
      - Quality setting to control filter length (and sinc window?)
      - Variable calculation resolution depending on quality setting
         - Single vs double in float mode
         - 16-bit vs 32-bit (sinc only) in fixed-point mode
      - Make it possible to change the filter length without major artifacts

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef OUTSIDE_SPEEX
#include <stdlib.h>
void *speex_alloc (int size) {return calloc(size,1);}
void *speex_realloc (void *ptr, int size) {return realloc(ptr, size);}
void speex_free (void *ptr) {free(ptr);}
#else
#include "misc.h"
#endif


#include <math.h>
#include "speex/speex_resampler.h"

#ifdef FIXED_POINT
#define WORD2INT(x) ((x) < -32767 ? -32768 : ((x) > 32766 ? 32767 : (x)))  
#else
#define WORD2INT(x) ((x) < -32767.5f ? -32768 : ((x) > 32766.5f ? 32767 : floor(.5+(x))))  
#endif
               
/*#define float double*/
#define FILTER_SIZE 64
#define OVERSAMPLE 8

#define IMAX(a,b) ((a) > (b) ? (a) : (b))

typedef enum {SPEEX_RESAMPLER_DIRECT=0, SPEEX_RESAMPLER_INTERPOLATE=1} SpeexSincType;

struct SpeexResamplerState_ {
   int    in_rate;
   int    out_rate;
   int    num_rate;
   int    den_rate;
   
   int    quality;
   int    nb_channels;
   int    last_sample;
   int    samp_frac_num;
   int    filt_len;
   int    int_advance;
   int    frac_advance;
   float  cutoff;
   int    oversample;
   int    initialised;
   
   
   spx_word16_t *mem;
   spx_word16_t *sinc_table;
   int    sinc_table_length;
   int    in_stride;
   int    out_stride;
   SpeexSincType type;
} ;

#ifdef FIXED_POINT
/* The slow way of computing a sinc for the table. Should improve that some day */
static spx_word16_t sinc(float cutoff, float x, int N)
{
   /*fprintf (stderr, "%f ", x);*/
   x *= cutoff;
   if (fabs(x)<1e-6f)
      return WORD2INT(32768.*cutoff);
   else if (fabs(x) > .5f*N)
      return 0;
   /*FIXME: Can it really be any slower than this? */
   return WORD2INT(32768.*cutoff*sin(M_PI*x)/(M_PI*x) * (.42+.5*cos(2*x*M_PI/N)+.08*cos(4*x*M_PI/N)));
}
#else
/* The slow way of computing a sinc for the table. Should improve that some day */
static spx_word16_t sinc(float cutoff, float x, int N)
{
   /*fprintf (stderr, "%f ", x);*/
   x *= cutoff;
   if (fabs(x)<1e-6)
      return cutoff;
   else if (fabs(x) > .5*N)
      return 0;
   /*FIXME: Can it really be any slower than this? */
   return cutoff*sin(M_PI*x)/(M_PI*x) * (.42+.5*cos(2*x*M_PI/N)+.08*cos(4*x*M_PI/N));
}
#endif

static void update_filter(SpeexResamplerState *st)
{
   int i;
   
   st->oversample = OVERSAMPLE;
   if (st->quality > 7)
      st->oversample *= 2;
   st->filt_len = 8 + 12*st->quality;
   
   if (st->num_rate > st->den_rate)
   {
      /* down-sampling */
      st->cutoff = .92f * st->den_rate / st->num_rate;
      st->filt_len *= st->num_rate / st->den_rate;
   } else {
      /* up-sampling */
      st->cutoff = .97f;
   }

   /* Choose the resampling type that requires the least amount of memory */
   if (st->den_rate <= st->oversample)
   {
      if (!st->sinc_table)
         st->sinc_table = (spx_word16_t *)speex_alloc(st->filt_len*st->den_rate*sizeof(spx_word16_t));
      else if (st->sinc_table_length < st->filt_len*st->den_rate)
      {
         st->sinc_table = (spx_word16_t *)speex_realloc(st->sinc_table,st->filt_len*st->den_rate*sizeof(spx_word16_t));
         st->sinc_table_length = st->filt_len*st->den_rate;
      }
      for (i=0;i<st->den_rate;i++)
      {
         int j;
         for (j=0;j<st->filt_len;j++)
         {
            st->sinc_table[i*st->filt_len+j] = sinc(st->cutoff,((j-st->filt_len/2+1)-((float)i)/st->den_rate), st->filt_len);
         }
      }
      st->type = SPEEX_RESAMPLER_DIRECT;
      /*fprintf (stderr, "resampler uses direct sinc table and normalised cutoff %f\n", cutoff);*/
   } else {
      if (!st->sinc_table)
         st->sinc_table = (spx_word16_t *)speex_alloc((st->filt_len*st->oversample+8)*sizeof(spx_word16_t));
      else if (st->sinc_table_length < st->filt_len*st->oversample+8)
      {
         st->sinc_table = (spx_word16_t *)speex_realloc(st->sinc_table,(st->filt_len*st->oversample+8)*sizeof(spx_word16_t));
         st->sinc_table_length = st->filt_len*st->oversample+8;
      }
      for (i=-4;i<st->oversample*st->filt_len+4;i++)
         st->sinc_table[i+4] = sinc(st->cutoff,(i/(float)st->oversample - st->filt_len/2), st->filt_len);
      st->type = SPEEX_RESAMPLER_INTERPOLATE;
      /*fprintf (stderr, "resampler uses interpolated sinc table and normalised cutoff %f\n", cutoff);*/
   }
   st->int_advance = st->num_rate/st->den_rate;
   st->frac_advance = st->num_rate%st->den_rate;

   if (!st->mem)
      st->mem = (spx_word16_t*)speex_alloc(st->nb_channels*(st->filt_len-1) * sizeof(spx_word16_t));
   else
      st->mem = (spx_word16_t*)speex_realloc(st->mem, st->nb_channels*(st->filt_len-1) * sizeof(spx_word16_t));
   for (i=0;i<st->nb_channels*(st->filt_len-1);i++)
      st->mem[i] = 0;

}


SpeexResamplerState *speex_resampler_init(int nb_channels, int ratio_num, int ratio_den, int in_rate, int out_rate, int quality)
{
   SpeexResamplerState *st = (SpeexResamplerState *)speex_alloc(sizeof(SpeexResamplerState));
   st->initialised = 0;
   st->in_rate = 0;
   st->out_rate = 0;
   st->num_rate = 0;
   st->den_rate = 0;
   st->quality = -1;
   st->sinc_table_length = 0;
   st->mem = 0;
   
   st->cutoff = 1.f;
   st->nb_channels = nb_channels;
   st->last_sample = 0;
   
   st->in_stride = 1;
   st->out_stride = 1;

   speex_resample_set_quality(st, quality);
   speex_resample_set_rate(st, ratio_num, ratio_den, in_rate, out_rate);

   
   update_filter(st);
   
   st->initialised = 1;
   return st;
}

void speex_resampler_destroy(SpeexResamplerState *st)
{
   speex_free(st->mem);
   speex_free(st->sinc_table);
   speex_free(st);
}

static void speex_resampler_process_native(SpeexResamplerState *st, int channel_index, const spx_word16_t *in, int *in_len, spx_word16_t *out, int *out_len)
{
   int j=0;
   int N = st->filt_len;
   int out_sample = 0;
   spx_word16_t *mem;
   mem = st->mem + channel_index * (N-1);
   while (!(st->last_sample >= *in_len || out_sample >= *out_len))
   {
      int j;
      spx_word32_t sum=0;
      
      if (st->type == SPEEX_RESAMPLER_DIRECT)
      {
         /* We already have all the filter coefficients pre-computed in the table */
         const spx_word16_t *ptr;
         /* Do the memory part */
         for (j=0;st->last_sample-N+1+j < 0;j++)
         {
            sum += MULT16_16(mem[st->last_sample+j],st->sinc_table[st->samp_frac_num*st->filt_len+j]);
         }
         
         /* Do the new part */
         ptr = in+st->last_sample-N+1+j;
         for (;j<N;j++)
         {
            sum += MULT16_16(*ptr,st->sinc_table[st->samp_frac_num*st->filt_len+j]);
            ptr += st->in_stride;
         }
      } else {
         /* We need to interpolate the sinc filter */
         spx_word32_t accum[4] = {0.f,0.f, 0.f, 0.f};
         float interp[4];
         const spx_word16_t *ptr;
         float alpha = ((float)st->samp_frac_num)/st->den_rate;
         int offset = st->samp_frac_num*st->oversample/st->den_rate;
         float frac = alpha*st->oversample - offset;
         /* This code is written like this to make it easy to optimise with SIMD.
            For most DSPs, it would be best to split the loops in two because most DSPs 
            have only two accumulators */
         for (j=0;st->last_sample-N+1+j < 0;j++)
         {
            spx_word16_t curr_mem = mem[st->last_sample+j];
            accum[0] += MULT16_16(curr_mem,st->sinc_table[4+(j+1)*st->oversample-offset-2]);
            accum[1] += MULT16_16(curr_mem,st->sinc_table[4+(j+1)*st->oversample-offset-1]);
            accum[2] += MULT16_16(curr_mem,st->sinc_table[4+(j+1)*st->oversample-offset]);
            accum[3] += MULT16_16(curr_mem,st->sinc_table[4+(j+1)*st->oversample-offset+1]);
         }
         ptr = in+st->last_sample-N+1+j;
         /* Do the new part */
         for (;j<N;j++)
         {
            spx_word16_t curr_in = *ptr;
            ptr += st->in_stride;
            accum[0] += MULT16_16(curr_in,st->sinc_table[4+(j+1)*st->oversample-offset-2]);
            accum[1] += MULT16_16(curr_in,st->sinc_table[4+(j+1)*st->oversample-offset-1]);
            accum[2] += MULT16_16(curr_in,st->sinc_table[4+(j+1)*st->oversample-offset]);
            accum[3] += MULT16_16(curr_in,st->sinc_table[4+(j+1)*st->oversample-offset+1]);
         }
         /* Compute interpolation coefficients. I'm not sure whether this corresponds to cubic interpolation
            but I know it's MMSE-optimal on a sinc */
         interp[0] =  -0.16667f*frac + 0.16667f*frac*frac*frac;
         interp[1] = frac + 0.5f*frac*frac - 0.5f*frac*frac*frac;
         /*interp[2] = 1.f - 0.5f*frac - frac*frac + 0.5f*frac*frac*frac;*/
         interp[3] = -0.33333f*frac + 0.5f*frac*frac - 0.16667f*frac*frac*frac;
         /* Just to make sure we don't have rounding problems */
         interp[2] = 1.f-interp[0]-interp[1]-interp[3];
         /*sum = frac*accum[1] + (1-frac)*accum[2];*/
         sum = interp[0]*accum[0] + interp[1]*accum[1] + interp[2]*accum[2] + interp[3]*accum[3];
      }
      *out = PSHR32(sum,15);
      out += st->out_stride;
      out_sample++;
      st->last_sample += st->int_advance;
      st->samp_frac_num += st->frac_advance;
      if (st->samp_frac_num >= st->den_rate)
      {
         st->samp_frac_num -= st->den_rate;
         st->last_sample++;
      }
   }
   if (st->last_sample < *in_len)
      *in_len = st->last_sample;
   *out_len = out_sample;
   st->last_sample -= *in_len;
   
   for (j=0;j<N-1-*in_len;j++)
      mem[j] = mem[j+*in_len];
   for (;j<N-1;j++)
      mem[j] = in[st->in_stride*(j+*in_len-N+1)];
   
}

#ifdef FIXED_POINT
void speex_resampler_process_float(SpeexResamplerState *st, int channel_index, const float *in, int *in_len, float *out, int *out_len)
{
   int i;
   spx_word16_t x[*in_len];
   spx_word16_t y[*out_len];
   for (i=0;i<*in_len;i++)
      x[i] = WORD2INT(in[i]);
   speex_resampler_process_native(st, channel_index, x, in_len, y, out_len);
   for (i=0;i<*out_len;i++)
      out[i] = y[i];

}
void speex_resampler_process_int(SpeexResamplerState *st, int channel_index, const spx_int16_t *in, int *in_len, spx_int16_t *out, int *out_len)
{
   speex_resampler_process_native(st, channel_index, in, in_len, out, out_len);
}
#else
void speex_resampler_process_float(SpeexResamplerState *st, int channel_index, const float *in, int *in_len, float *out, int *out_len)
{
   speex_resampler_process_native(st, channel_index, in, in_len, out, out_len);
}
void speex_resampler_process_int(SpeexResamplerState *st, int channel_index, const spx_int16_t *in, int *in_len, spx_int16_t *out, int *out_len)
{
   int i;
   spx_word16_t x[*in_len];
   spx_word16_t y[*out_len];
   for (i=0;i<*out_len;i++)
      x[i] = in[i];
   speex_resampler_process_native(st, channel_index, x, in_len, y, out_len);
   for (i=0;i<*in_len;i++)
      out[i] = WORD2INT(y[i]);
}
#endif

void speex_resampler_process_interleaved_float(SpeexResamplerState *st, const float *in, int *in_len, float *out, int *out_len)
{
   int i;
   int istride_save, ostride_save;
   istride_save = st->in_stride;
   ostride_save = st->out_stride;
   st->in_stride = st->out_stride = st->nb_channels;
   for (i=0;i<st->nb_channels;i++)
   {
      speex_resampler_process_float(st, i, in+i, in_len, out+i, out_len);
   }
   st->in_stride = istride_save;
   st->out_stride = ostride_save;
}


void speex_resample_set_rate(SpeexResamplerState *st, int ratio_num, int ratio_den, int in_rate, int out_rate)
{
   int fact;
   if (st->in_rate == in_rate && st->out_rate == out_rate && st->num_rate == ratio_num && st->den_rate == ratio_den)
      return;
   
   st->in_rate = in_rate;
   st->out_rate = out_rate;
   st->num_rate = ratio_num;
   st->den_rate = ratio_den;
   /* FIXME: This is terribly inefficient, but who cares (at least for now)? */
   for (fact=2;fact<=sqrt(IMAX(in_rate, out_rate));fact++)
   {
      while ((st->num_rate % fact == 0) && (st->den_rate % fact == 0))
      {
         st->num_rate /= fact;
         st->den_rate /= fact;
      }
   }
      
   if (st->initialised)
      update_filter(st);
}

void speex_resample_set_quality(SpeexResamplerState *st, int quality)
{
   if (st->quality == quality)
      return;
   st->quality = quality;
   if (st->initialised)
      update_filter(st);
}

void speex_resample_set_input_stride(SpeexResamplerState *st, int stride)
{
   st->in_stride = stride;
}

void speex_resample_set_output_stride(SpeexResamplerState *st, int stride)
{
   st->out_stride = stride;
}

void speex_resample_skip_zeros(SpeexResamplerState *st)
{
   st->last_sample = st->filt_len/2;
}

void speex_resample_reset_mem(SpeexResamplerState *st)
{
   int i;
   for (i=0;i<st->nb_channels*(st->filt_len-1);i++)
      st->mem[i] = 0;
}

