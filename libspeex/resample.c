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
      - Variable calculation resolution depending on quality setting
         - Single vs double in float mode
         - 16-bit vs 32-bit (sinc only) in fixed-point mode
      - Make sure the filter update works even when changing params 
             after only a few samples procesed
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

#ifndef M_PI
#define M_PI 3.14159263
#endif

#ifdef FIXED_POINT
#define WORD2INT(x) ((x) < -32767 ? -32768 : ((x) > 32766 ? 32767 : (x)))  
#else
#define WORD2INT(x) ((x) < -32767.5f ? -32768 : ((x) > 32766.5f ? 32767 : floor(.5+(x))))  
#endif
               
/*#define float double*/
#define FILTER_SIZE 64
#define OVERSAMPLE 8

#define IMAX(a,b) ((a) > (b) ? (a) : (b))

struct QualityMapping {
   int base_length;
   int oversample;
   float downsample_bandwidth;
   float upsample_bandwidth;
};

/* This table maps conversion quality to internal parameters. There are two
   reasons that explain why the up-sampling bandwidth is larger than the 
   down-sampling bandwidth:
   1) When up-sampling, we can assume that the spectrum is already attenuated
      close to the Nyquist rate (from an A/D or a previous resampling filter)
   2) Any aliasing that occurs very close to the Nyquist rate will be masked
      by the sinusoids/noise just below the Nyquist rate (guaranteed only for
      up-sampling).
*/
const struct QualityMapping quality_map[11] = {
   {  8,  4, 0.70f, 0.80f}, /* 0 */
   { 16,  4, 0.74f, 0.83f}, /* 1 */
   { 32,  4, 0.77f, 0.87f}, /* 2 */
   { 48,  8, 0.84f, 0.90f}, /* 3 */
   { 64,  8, 0.88f, 0.92f}, /* 4 */
   { 80,  8, 0.90f, 0.94f}, /* 5 */
   { 96,  8, 0.91f, 0.94f}, /* 6 */
   {128, 16, 0.93f, 0.95f}, /* 7 */
   {160, 16, 0.94f, 0.96f}, /* 8 */
   {192, 16, 0.95f, 0.96f}, /* 9 */
   {256, 16, 0.96f, 0.97f}, /* 10 */
};

typedef enum {SPEEX_RESAMPLER_DIRECT_SINGLE=0, SPEEX_RESAMPLER_INTERPOLATE_SINGLE=1} SpeexSincType;

typedef int (*resampler_basic_func)(SpeexResamplerState *, int , const spx_word16_t *, int *, spx_word16_t *, int *);

struct SpeexResamplerState_ {
   int    in_rate;
   int    out_rate;
   int    num_rate;
   int    den_rate;
   
   int    quality;
   int    nb_channels;
   int    filt_len;
   int    mem_alloc_size;
   int    int_advance;
   int    frac_advance;
   float  cutoff;
   int    oversample;
   int    initialised;
   int    started;
   
   /* These are per-channel */
   int    *last_sample;
   int    *samp_frac_num;
   int    *magic_samples;
   
   spx_word16_t *mem;
   spx_word16_t *sinc_table;
   int    sinc_table_length;
   resampler_basic_func resampler_ptr;
         
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

static int resampler_basic_direct_single(SpeexResamplerState *st, int channel_index, const spx_word16_t *in, int *in_len, spx_word16_t *out, int *out_len)
{
   int N = st->filt_len;
   int out_sample = 0;
   spx_word16_t *mem;
   int last_sample = st->last_sample[channel_index];
   int samp_frac_num = st->samp_frac_num[channel_index];
   mem = st->mem + channel_index * st->mem_alloc_size;
   while (!(last_sample >= *in_len || out_sample >= *out_len))
   {
      int j;
      spx_word32_t sum=0;
      
      /* We already have all the filter coefficients pre-computed in the table */
      const spx_word16_t *ptr;
      /* Do the memory part */
      for (j=0;last_sample-N+1+j < 0;j++)
      {
         sum += MULT16_16(mem[last_sample+j],st->sinc_table[samp_frac_num*st->filt_len+j]);
      }
      
      /* Do the new part */
      ptr = in+st->in_stride*(last_sample-N+1+j);
      for (;j<N;j++)
      {
         sum += MULT16_16(*ptr,st->sinc_table[samp_frac_num*st->filt_len+j]);
         ptr += st->in_stride;
      }
   
      *out = PSHR32(sum,15);
      out += st->out_stride;
      out_sample++;
      last_sample += st->int_advance;
      samp_frac_num += st->frac_advance;
      if (samp_frac_num >= st->den_rate)
      {
         samp_frac_num -= st->den_rate;
         last_sample++;
      }
   }
   st->last_sample[channel_index] = last_sample;
   st->samp_frac_num[channel_index] = samp_frac_num;
   return out_sample;
}

static int resampler_basic_interpolate_single(SpeexResamplerState *st, int channel_index, const spx_word16_t *in, int *in_len, spx_word16_t *out, int *out_len)
{
   int N = st->filt_len;
   int out_sample = 0;
   spx_word16_t *mem;
   int last_sample = st->last_sample[channel_index];
   int samp_frac_num = st->samp_frac_num[channel_index];
   mem = st->mem + channel_index * st->mem_alloc_size;
   while (!(last_sample >= *in_len || out_sample >= *out_len))
   {
      int j;
      spx_word32_t sum=0;
      
      /* We need to interpolate the sinc filter */
      spx_word32_t accum[4] = {0.f,0.f, 0.f, 0.f};
      float interp[4];
      const spx_word16_t *ptr;
      float alpha = ((float)samp_frac_num)/st->den_rate;
      int offset = samp_frac_num*st->oversample/st->den_rate;
      float frac = alpha*st->oversample - offset;
         /* This code is written like this to make it easy to optimise with SIMD.
      For most DSPs, it would be best to split the loops in two because most DSPs 
      have only two accumulators */
      for (j=0;last_sample-N+1+j < 0;j++)
      {
         spx_word16_t curr_mem = mem[last_sample+j];
         accum[0] += MULT16_16(curr_mem,st->sinc_table[4+(j+1)*st->oversample-offset-2]);
         accum[1] += MULT16_16(curr_mem,st->sinc_table[4+(j+1)*st->oversample-offset-1]);
         accum[2] += MULT16_16(curr_mem,st->sinc_table[4+(j+1)*st->oversample-offset]);
         accum[3] += MULT16_16(curr_mem,st->sinc_table[4+(j+1)*st->oversample-offset+1]);
      }
      ptr = in+st->in_stride*(last_sample-N+1+j);
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
   
      *out = PSHR32(sum,15);
      out += st->out_stride;
      out_sample++;
      last_sample += st->int_advance;
      samp_frac_num += st->frac_advance;
      if (samp_frac_num >= st->den_rate)
      {
         samp_frac_num -= st->den_rate;
         last_sample++;
      }
   }
   st->last_sample[channel_index] = last_sample;
   st->samp_frac_num[channel_index] = samp_frac_num;
   return out_sample;
}


static void update_filter(SpeexResamplerState *st)
{
   int i;
   int old_length;
   
   old_length = st->filt_len;
   st->oversample = quality_map[st->quality].oversample;
   st->filt_len = quality_map[st->quality].base_length;
   
   if (st->num_rate > st->den_rate)
   {
      /* down-sampling */
      st->cutoff = quality_map[st->quality].downsample_bandwidth * st->den_rate / st->num_rate;
      /* FIXME: divide the numerator and denominator by a certain amount if they're too large */
      st->filt_len = st->filt_len*st->num_rate / st->den_rate;
      /* Round down to make sure we have a multiple of 4 */
      st->filt_len &= (~0x3);
   } else {
      /* up-sampling */
      st->cutoff = quality_map[st->quality].upsample_bandwidth;
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
      st->type = SPEEX_RESAMPLER_DIRECT_SINGLE;
      st->resampler_ptr = resampler_basic_direct_single;
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
      st->type = SPEEX_RESAMPLER_INTERPOLATE_SINGLE;
      st->resampler_ptr = resampler_basic_interpolate_single;
      /*fprintf (stderr, "resampler uses interpolated sinc table and normalised cutoff %f\n", cutoff);*/
   }
   st->int_advance = st->num_rate/st->den_rate;
   st->frac_advance = st->num_rate%st->den_rate;

   if (!st->mem)
   {
      st->mem = (spx_word16_t*)speex_alloc(st->nb_channels*(st->filt_len-1) * sizeof(spx_word16_t));
      for (i=0;i<st->nb_channels*(st->filt_len-1);i++)
         st->mem[i] = 0;
      st->mem_alloc_size = st->filt_len-1;
      /*speex_warning("init filter");*/
   } else if (!st->started)
   {
      st->mem = (spx_word16_t*)speex_realloc(st->mem, st->nb_channels*(st->filt_len-1) * sizeof(spx_word16_t));
      for (i=0;i<st->nb_channels*(st->filt_len-1);i++)
         st->mem[i] = 0;
      st->mem_alloc_size = st->filt_len-1;
      /*speex_warning("reinit filter");*/
   } else if (st->filt_len > old_length)
   {
      /* Increase the filter length */
      /*speex_warning("increase filter size");*/
      int old_alloc_size = st->mem_alloc_size;
      if (st->filt_len-1 > st->mem_alloc_size)
      {
         st->mem = (spx_word16_t*)speex_realloc(st->mem, st->nb_channels*(st->filt_len-1) * sizeof(spx_word16_t));
         st->mem_alloc_size = st->filt_len-1;
      }
      for (i=0;i<st->nb_channels;i++)
      {
         int j;
         /* Copy data going backward */
         for (j=0;j<old_length-1;j++)
            st->mem[i*st->mem_alloc_size+(st->filt_len-2-j)] = st->mem[i*old_alloc_size+(old_length-2-j)];
         /* Then put zeros for lack of anything better */
         for (;j<st->filt_len-1;j++)
            st->mem[i*st->mem_alloc_size+(st->filt_len-2-j)] = 0;
         /* Adjust last_sample */
         st->last_sample[i] += (st->filt_len - old_length)/2;
      }
   } else if (st->filt_len < old_length)
   {
      /* Reduce filter length, this a bit tricky */
      /*speex_warning("decrease filter size (unimplemented)");*/
      /* Adjust last_sample (which will likely end up negative) */
      /*st->last_sample += (st->filt_len - old_length)/2;*/
      for (i=0;i<st->nb_channels;i++)
      {
         int j;
         st->magic_samples[i] = (old_length - st->filt_len)/2;
         /* Copy data going backward */
         for (j=0;j<st->filt_len-1+st->magic_samples[i];j++)
            st->mem[i*st->mem_alloc_size+j] = st->mem[i*st->mem_alloc_size+j+st->magic_samples[i]];
      }
   }

}


SpeexResamplerState *speex_resampler_init(int nb_channels, int ratio_num, int ratio_den, int in_rate, int out_rate, int quality)
{
   int i;
   SpeexResamplerState *st = (SpeexResamplerState *)speex_alloc(sizeof(SpeexResamplerState));
   st->initialised = 0;
   st->started = 0;
   st->in_rate = 0;
   st->out_rate = 0;
   st->num_rate = 0;
   st->den_rate = 0;
   st->quality = -1;
   st->sinc_table_length = 0;
   st->mem_alloc_size = 0;
   st->filt_len = 0;
   st->mem = 0;
   st->resampler_ptr = 0;
         
   st->cutoff = 1.f;
   st->nb_channels = nb_channels;
   st->in_stride = 1;
   st->out_stride = 1;
   
   /* Per channel data */
   st->last_sample = (int*)speex_alloc(nb_channels*sizeof(int));
   st->magic_samples = (int*)speex_alloc(nb_channels*sizeof(int));
   st->samp_frac_num = (int*)speex_alloc(nb_channels*sizeof(int));
   for (i=0;i<nb_channels;i++)
   {
      st->last_sample[i] = 0;
      st->magic_samples[i] = 0;
      st->samp_frac_num[i] = 0;
   }

   speex_resampler_set_quality(st, quality);
   speex_resampler_set_rate(st, ratio_num, ratio_den, in_rate, out_rate);

   
   update_filter(st);
   
   st->initialised = 1;
   return st;
}

void speex_resampler_destroy(SpeexResamplerState *st)
{
   speex_free(st->mem);
   speex_free(st->sinc_table);
   speex_free(st->last_sample);
   speex_free(st->magic_samples);
   speex_free(st->samp_frac_num);
   speex_free(st);
}



static void speex_resampler_process_native(SpeexResamplerState *st, int channel_index, const spx_word16_t *in, int *in_len, spx_word16_t *out, int *out_len)
{
   int j=0;
   int N = st->filt_len;
   int out_sample = 0;
   spx_word16_t *mem;
   int tmp_out_len = 0;
   mem = st->mem + channel_index * st->mem_alloc_size;
   st->started = 1;
   
   /* Handle the case where we have samples left from a reduction in filter length */
   if (st->magic_samples)
   {
      int tmp_in_len;
      tmp_in_len = st->magic_samples[channel_index];
      tmp_out_len = *out_len;
      /* FIXME: Need to handle the case where the out array is too small */
      /* magic_samples needs to be set to zero to avoid infinite recursion */
      st->magic_samples = 0;
      speex_resampler_process_native(st, channel_index, mem+N-1, &tmp_in_len, out, &tmp_out_len);
      /*speex_warning_int("extra samples:", tmp_out_len);*/
      out += tmp_out_len;
   }
   
   /* Call the right resampler through the function ptr */
   out_sample = st->resampler_ptr(st, channel_index, in, in_len, out, out_len);
   
   if (st->last_sample[channel_index] < *in_len)
      *in_len = st->last_sample[channel_index];
   *out_len = out_sample+tmp_out_len;
   st->last_sample[channel_index] -= *in_len;
   
   for (j=0;j<N-1-*in_len;j++)
      mem[j] = mem[j+*in_len];
   for (;j<N-1;j++)
      mem[j] = in[st->in_stride*(j+*in_len-N+1)];
   
}

#ifdef FIXED_POINT
void speex_resampler_process_float(SpeexResamplerState *st, int channel_index, const float *in, int *in_len, float *out, int *out_len)
{
   int i;
   int istride_save, ostride_save;
   spx_word16_t x[*in_len];
   spx_word16_t y[*out_len];
   istride_save = st->in_stride;
   ostride_save = st->out_stride;
   for (i=0;i<*in_len;i++)
      x[i] = WORD2INT(in[i*st->in_stride]);
   st->in_stride = st->out_stride = 1;
   speex_resampler_process_native(st, channel_index, x, in_len, y, out_len);
   st->in_stride = istride_save;
   st->out_stride = ostride_save;
   for (i=0;i<*out_len;i++)
      out[i*st->out_stride] = y[i];
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
   int istride_save, ostride_save;
   spx_word16_t x[*in_len];
   spx_word16_t y[*out_len];
   istride_save = st->in_stride;
   ostride_save = st->out_stride;
   for (i=0;i<*in_len;i++)
      x[i] = in[i+st->in_stride];
   st->in_stride = st->out_stride = 1;
   speex_resampler_process_native(st, channel_index, x, in_len, y, out_len);
   st->in_stride = istride_save;
   st->out_stride = ostride_save;
   for (i=0;i<*out_len;i++)
      out[i+st->out_stride] = WORD2INT(y[i]);
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


void speex_resampler_set_rate(SpeexResamplerState *st, int ratio_num, int ratio_den, int in_rate, int out_rate)
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

void speex_resampler_set_quality(SpeexResamplerState *st, int quality)
{
   if (quality < 0)
      quality = 0;
   if (quality > 10)
      quality = 10;
   if (st->quality == quality)
      return;
   st->quality = quality;
   if (st->initialised)
      update_filter(st);
}

void speex_resampler_set_input_stride(SpeexResamplerState *st, int stride)
{
   st->in_stride = stride;
}

void speex_resampler_set_output_stride(SpeexResamplerState *st, int stride)
{
   st->out_stride = stride;
}

void speex_resampler_skip_zeros(SpeexResamplerState *st)
{
   int i;
   for (i=0;i<st->nb_channels;i++)
      st->last_sample[i] = st->filt_len/2;
}

void speex_resampler_reset_mem(SpeexResamplerState *st)
{
   int i;
   for (i=0;i<st->nb_channels*(st->filt_len-1);i++)
      st->mem[i] = 0;
}

