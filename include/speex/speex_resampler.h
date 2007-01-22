/* Copyright (C) 2007 Jean-Marc Valin
      
   File: speex_resampler.h
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


#ifndef SPEEX_RESAMPLER_H
#define SPEEX_RESAMPLER_H

#ifdef OUTSIDE_SPEEX

#define spx_int16_t short
#ifdef FIXED_POINT
#define spx_word16_t short
#define spx_word32_t int
#else
#define spx_word16_t float
#define spx_word32_t float
#define MULT16_16(a,b) ((a)*(b))
#define PSHR32(a,b) (a)
#endif

#else

#include "speex/speex_types.h"

#endif

#ifdef __cplusplus
extern "C" {
#endif

struct SpeexResamplerState_;
typedef struct SpeexResamplerState_ SpeexResamplerState;
//typedef SpeexResamplerState;

SpeexResamplerState *speex_resampler_init(int nb_channels, int in_rate, int out_rate, int in_rate_den, int out_rate_den);

void speex_resampler_destroy(SpeexResamplerState *st);

void speex_resampler_process_float(SpeexResamplerState *st, int channel_index, const float *in, int *in_len, float *out, int *out_len);

void speex_resampler_process_int(SpeexResamplerState *st, int channel_index, const spx_int16_t *in, int *in_len, spx_int16_t *out, int *out_len);

void speex_resampler_process_interleaved_float(SpeexResamplerState *st, const float *in, int *in_len, float *out, int *out_len);

void speex_resample_set_rate(SpeexResamplerState *st, int in_rate, int out_rate, int in_rate_den, int out_rate_den);

void speex_resample_set_input_stride(SpeexResamplerState *st, int stride);

void speex_resample_set_output_stride(SpeexResamplerState *st, int stride);

void speex_resample_skip_zeros(SpeexResamplerState *st);

void speex_resample_reset_mem(SpeexResamplerState *st);

#ifdef __cplusplus
}
#endif

#endif
