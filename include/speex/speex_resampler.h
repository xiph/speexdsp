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

#define SPEEX_RESAMPLER_QUALITY_MAX 10
#define SPEEX_RESAMPLER_QUALITY_MIN 0
#define SPEEX_RESAMPLER_QUALITY_DEFAULT 4
#define SPEEX_RESAMPLER_QUALITY_VOIP 3
#define SPEEX_RESAMPLER_QUALITY_DESKTOP 5
   
struct SpeexResamplerState_;
typedef struct SpeexResamplerState_ SpeexResamplerState;

/** Create a new resampler. The sampling rate ratio is an arbitrary rational number 
 * with both the numerator and denominator being 32-bit integers.
 * @param nb_channels Number of channels to be processed
 * @param ratio_num Numerator of the sampling rate ratio
 * @param ratio_den Denominator of the sampling rate ratio
 * @param in_rate Nominal input sampling rate rounded to the nearest integer (in Hz)
 * @param out_rate Nominal output sampling rate rounded to the nearest integer (in Hz)
 * @param quality Resampling quality (0-10)
 * @return Newly created resampler state
 */
SpeexResamplerState *speex_resampler_init(int nb_channels, int ratio_num, int ratio_den, int in_rate, int out_rate, int quality);

/** Destroy a resampler state.
 * @param st Resampler state
 */
void speex_resampler_destroy(SpeexResamplerState *st);

/** Resample a float array.
 * @param st Resampler state
 * @param channel_index Index of the channel to process for the multi-channel base (0 otherwise)
 * @param in Input buffer
 * @param in_len Number of input samples in the input buffer. Returns the number of samples processed
 * @param out Output buffer
 * @param out_len Size of the output buffer. Returns the number of samples written
 */
void speex_resampler_process_float(SpeexResamplerState *st, int channel_index, const float *in, int *in_len, float *out, int *out_len);

/** Resample an int array.
 * @param st Resampler state
 * @param channel_index Index of the channel to process for the multi-channel base (0 otherwise)
 * @param in Input buffer
 * @param in_len Number of input samples in the input buffer. Returns the number of samples processed
 * @param out Output buffer
 * @param out_len Size of the output buffer. Returns the number of samples written
 */
void speex_resampler_process_int(SpeexResamplerState *st, int channel_index, const spx_int16_t *in, int *in_len, spx_int16_t *out, int *out_len);

/** Resample an interleaved float array.
 * @param st Resampler state
 * @param in Input buffer
 * @param in_len Number of input samples in the input buffer. Returns the number of samples processed. This is all per-channel.
 * @param out Output buffer
 * @param out_len Size of the output buffer. Returns the number of samples written. This is all per-channel.
 */
void speex_resampler_process_interleaved_float(SpeexResamplerState *st, const float *in, int *in_len, float *out, int *out_len);

/** Set (change) the input/output sampling rates and resampling ratio.
 * @param st Resampler state
 * @param ratio_num Numerator of the sampling rate ratio
 * @param ratio_den Denominator of the sampling rate ratio
 * @param in_rate Nominal input sampling rate rounded to the nearest integer (in Hz)
 * @param out_rate Nominal output sampling rate rounded to the nearest integer (in Hz)
 */
void speex_resampler_set_rate(SpeexResamplerState *st, int ratio_num, int ratio_den, int in_rate, int out_rate);

/** Set (change) the conversion quality
 * @param st Resampler state
 * @param quality Resampling quality (0-10)
 */
void speex_resampler_set_quality(SpeexResamplerState *st, int quality);

/** Set (change) the input stride
 * @param st Resampler state
 * @param stride Input stride
 */
void speex_resampler_set_input_stride(SpeexResamplerState *st, int stride);

/** Set (change) the output stride
 * @param st Resampler state
 * @param stride Output stride
 */
void speex_resample_set_output_stride(SpeexResamplerState *st, int stride);

/** Make sure that the first samples to go out of the resamplers don't have leading zeros.
 * This is only useful before starting to use a newly created resampler.
 * @param st Resampler state
 */
void speex_resampler_skip_zeros(SpeexResamplerState *st);

/** Reset a resampler so a new (unrelated) stream can be processed
 * @param st Resampler state
 */
void speex_resampler_reset_mem(SpeexResamplerState *st);

#ifdef __cplusplus
}
#endif

#endif
