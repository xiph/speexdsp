/* Copyright (C) 2003 Epic Games 
   Written by Jean-Marc Valin

   File: speex_denoise.h


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


#include "smallft.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SpeexDenoiseState {
   int frame_size;           /**< Number of samples processed each time */
   int ps_size;              /**< Number of points in the power spectrum */

   float *frame;             /**< Processing frame (2*ps_size) */
   float *ps;                /**< Current power spectrum */
   float *gain2;             /**< Adjusted gains */
   float *window;            /**< Analysis/Synthesis window */
   float *noise;             /**< Noise estimate */
   float *old_ps;            /**< Power spectrum for last frame */
   float *gain;              /**< Ephraim Malah gain */
   float *prior;             /**< A-priori SNR */
   float *post;              /**< A-posteriori SNR */
   float *min_ps;            /**< */
   float *last_energy;       /**< Energy of the previous frames */
   float *last_ps;           /**< Power spectrum of the past frames */
   float *loudness_weight;   /**< */
   int    last_id;           /**< */

   float *inbuf;             /**< Input buffer (overlapped analysis) */
   float *outbuf;            /**< Output buffer (for overlap and add) */

   float  speech_prob;
   int    last_speech;
   float  loudness;          /**< loudness estimate */
   float  loudness2;          /**< loudness estimate */
   int    nb_adapt;          /**< Number of frames used for adaptation so far */
   int    nb_loudness_adapt; /**< Number of frames used for loudness adaptation so far */
   int    consec_noise;      /**< Number of consecutive noise frames */
   int    nb_denoise;        /**< Number of frames processed so far */
   int    nb_min_estimate;   /**< */
   int    last_update;       /**< */
   float  min_ener;          /**< */
   drft_lookup fft_lookup;   /**< */

} SpeexDenoiseState;

/** Creates a new denoising state */
SpeexDenoiseState *speex_denoise_state_init(int frame_size);

/** Destroys a denoising state */
void speex_denoise_state_destroy(SpeexDenoiseState *st);

/** Denoise a frame */
int speex_denoise(SpeexDenoiseState *st, float *x);

#ifdef __cplusplus
}
#endif
