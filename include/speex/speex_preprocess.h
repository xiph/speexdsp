/* Copyright (C) 2003 Epic Games
   Written by Jean-Marc Valin */
/**
   @file speex_preprocess.h
   @brief Speex preprocessor
*/
/*
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

#ifndef SPEEX_PREPROCESS_H
#define SPEEX_PREPROCESS_H

#include "speex/speex_types.h"

#ifdef __cplusplus
extern "C" {
#endif

struct SpeexPreprocessState_;

typedef struct SpeexPreprocessState_ SpeexPreprocessState;


/** Creates a new preprocessing state */
SpeexPreprocessState *speex_preprocess_state_init(int frame_size, int sampling_rate);

/** Destroys a denoising state */
void speex_preprocess_state_destroy(SpeexPreprocessState *st);

/** Preprocess a frame */
int speex_preprocess(SpeexPreprocessState *st, spx_int16_t *x, spx_int32_t *echo);

/** Preprocess a frame */
void speex_preprocess_estimate_update(SpeexPreprocessState *st, spx_int16_t *x, spx_int32_t *echo);

/** Used like the ioctl function to control the preprocessor parameters */
int speex_preprocess_ctl(SpeexPreprocessState *st, int request, void *ptr);



/** Set preprocessor denoiser state */
#define SPEEX_PREPROCESS_SET_DENOISE 0
/** Get preprocessor denoiser state */
#define SPEEX_PREPROCESS_GET_DENOISE 1

/** Set preprocessor Automatic Gain Control state */
#define SPEEX_PREPROCESS_SET_AGC 2
/** Get preprocessor Automatic Gain Control state */
#define SPEEX_PREPROCESS_GET_AGC 3

/** Set preprocessor Voice Activity Detection state */
#define SPEEX_PREPROCESS_SET_VAD 4
/** Get preprocessor Voice Activity Detection state */
#define SPEEX_PREPROCESS_GET_VAD 5

/** Set preprocessor Automatic Gain Control level */
#define SPEEX_PREPROCESS_SET_AGC_LEVEL 6
/** Get preprocessor Automatic Gain Control level */
#define SPEEX_PREPROCESS_GET_AGC_LEVEL 7

/** Set preprocessor dereverb state */
#define SPEEX_PREPROCESS_SET_DEREVERB 8
/** Get preprocessor dereverb state */
#define SPEEX_PREPROCESS_GET_DEREVERB 9

/** Set preprocessor dereverb level */
#define SPEEX_PREPROCESS_SET_DEREVERB_LEVEL 10
/** Get preprocessor dereverb level */
#define SPEEX_PREPROCESS_GET_DEREVERB_LEVEL 11

/** Set preprocessor dereverb decay */
#define SPEEX_PREPROCESS_SET_DEREVERB_DECAY 12
/** Get preprocessor dereverb decay */
#define SPEEX_PREPROCESS_GET_DEREVERB_DECAY 13

#define SPEEX_PREPROCESS_SET_PROB_START 14
#define SPEEX_PREPROCESS_GET_PROB_START 15

#define SPEEX_PREPROCESS_SET_PROB_CONTINUE 16
#define SPEEX_PREPROCESS_GET_PROB_CONTINUE 17

#ifdef __cplusplus
}
#endif

#endif
