/* Copyright (C) 2002 Jean-Marc Valin*/
/**
  @file speex_callbacks.h
  @brief Describes callback handling and in-band signalling
*/
/*
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:
   
   - Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
   
   - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
   
   - Neither the name of the Xiph.org Foundation nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.
   
   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE FOUNDATION OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef SPEEX_CALLBACKS_H
#define SPEEX_CALLBACKS_H

#include "speex.h"

#ifdef __cplusplus
extern "C" {
#endif

#define SPEEX_MAX_CALLBACKS 16

/* Describes all the in-band requests */

/*These are 1-bit requests*/
#define SPEEX_INBAND_ENH_REQUEST         0
#define SPEEX_INBAND_VBR_REQUEST         1

/*These are 4-bit requests*/
#define SPEEX_INBAND_MODE_REQUEST        2
#define SPEEX_INBAND_LOW_MODE_REQUEST    3
#define SPEEX_INBAND_HIGH_MODE_REQUEST   4
#define SPEEX_INBAND_VBR_QUALITY_REQUEST 5
#define SPEEX_INBAND_ACKNOWLEDGE_REQUEST 6

/*These are 8-bit requests*/
/** Send a character in-band */
#define SPEEX_INBAND_CHAR                8

#define SPEEX_INBAND_MAX_BITRATE         10

#define SPEEX_INBAND_ACKNOWLEDGE         12

/** Callback function type */
typedef int (*speex_callback_func)(SpeexBits *bits, void *state, void *data);

typedef struct SpeexCallback {
   int callback_id;
   speex_callback_func func;
   void *data;
} SpeexCallback;

int speex_inband_handler(SpeexBits *bits, SpeexCallback *callback_list, void *state);

int speex_std_mode_request_handler(SpeexBits *bits, void *state, void *data);

int speex_std_high_mode_request_handler(SpeexBits *bits, void *state, void *data);

int speex_std_char_handler(SpeexBits *bits, void *state, void *data);


int speex_default_user_handler(SpeexBits *bits, void *state, void *data);

#ifdef __cplusplus
}
#endif


#endif
