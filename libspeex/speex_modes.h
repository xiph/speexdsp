/* Copyright (C) 2002 Jean-Marc Valin 
   File: modes.h

   Describes the different modes of the codec

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.
   
   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

#ifndef SPEEX_MODES_H
#define SPEEX_MODES_H

#ifdef __cplusplus
extern "C" {
#endif

#include "speex_bits.h"

struct SpeexMode;

typedef void *(*encoder_init_func)(struct SpeexMode *mode);
typedef void (*encoder_destroy_func)(void *st);
typedef void (*encode_func)(void *state, float *in, FrameBits *bits);
typedef void *(*decoder_init_func)(struct SpeexMode *mode);
typedef void (*decoder_destroy_func)(void *st);
typedef void (*decode_func)(void *state, FrameBits *bits, float *out);

typedef struct SpeexMode {
   void *mode;
   encoder_init_func enc_init;
   encoder_destroy_func enc_destroy;
   encode_func enc;
   decoder_init_func dec_init;
   decoder_destroy_func dec_destroy;
   decode_func dec;
   int frameSize;
} SpeexMode;

void *encoder_init(SpeexMode *mode);
void encoder_destroy(void *state);
void encode(void *state, float *in, FrameBits *bits);
void *decoder_init(SpeexMode *mode);
void decoder_destroy(void *state);
void decode(void *state, FrameBits *bits, float *out);

extern SpeexMode nb_mode;
extern SpeexMode wb_mode;


#ifdef __cplusplus
}
#endif


#endif
