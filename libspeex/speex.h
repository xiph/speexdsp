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

/** Struct defining a Speex mode */ 
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

   /** Creates an encoder state ("object") from a mode */ 
void *encoder_init(SpeexMode *mode);

   /** Destroy a Speex encoder state */
void encoder_destroy(void *state);

   /** Encode a frame */
void encode(void *state, float *in, FrameBits *bits);

   /** Creates a decoder state ("object") from a mode */ 
void *decoder_init(SpeexMode *mode);

   /** Destroy a Speex decoder state */
void decoder_destroy(void *state);

   /** Decode a frame */
void decode(void *state, FrameBits *bits, float *out);

   /** Default narrowband mode */
extern SpeexMode speex_nb_mode;

   /** Default wideband mode */
extern SpeexMode speex_wb_mode;


#ifdef __cplusplus
}
#endif


#endif
