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

#include "speex_bits.h"

#ifdef __cplusplus
extern "C" {
#endif

#define SPEEX_SET_PF 0

struct SpeexMode;

typedef void *(*encoder_init_func)(struct SpeexMode *mode);
typedef void (*encoder_destroy_func)(void *st);
typedef void (*encode_func)(void *state, float *in, SpeexBits *bits);
typedef void *(*decoder_init_func)(struct SpeexMode *mode);
typedef void (*decoder_destroy_func)(void *st);
typedef void (*decode_func)(void *state, SpeexBits *bits, float *out, int lost);
typedef void (*ctl_func)(void *state, int request, void *ptr);

/** Struct defining a Speex mode */ 
typedef struct SpeexMode {
   /** Pointer to the low-level mode data */
   void *mode;

   /** Pointer to encoder initialization function */
   encoder_init_func enc_init;

   /** Pointer to encoder destruction function */
   encoder_destroy_func enc_destroy;

   /** Pointer to frame encoding function */
   encode_func enc;

   /** Pointer to decoder initialization function */
   decoder_init_func dec_init;

   /** Pointer to decoder destruction function */
   decoder_destroy_func dec_destroy;

   /** Pointer to frame decoding function */
   decode_func dec;

   /** ioctl-like requests for codec state */
   ctl_func ctl;

   /** Frame size used for the current mode */
   int frameSize;

} SpeexMode;

/**Returns a handle to a newly created Speex encoder state structure. For now, the 
   "mode" arguent can be &nb_mode or &wb_mode . In the future, more modes may be 
   added. Note that for now if you have more than one channels to encode, you need 
   one state per channel.*/
void *speex_encoder_init(SpeexMode *mode);

/** Frees all resources associated to an existing Speex encoder state. */
void speex_encoder_destroy(void *state);

/** Uses an existing encoder state to encode one frame of speech pointed to by
    "in". The encoded bit-stream is saved in "bits".*/
void speex_encode(void *state, float *in, SpeexBits *bits);

/** Returns a handle to a newly created decoder state structure. For now, the mode
    arguent can be &nb_mode or &wb_mode . In the future, more modes may be added. 
    Note that for now if you have more than one channels to decode, you need one 
    state per channel. */ 
void *speex_decoder_init(SpeexMode *mode);

/** Frees all resources associated to an existing decoder state. */
void speex_decoder_destroy(void *state);

/** Uses an existing decoder state to decode one frame of speech from bit-stream 
    bits. The output speech is saved written to out. */
void speex_decode(void *state, SpeexBits *bits, float *out, int lost);


void speex_ctl(void *state, int request, void *ptr);

/** Default narrowband mode */
extern SpeexMode speex_nb_mode;

/** Default wideband mode */
extern SpeexMode speex_wb_mode;


#ifdef __cplusplus
}
#endif


#endif
