/* Copyright (C) 2002 Jean-Marc Valin*/
/**
  @file speex.h
  @brief Describes the different modes of the codec
*/
/*
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

#ifndef SPEEX_H
#define SPEEX_H

#include "speex_bits.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Values allowed for *ctl() requests */

/** Set enhancement on/off (decoder only) */
#define SPEEX_SET_ENH 0
/** Get enhancement state (decoder only) */
#define SPEEX_GET_ENH 1

/*Would be SPEEX_SET_FRAME_SIZE, but it's (currently) invalid*/
/** Obtain frame size used by encoder/decoder */
#define SPEEX_GET_FRAME_SIZE 3

/** Set quality value */
#define SPEEX_SET_QUALITY 4
/** Get current quality setting */
#define SPEEX_GET_QUALITY 5

/** Set sub-mode to use */
#define SPEEX_SET_MODE 6
/** Get current sub-mode in use */
#define SPEEX_GET_MODE 7

/** Set low-band sub-mode to use (wideband only)*/
#define SPEEX_SET_LOW_MODE 8
/** Get current low-band mode in use (wideband only)*/
#define SPEEX_GET_LOW_MODE 9

/** Set high-band sub-mode to use (wideband only)*/
#define SPEEX_SET_HIGH_MODE 10
/** Get current high-band mode in use (wideband only)*/
#define SPEEX_GET_HIGH_MODE 11

/** Set VBR on (1) or off (0) */
#define SPEEX_SET_VBR 12
/** Get VBR status (1 for on, 0 for off) */
#define SPEEX_GET_VBR 13

/** Set quality value for VBR encoding (0-10) */
#define SPEEX_SET_VBR_QUALITY 14
/** Get current quality value for VBR encoding (0-10) */
#define SPEEX_GET_VBR_QUALITY 15

/** Set complexity of the encoder (0-10) */
#define SPEEX_SET_COMPLEXITY 16
/** Get current complexity of the encoder (0-10) */
#define SPEEX_GET_COMPLEXITY 17

/*Would be SPEEX_SET_FRAME_SIZE, but it's (currently) invalid*/
/** Get current bit-rate used by the encoder or decoder */
#define SPEEX_GET_BITRATE 19


/* Preserving compatibility:*/
/** Equivalent to SPEEX_SET_ENH */
#define SPEEX_SET_PF 0
/** Equivalent to SPEEX_GET_ENH */
#define SPEEX_GET_PF 1


/* Values allowed for mode queries */
/** Query the frame size of a mode */
#define SPEEX_MODE_FRAME_SIZE 0

/** Query the size of an encoded frame for a particular sub-mode */
#define SPEEX_SUBMODE_BITS_PER_FRAME 1


/** Number of defined modes in Speex */
#define SPEEX_NB_MODES 2

struct SpeexMode;


/* Prototypes for mode function pointers */

/** Encoder state initialization function */
typedef void *(*encoder_init_func)(struct SpeexMode *mode);

/** Encoder state destruction function */
typedef void (*encoder_destroy_func)(void *st);

/** Main encoding function */
typedef void (*encode_func)(void *state, float *in, SpeexBits *bits);

/** Function for controlling the encoder options */
typedef void (*encoder_ctl_func)(void *state, int request, void *ptr);

/** Decoder state initialization function */
typedef void *(*decoder_init_func)(struct SpeexMode *mode);

/** Decoder state destruction function */
typedef void (*decoder_destroy_func)(void *st);

/** Main decoding function */
typedef int  (*decode_func)(void *state, SpeexBits *bits, float *out);

/** Function for controlling the decoder options */
typedef void (*decoder_ctl_func)(void *state, int request, void *ptr);


/** Query function for a mode */
typedef void (*mode_query_func)(void *mode, int request, void *ptr);

/** Struct defining a Speex mode */ 
typedef struct SpeexMode {
   /** Pointer to the low-level mode data */
   void *mode;

   /** Pointer to the mode query function */
   mode_query_func query;
   
   /** The name of the mode (you should not rely on this to identify the mode)*/
   char *modeName;

   /**ID of the mode*/
   int modeID;

   /**Version number of the bitstream (incremented every time we break
    bitstream compatibility*/
   int bitstream_version;

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

   /** ioctl-like requests for encoder */
   encoder_ctl_func enc_ctl;

   /** ioctl-like requests for decoder */
   decoder_ctl_func dec_ctl;

} SpeexMode;

/**
 * Returns a handle to a newly created Speex encoder state structure. For now, 
 * the "mode" arguent can be &nb_mode or &wb_mode . In the future, more modes 
 * may be added. Note that for now if you have more than one channels to 
 * encode, you need one state per channel.
 *
 * @param mode The mode to use (either speex_nb_mode or speex_wb.mode) 
 * @return A newly created encoder
 */
void *speex_encoder_init(SpeexMode *mode);

/** Frees all resources associated to an existing Speex encoder state. 
 * @param state Encoder state to be destroyed */
void speex_encoder_destroy(void *state);

/** Uses an existing encoder state to encode one frame of speech pointed to by
    "in". The encoded bit-stream is saved in "bits".
 @param state Encoder state
 @param in Frame that will be encoded with a +-2^16 range
 @param bits Bit-stream where the data will be written
 */
void speex_encode(void *state, float *in, SpeexBits *bits);

/** Used like the ioctl function to control the encoder parameters
 *
 * @param state Encoder state
 * @param request ioctl-type request (one of the SPEEX_* macros)
 * @param ptr Data exchanged to-from function
 */
void speex_encoder_ctl(void *state, int request, void *ptr);


/** Returns a handle to a newly created decoder state structure. For now, 
 * the mode arguent can be &nb_mode or &wb_mode . In the future, more modes
 * may be added.  Note that for now if you have more than one channels to
 * decode, you need one state per channel.
 *
 * @param mode Speex mode (one of speex_nb_mode or speex_wb_mode)
 * @return A newly created decoder state
 */ 
void *speex_decoder_init(SpeexMode *mode);

/** Frees all resources associated to an existing decoder state.
 *
 * @param state State to be destroyed
 */
void speex_decoder_destroy(void *state);

/** Uses an existing decoder state to decode one frame of speech from
 * bit-stream bits. The output speech is saved written to out.
 *
 * @param state Decoder state
 * @param bits Bit-stream from which to decode the frame (NULL if the packet was lost)
 * @param out Where to write the decoded frame
 * @return return status (0 for no error, -1 for end of stream, -2 other)
 */
int speex_decode(void *state, SpeexBits *bits, float *out);

/** Used like the ioctl function to control the encoder parameters
 *
 * @param state Decoder state
 * @param request ioctl-type request (one of the SPEEX_* macros)
 * @param ptr Data exchanged to-from function
 */
void speex_decoder_ctl(void *state, int request, void *ptr);


/** Query function for mode information
 *
 * @param mode Speex mode
 * @param request ioctl-type request (one of the SPEEX_* macros)
 * @param ptr Data exchanged to-from function
 */
void speex_mode_query(SpeexMode *mode, int request, void *ptr);


/** Default narrowband mode */
extern SpeexMode speex_nb_mode;

/** Default wideband mode */
extern SpeexMode speex_wb_mode;

/** List of all modes availavle */
extern SpeexMode *speex_mode_list[SPEEX_NB_MODES];

#ifdef __cplusplus
}
#endif


#endif
