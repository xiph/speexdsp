/* Copyright (C) 2002 Jean-Marc Valin 
   File: speex.h

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

/**Structure representing the full state of the encoder*/
typedef struct EncState {
   int    first;          /* Is this the first frame? */
   int    frameSize;      /* Size of frames */
   int    subframeSize;   /* Size of sub-frames */
   int    nbSubframes;    /* Number of sub-frames */
   int    windowSize;     /* Analysis (LPC) window length */
   int    lpcSize;        /* LPC order */
   int    bufSize;        /* Buffer size */
   float  gamma;          /* Perceptual filter: A(z/gamma) */
   float *inBuf;          /* Input buffer (original signal) */
   float *frame;          /* Start of original frame */
   float *wBuf;           /* "Weighted" buffer */
   float *wframe;         /* Start of "weighted" frame */
   float *excBuf;         /* Excitation buffer */
   float *exc_frame;      /* Start of excitation frame */
   float *tBuf;           /* "weighted target" buffer */
   float *tframe;         /* Start of "weighted target" frame */
   float *window;         /* Temporary (Hanning) window */
   float *buf2;           /* 2nd temporary buffer */
   float *autocorr;       /* auto-correlation */
   float *lagWindow;      /* Window applied to auto-correlation */
   float *lpc;            /* LPCs for current frame */
   float *lsp;            /* LSPs for current frame */
   float *old_lsp;        /* LSPs for previous frame */
   float *interp_lsp;     /* Interpolated LSPs */
   float *interp_lpc;     /* Interpolated LPCs */
   float *bw_lpc;         /* LPCs after bandwidth expansion for perceptual weighting*/
   float *rc;             /* Reflection coefficients */
} EncState;

typedef struct DecState {
} DecState;

/**Initializes encoder state*/
void encoder_init(EncState *st);
/**De-allocates encoder state resources*/
void encoder_destroy(EncState *st);
/**Encodes one frame*/
void encode(EncState *st, float *in, int *outSize, void *bits);

void decoder_init(DecState *st);
void decoder_destroy(DecState *st);
void decode(DecState *st, float *bits, float *out);



#endif
