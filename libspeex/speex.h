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

#include "modes.h"
#include "bits.h"

/**Structure representing the full state of the encoder*/
typedef struct EncState {
   int    first;          /* Is this the first frame? */
   int    frameSize;      /* Size of frames */
   int    subframeSize;   /* Size of sub-frames */
   int    nbSubframes;    /* Number of sub-frames */
   int    windowSize;     /* Analysis (LPC) window length */
   int    lpcSize;        /* LPC order */
   int    bufSize;        /* Buffer size */
   int    min_pitch;      /* Minimum pitch value allowed */
   int    max_pitch;      /* Maximum pitch value allowed */
   int    ol_pitch;       /* Open-loop pitch */
   int    ol_voiced;      /* Open-loop voiced/non-voiced decision */
   float  gamma1;         /* Perceptual filter: A(z/gamma1) */
   float  gamma2;         /* Perceptual filter: A(z/gamma2) */
   float  lag_factor;     /* Lag windowing gaussian width */
   float  lpc_floor;      /* Noise floor multiplier for A[0] in LPC analysis*/
   float *stack;          /* Pseudo-stack allocation */
   int    os_fact;        /* Over-sampling factor for fractional pitch */
   int    os_filt_ord2;   /* Over-sampling filter size for fractional pitch */
   float *os_exc;         /* Over-sampled excitation for fractional pitch */
   float *os_filt;        /* Over-sampling filter for fractional pitch */
   float *inBuf;          /* Input buffer (original signal) */
   float *frame;          /* Start of original frame */
   float *excBuf;         /* Excitation buffer */
   float *exc;            /* Start of excitation frame */
   float *swBuf;          /* Weighted signal buffer */
   float *sw;             /* Start of weighted signal frame */
   float *window;         /* Temporary (Hanning) window */
   float *buf2;           /* 2nd temporary buffer */
   float *autocorr;       /* auto-correlation */
   float *lagWindow;      /* Window applied to auto-correlation */
   float *lpc;            /* LPCs for current frame */
   float *lsp;            /* LSPs for current frame */
   float *qlsp;           /* Quantized LSPs for current frame */
   float *old_lsp;        /* LSPs for previous frame */
   float *old_qlsp;       /* Quantized LSPs for previous frame */
   float *interp_lsp;     /* Interpolated LSPs */
   float *interp_qlsp;    /* Interpolated quantized LSPs */
   float *interp_lpc;     /* Interpolated LPCs */
   float *interp_qlpc;    /* Interpolated quantized LPCs */
   float *bw_lpc1;        /* LPCs after bandwidth expansion by gamma1 for perceptual weighting*/
   float *bw_lpc2;        /* LPCs after bandwidth expansion by gamma2 for perceptual weighting*/
   float *rc;             /* Reflection coefficients */
   float *mem_sp, *mem_sw;
   float *dmem1, *dmem2;

   lsp_quant_func    lsp_quant;
   ltp_quant_func    ltp_quant;
   void             *ltp_params;
   innovation_quant_func innovation_quant;
   void             *innovation_params;
} EncState;

typedef struct DecState {
   int    first;          /* Is this the first frame? */
   int    frameSize;      /* Size of frames */
   int    subframeSize;   /* Size of sub-frames */
   int    nbSubframes;    /* Number of sub-frames */
   int    windowSize;     /* Analysis (LPC) window length */
   int    lpcSize;        /* LPC order */
   int    bufSize;        /* Buffer size */
   int    min_pitch;      /* Minimum pitch value allowed */
   int    max_pitch;      /* Maximum pitch value allowed */
   float  gamma1;         /* Perceptual filter: A(z/gamma1) */
   float  gamma2;         /* Perceptual filter: A(z/gamma2) */
   float *stack;
   float *inBuf;          /* Input buffer (original signal) */
   float *frame;          /* Start of original frame */
   float *excBuf;         /* Excitation buffer */
   float *exc;            /* Start of excitation frame */
   float *qlsp;           /* Quantized LSPs for current frame */
   float *old_qlsp;       /* Quantized LSPs for previous frame */
   float *interp_qlsp;    /* Interpolated quantized LSPs */
   float *interp_qlpc;    /* Interpolated quantized LPCs */
   float *mem_sp;

   lsp_unquant_func  lsp_unquant;
   ltp_unquant_func  ltp_unquant;
   void             *ltp_params;
   innovation_unquant_func innovation_unquant;
   void             *innovation_params;
} DecState;

/**Initializes encoder state*/
void encoder_init(EncState *st, SpeexMode *mode);

/**De-allocates encoder state resources*/
void encoder_destroy(EncState *st);

/**Encodes one frame*/
void encode(EncState *st, float *in, FrameBits *bits);


/**Initializes decoder state*/
void decoder_init(DecState *st, SpeexMode *mode);

/**De-allocates decoder state resources*/
void decoder_destroy(DecState *st);

/**Decodes one frame*/
void decode(DecState *st, FrameBits *bits, float *out);



#endif
