/* Copyright (C) 2002 Jean-Marc Valin 
   File: sb_celp.h

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

#ifndef SB_CELP_H
#define SB_CELP_H

#include "modes.h"
#include "bits.h"
#include "speex.h"

/**Structure representing the full state of the encoder*/
typedef struct SBEncState {
   EncState st_low;
   int    full_frame_size;     /* Length of full-band frames*/
   int    frame_size;          /* Length of high-band frames*/
   int    subframeSize;        /* Length of high-band sub-frames*/
   int    nbSubframes;         /* Number of high-band sub-frames*/
   int    windowSize;          /* Length of high-band LPC window*/
   int    lpcSize;             /* Order of high-band LPC analysis */
   int    first;               /* First frame? */
   float  lag_factor;          /* Lag-windowing control parameter */
   float  lpc_floor;           /* Controls LPC analysis noise floor */
   float  gamma1;              /* Perceptual weighting coef 1 */
   float  gamma2;              /* Perceptual weighting coef 2 */

   float *stack;               /* Temporary allocation stack */
   float *x0, *x0d, *x1, *x1d; /* QMF filter signals*/
   float *high;                /* High-band signal (buffer) */
   float *y0, *y1;             /* QMF synthesis signals */
   float *h0_mem, *h1_mem, *g0_mem, *g1_mem; /* QMF memories */

   float *excBuf;              /* High-band excitation */
   float *exc;                 /* High-band excitation (for QMF only)*/
   float *buf;                 /* Temporary buffer */
   float *res;                 /* Zero-input response (ringing) */
   float *sw;                  /* Perceptually weighted signal */
   float *target;              /* Weighted target signal (analysis by synthesis) */
   float *window;              /* LPC analysis window */
   float *lagWindow;           /* Auto-correlation window */
   float *autocorr;            /* Auto-correlation (for LPC analysis) */
   float *rc;                  /* Reflection coefficients (unused) */
   float *lpc;                 /* LPC coefficients */
   float *lsp;                 /* LSP coefficients */
   float *qlsp;                /* Quantized LSPs */
   float *old_lsp;             /* LSPs of previous frame */
   float *old_qlsp;            /* Quantized LSPs of previous frame */
   float *interp_lsp;          /* Interpolated LSPs for current sub-frame */
   float *interp_qlsp;         /* Interpolated quantized LSPs for current sub-frame */
   float *interp_lpc;          /* Interpolated LPCs for current sub-frame */
   float *interp_qlpc;         /* Interpolated quantized LPCs for current sub-frame */
   float *bw_lpc1;             /* Bandwidth-expanded version of LPCs (#1) */
   float *bw_lpc2;             /* Bandwidth-expanded version of LPCs (#2) */

   float *mem_sp;              /* Synthesis signal memory */
   float *mem_sp2;
   float *mem_sw;              /* Perceptual signal memory */
} SBEncState;


/**Structure representing the full state of the decoder*/
typedef struct SBDecState {
   DecState st_low;
   int    full_frame_size;
   int    frame_size;
   int    subframeSize;
   int    nbSubframes;
   int    lpcSize;
   int    first;

   float *stack;
   float *x0, *x0d, *x1, *x1d;
   float *high;
   float *y0, *y1;
   float *h0_mem, *h1_mem, *g0_mem, *g1_mem;

   float *exc;
   float *qlsp;
   float *old_qlsp;
   float *interp_qlsp;
   float *interp_qlpc;

   float *mem_sp;
   float *mem_sw;
} SBDecState;


/**Initializes encoder state*/
void sb_encoder_init(SBEncState *st, SpeexMode *mode);

/**De-allocates encoder state resources*/
void sb_encoder_destroy(SBEncState *st);

/**Encodes one frame*/
void sb_encode(SBEncState *st, float *in, FrameBits *bits);


/**Initializes decoder state*/
void sb_decoder_init(SBDecState *st, SpeexMode *mode);

/**De-allocates decoder state resources*/
void sb_decoder_destroy(SBDecState *st);

/**Decodes one frame*/
void sb_decode(SBDecState *st, FrameBits *bits, float *out);


#endif
