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

#ifndef MODES_H
#define MODES_H

#include "speex.h"
#include "speex_bits.h"


#define NB_SUBMODES 16
#define NB_SUBMODE_BITS 4

#define SB_SUBMODES 8
#define SB_SUBMODE_BITS 3


/* Quantizes LSPs */
typedef void (*lsp_quant_func)(float *, float *, int, SpeexBits *);

/* Decodes quantized LSPs */
typedef void (*lsp_unquant_func)(float *, int, SpeexBits *);


/*Long-term predictor quantization*/
typedef int (*ltp_quant_func)(float *, float *, float *, float *, 
                                float *, float *, void *, int, int, 
                                int, int, SpeexBits*, float *, float *, int);

/*Long-term un-quantize*/
typedef void (*ltp_unquant_func)(float *, int, int, void *, int, int *, float *, SpeexBits*, float*, int);


typedef void (*innovation_quant_func)(float *, float *, float *, float *, void *, int, int, 
                                      float *, SpeexBits *, float *, int);

typedef void (*innovation_unquant_func)(float *, void *, int, SpeexBits*, float *);

typedef struct SpeexSubmode {
   int     lbr_pitch;
   int     have_subframe_gain;
   /*LSP functions*/
   lsp_quant_func    lsp_quant;
   lsp_unquant_func  lsp_unquant;

   /*Lont-term predictor functions*/
   ltp_quant_func    ltp_quant;
   ltp_unquant_func  ltp_unquant;
   void             *ltp_params;

   /*Quantization of innovation*/
   innovation_quant_func innovation_quant;
   innovation_unquant_func innovation_unquant;
   void             *innovation_params;

   /*Synthesis filter enhancement*/
   float lpc_enh_k1, lpc_enh_k2, comb_gain;

   int               bits_per_frame;
} SpeexSubmode;

/*Struct defining the encoding/decoding mode*/
typedef struct SpeexNBMode {
   int     frameSize;
   int     subframeSize;
   int     lpcSize;
   int     bufSize;
   int     pitchStart;
   int     pitchEnd;

   float   gamma1;
   float   gamma2;
   float   lag_factor;
   float   lpc_floor;
   float   preemph;

   SpeexSubmode *submodes[NB_SUBMODES];
   int     defaultSubmode;

} SpeexNBMode;


/*Struct defining the encoding/decoding mode*/
typedef struct SpeexSBMode {
   SpeexMode *nb_mode;
   int     frameSize;
   int     subframeSize;
   int     lpcSize;
   int     bufSize;
   float   gamma1;
   float   gamma2;
   float   lag_factor;
   float   lpc_floor;
   float   preemph;

   SpeexSubmode *submodes[SB_SUBMODES];
   int     defaultSubmode;

} SpeexSBMode;


extern SpeexSBMode sb_wb_mode;

#endif
