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

#include "bits.h"

/* Quantizes LSPs */
typedef void (*lsp_quant_func)(float *, float *, int, FrameBits *);

/* Decodes quantized LSPs */
typedef void (*lsp_unquant_func)(float *, int, FrameBits *);


/*Long-term predictor quantization*/
typedef void (*ltp_quant_func)(float *, float *, float *, 
                                float *, float *, void *, int, int, 
                                int, int, FrameBits*, float *);

/*Long-term un-quantize*/
typedef void (*ltp_unquant_func)(float *, int, int, void *, int, FrameBits*, float*);


typedef void (*innovation_quant_func)(float *, float *, float *, float *, void *, int, int, 
                                      float *, FrameBits *, float *);

typedef void (*innovation_unquant_func)(float *, void *, int, FrameBits*, float *);

/*Struct defining the encoding/decoding mode*/
typedef struct SpeexMode {
   int     frameSize;
   int     subframeSize;
   int     windowSize;
   int     lpcSize;
   int     bufSize;
   int     pitchStart;
   int     pitchEnd;
   float   gamma1;
   float   gamma2;
   
   /*LSP functions*/
   lsp_quant_func    lsp_quant;
   lsp_unquant_func  lsp_unquant;

   /*Lont-term predictor functions*/
   ltp_quant_func    ltp_quant;
   ltp_unquant_func  ltp_unquant;
   void             *ltp_params;
   /*Quantization of innovation */
   innovation_quant_func innovation_quant;
   innovation_unquant_func innovation_unquant;
   void             *innovation_params;

} SpeexMode;

extern SpeexMode nb_mode;

#endif
