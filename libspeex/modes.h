/* Copyright (C) 2002 Jean-Marc Valin */
/**
   @file modes.h
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

#ifndef MODES_H
#define MODES_H

#include "speex.h"
#include "speex_bits.h"


#define NB_SUBMODES 16
#define NB_SUBMODE_BITS 4

#define SB_SUBMODES 8
#define SB_SUBMODE_BITS 3


/** Quantizes LSPs */
typedef void (*lsp_quant_func)(float *, float *, int, SpeexBits *);

/** Decodes quantized LSPs */
typedef void (*lsp_unquant_func)(float *, int, SpeexBits *);


/** Long-term predictor quantization */
typedef int (*ltp_quant_func)(float *, float *, float *, float *, 
                              float *, float *, void *, int, int, float, 
                              int, int, SpeexBits*, float *, float *, int);

/** Long-term un-quantize */
typedef void (*ltp_unquant_func)(float *, int, int, float, void *, int, int *,
                                 float *, SpeexBits*, float*, int);


/** Innovation quantization function */
typedef void (*innovation_quant_func)(float *, float *, float *, float *, void *, int, int, 
                                      float *, SpeexBits *, float *, int);

/** Innovation unquantization function */
typedef void (*innovation_unquant_func)(float *, void *, int, SpeexBits*, float *);

/** Description of a Speex sub-mode (wither narrowband or wideband */
typedef struct SpeexSubmode {
   int     lbr_pitch; /**< Set to -1 for "normal" modes, otherwise encode pitch using a global pitch and allowing a +- lbr_pitch variation (for low not-rates)*/
   int     forced_pitch_gain; /**< Use the same (forced) pitch gain for all sub-frames */
   int     have_subframe_gain; /**< Number of bits to use as sub-frame innovation gain */
   int     double_codebook; /**< Apply innovation quantization twice for higher quality (and higher bit-rate)*/
   /*LSP functions*/
   lsp_quant_func    lsp_quant; /**< LSP quantization function */
   lsp_unquant_func  lsp_unquant; /**< LSP unquantization function */

   /*Lont-term predictor functions*/
   ltp_quant_func    ltp_quant; /**< Long-term predictor (pitch) quantizer */
   ltp_unquant_func  ltp_unquant; /**< Long-term predictor (pitch) un-quantizer */
   void             *ltp_params; /**< Pitch parameters (options) */

   /*Quantization of innovation*/
   innovation_quant_func innovation_quant; /**< Innovation quantization */
   innovation_unquant_func innovation_unquant; /**< Innovation un-quantization */
   void             *innovation_params; /**< Innovation quantization parameters*/

   /*Synthesis filter enhancement*/
   float             lpc_enh_k1; /**< Enhancer constant */
   float             lpc_enh_k2; /**< Enhancer constant */
   float             comb_gain;  /**< Gain of enhancer comb filter */

   int               bits_per_frame; /**< Number of bits per frame after encoding*/
} SpeexSubmode;

/** Struct defining the encoding/decoding mode*/
typedef struct SpeexNBMode {
   int     frameSize; /**< Size of frames used for encoding */
   int     subframeSize; /**< Size of sub-frames used for encoding */
   int     lpcSize; /**< Order of LPC filter */
   int     bufSize; /**< Size of signal buffer to use in encoder */
   int     pitchStart; /**< Smallest pitch value allowed */
   int     pitchEnd; /**< Largest pitch value allowed */

   float   gamma1; /**< Perceptual filter parameter #1 */
   float   gamma2; /**< Perceptual filter parameter #2 */
   float   lag_factor; /**< Lag-windowing parameter */
   float   lpc_floor; /**< Noise floor for LPC analysis */
   float   preemph; /**< Pre-emphasis */

   SpeexSubmode *submodes[NB_SUBMODES]; /**< Sub-mode data for the mode */
   int     defaultSubmode; /**< Default sub-mode to use when encoding */

} SpeexNBMode;


/** Struct defining the encoding/decoding mode for SB-CELP (wideband) */
typedef struct SpeexSBMode {
   SpeexMode *nb_mode; /**< Embedded narrowband mode */
   int     frameSize; /**< Size of frames used for encoding */
   int     subframeSize; /**< Size of sub-frames used for encoding */
   int     lpcSize; /**< Order of LPC filter */
   int     bufSize; /**< Signal buffer size in encoder */
   float   gamma1; /**< Perceptual filter parameter #1 */
   float   gamma2; /**< Perceptual filter parameter #1 */
   float   lag_factor; /**< Lag-windowing parameter */
   float   lpc_floor; /**< Noise floor for LPC analysis */
   float   preemph; /**< Pre-emphasis */

   SpeexSubmode *submodes[SB_SUBMODES]; /**< Sub-mode data for the mode */
   int     defaultSubmode; /**< Default sub-mode to use when encoding */

} SpeexSBMode;


#endif
