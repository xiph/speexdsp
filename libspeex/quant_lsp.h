/* Copyright (C) 2002 Jean-Marc Valin 
   File: quant_lsp.h
   LSP vector quantization

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

#ifndef QUANT_LSP_H
#define QUANT_LSP_H

#include "speex_bits.h"

#define MAX_LSP_SIZE 20

#define NB_CDBK_SIZE 64
#define NB_CDBK_SIZE_LOW1 64
#define NB_CDBK_SIZE_LOW2 64
#define NB_CDBK_SIZE_HIGH1 64
#define NB_CDBK_SIZE_HIGH2 64

/*Narrowband codebooks*/
extern float cdbk_nb[];
extern float cdbk_nb_low1[];
extern float cdbk_nb_low2[];
extern float cdbk_nb_high1[];
extern float cdbk_nb_high2[];

/* Quantizes narrowband LSPs with 30 bits */
void lsp_quant_nb(float *lsp, float *qlsp, int order, FrameBits *bits);

/* Decodes quantized narrowband LSPs */
void lsp_unquant_nb(float *lsp, int order, FrameBits *bits);

/* Quantizes wideband LSPs with 50 bits */
void lsp_quant_wb(float *lsp, float *qlsp, int order, FrameBits *bits);

/* Decodes quantized wideband LSPs */
void lsp_unquant_wb(float *lsp, int order, FrameBits *bits);

/* Quantizes high-band LSPs with 12 bits */
void lsp_quant_high(float *lsp, float *qlsp, int order, FrameBits *bits);

/* Decodes high-band LSPs */
void lsp_unquant_high(float *lsp, int order, FrameBits *bits);

#endif
