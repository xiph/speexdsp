/* Copyright (C) 2002 Jean-Marc Valin 
   File: quant_lsp.h
   LSP vector quantization
*/

#ifndef QUANT_LSP_H
#define QUANT_LSP_H

#define MAX_LSP_SIZE 20

#define NB_CDBK_SIZE 64
#define NB_CDBK_SIZE_LOW1 64
#define NB_CDBK_SIZE_LOW2 64
#define NB_CDBK_SIZE_HIGH1 64
#define NB_CDBK_SIZE_HIGH2 64

extern float cdbk_nb[];
extern float cdbk_nb_low1[];
extern float cdbk_nb_low2[];
extern float cdbk_nb_high1[];
extern float cdbk_nb_high2[];


unsigned int lsp_quant_nb(float *lsp, int order);
void lsp_unquant_nb(float *lsp, int order, unsigned int id);

#endif
