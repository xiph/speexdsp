/* Copyright (C) 2002 Jean-Marc Valin 
   File: quant_lsp.c
   LSP vector quantization

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:
   
   - Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
   
   - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
   
   - Neither the name of the Xiph.org Foundation nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.
   
   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE FOUNDATION OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "quant_lsp.h"
#include <math.h>
#include "misc.h"

extern int lsp_nb_vqid[64];

/* FIXME: Get rid of this kludge quick before someone gets hurt */
static spx_word16_t quant_weight[MAX_LSP_SIZE];

#ifdef FIXED_POINT
#define LSP_SCALE 8192
#define LSP_OVERSCALE 32
#else
#define LSP_SCALE 256
#define LSP_OVERSCALE 1
#endif

/* Note: x is modified*/
static int lsp_quant(spx_word16_t *x, signed char *cdbk, int nbVec, int nbDim)
{
   int i,j;
   spx_word32_t dist;
   spx_word16_t tmp;
   spx_word32_t best_dist=0;
   int best_id=0;
   signed char *ptr=cdbk;
   for (i=0;i<nbVec;i++)
   {
      dist=0;
      for (j=0;j<nbDim;j++)
      {
         tmp=(x[j]-SHL((spx_word16_t)*ptr++,5));
         dist+=MULT16_16(tmp,tmp);
      }
      if (dist<best_dist || i==0)
      {
         best_dist=dist;
         best_id=i;
      }
   }

   for (j=0;j<nbDim;j++)
      x[j] -= SHL((spx_word16_t)cdbk[best_id*nbDim+j],5);
    
   return best_id;
}

/* Note: x is modified*/
static int lsp_weight_quant(spx_word16_t *x, spx_word16_t *weight, signed char *cdbk, int nbVec, int nbDim)
{
   int i,j;
   float dist;
   spx_word16_t tmp;
   float best_dist=0;
   int best_id=0;
   signed char *ptr=cdbk;
   for (i=0;i<nbVec;i++)
   {
      dist=0;
      for (j=0;j<nbDim;j++)
      {
         tmp=(x[j]-SHL((spx_word16_t)*ptr++,5));
         dist+=MULT16_32_Q15(weight[j],MULT16_16(tmp,tmp));
      }
      if (dist<best_dist || i==0)
      {
         best_dist=dist;
         best_id=i;
      }
   }
   
   for (j=0;j<nbDim;j++)
      x[j] -= SHL((spx_word16_t)cdbk[best_id*nbDim+j],5);
   return best_id;
}


void lsp_quant_nb(float *lsp, float *qlsp, int order, SpeexBits *bits)
{
   int i;
   float tmp1, tmp2;
   int id;
   /* FIXME: get rid of that static allocation */
   spx_word16_t nlsp[10];
   
   for (i=0;i<order;i++)
      qlsp[i]=lsp[i];

   quant_weight[0] = 10/(qlsp[1]-qlsp[0]);
   quant_weight[order-1] = 10/(qlsp[order-1]-qlsp[order-2]);
   for (i=1;i<order-1;i++)
   {
#if 1
      tmp1 = 10/((.15+qlsp[i]-qlsp[i-1])*(.15+qlsp[i]-qlsp[i-1]));
      tmp2 = 10/((.15+qlsp[i+1]-qlsp[i])*(.15+qlsp[i+1]-qlsp[i]));
#else
      tmp1 = 10/(qlsp[i]-qlsp[i-1]);
      tmp2 = 10/(qlsp[i+1]-qlsp[i]);
#endif
      quant_weight[i] = tmp1 > tmp2 ? tmp1 : tmp2;
   }
   for (i=0;i<order;i++)
      qlsp[i]-=(.25*i+.25);

   for (i=0;i<order;i++)
      nlsp[i] = LSP_SCALE*qlsp[i];

   id = lsp_quant(nlsp, cdbk_nb, NB_CDBK_SIZE, order);
   speex_bits_pack(bits, id, 6);

   for (i=0;i<order;i++)
      nlsp[i]*=2;
 
   id = lsp_weight_quant(nlsp, quant_weight, cdbk_nb_low1, NB_CDBK_SIZE_LOW1, 5);
   speex_bits_pack(bits, id, 6);

   for (i=0;i<5;i++)
      nlsp[i]*=2;

   id = lsp_weight_quant(nlsp, quant_weight, cdbk_nb_low2, NB_CDBK_SIZE_LOW2, 5);
   speex_bits_pack(bits, id, 6);

   id = lsp_weight_quant(nlsp+5, quant_weight+5, cdbk_nb_high1, NB_CDBK_SIZE_HIGH1, 5);
   speex_bits_pack(bits, id, 6);

   for (i=5;i<10;i++)
      nlsp[i]*=2;

   id = lsp_weight_quant(nlsp+5, quant_weight+5, cdbk_nb_high2, NB_CDBK_SIZE_HIGH2, 5);
   speex_bits_pack(bits, id, 6);

   for (i=0;i<order;i++)
      qlsp[i]=nlsp[i] * (.00097656/LSP_OVERSCALE);

   for (i=0;i<order;i++)
      qlsp[i]=lsp[i]-qlsp[i];
}

void lsp_unquant_nb(float *lsp, int order, SpeexBits *bits)
{
   int i, id;
   for (i=0;i<order;i++)
      lsp[i]=.25*i+.25;


   id=speex_bits_unpack_unsigned(bits, 6);
   for (i=0;i<10;i++)
      lsp[i] += 0.0039062*cdbk_nb[id*10+i];

   id=speex_bits_unpack_unsigned(bits, 6);
   for (i=0;i<5;i++)
      lsp[i] += 0.0019531 * cdbk_nb_low1[id*5+i];

   id=speex_bits_unpack_unsigned(bits, 6);
   for (i=0;i<5;i++)
      lsp[i] +=  0.00097656 * cdbk_nb_low2[id*5+i];

   id=speex_bits_unpack_unsigned(bits, 6);
   for (i=0;i<5;i++)
      lsp[i+5] += 0.0019531 * cdbk_nb_high1[id*5+i];
   
   id=speex_bits_unpack_unsigned(bits, 6);
   for (i=0;i<5;i++)
      lsp[i+5] += 0.00097656 * cdbk_nb_high2[id*5+i];
}


void lsp_quant_lbr(float *lsp, float *qlsp, int order, SpeexBits *bits)
{
   int i;
   float tmp1, tmp2;
   int id;
   spx_word16_t nlsp[10];

   for (i=0;i<order;i++)
      qlsp[i]=lsp[i];

   quant_weight[0] = 10/(qlsp[1]-qlsp[0]);
   quant_weight[order-1] = 10/(qlsp[order-1]-qlsp[order-2]);
   for (i=1;i<order-1;i++)
   {
#if 1
      tmp1 = 10/((.15+qlsp[i]-qlsp[i-1])*(.15+qlsp[i]-qlsp[i-1]));
      tmp2 = 10/((.15+qlsp[i+1]-qlsp[i])*(.15+qlsp[i+1]-qlsp[i]));
#else
      tmp1 = 10/(qlsp[i]-qlsp[i-1]);
      tmp2 = 10/(qlsp[i+1]-qlsp[i]);
#endif
      quant_weight[i] = tmp1 > tmp2 ? tmp1 : tmp2;
   }

   for (i=0;i<order;i++)
      qlsp[i]-=(.25*i+.25);
   for (i=0;i<order;i++)
      nlsp[i]=qlsp[i]*LSP_SCALE;
   
   id = lsp_quant(nlsp, cdbk_nb, NB_CDBK_SIZE, order);
   speex_bits_pack(bits, id, 6);

   for (i=0;i<order;i++)
      nlsp[i]*=2;
   
   id = lsp_weight_quant(nlsp, quant_weight, cdbk_nb_low1, NB_CDBK_SIZE_LOW1, 5);
   speex_bits_pack(bits, id, 6);

   id = lsp_weight_quant(nlsp+5, quant_weight+5, cdbk_nb_high1, NB_CDBK_SIZE_HIGH1, 5);
   speex_bits_pack(bits, id, 6);

   for (i=0;i<order;i++)
      qlsp[i] = nlsp[i]*(0.0019531/LSP_OVERSCALE);

   for (i=0;i<order;i++)
      qlsp[i]=lsp[i]-qlsp[i];
}

void lsp_unquant_lbr(float *lsp, int order, SpeexBits *bits)
{
   int i, id;
   for (i=0;i<order;i++)
      lsp[i]=.25*i+.25;


   id=speex_bits_unpack_unsigned(bits, 6);
   for (i=0;i<10;i++)
      lsp[i] += 0.0039062*cdbk_nb[id*10+i];

   id=speex_bits_unpack_unsigned(bits, 6);
   for (i=0;i<5;i++)
      lsp[i] += 0.0019531*cdbk_nb_low1[id*5+i];

   id=speex_bits_unpack_unsigned(bits, 6);
   for (i=0;i<5;i++)
      lsp[i+5] += 0.0019531*cdbk_nb_high1[id*5+i];
   
}


extern signed char high_lsp_cdbk[];
extern signed char high_lsp_cdbk2[];


void lsp_quant_high(float *lsp, float *qlsp, int order, SpeexBits *bits)
{
   int i;
   float tmp1, tmp2;
   int id;
   spx_word16_t nlsp[10];

   for (i=0;i<order;i++)
      qlsp[i]=lsp[i];

   quant_weight[0] = 10/(qlsp[1]-qlsp[0]);
   quant_weight[order-1] = 10/(qlsp[order-1]-qlsp[order-2]);
   for (i=1;i<order-1;i++)
   {
      tmp1 = 10/(qlsp[i]-qlsp[i-1]);
      tmp2 = 10/(qlsp[i+1]-qlsp[i]);
      quant_weight[i] = tmp1 > tmp2 ? tmp1 : tmp2;
   }

   for (i=0;i<order;i++)
      qlsp[i]-=(.3125*i+.75);
   for (i=0;i<order;i++)
      nlsp[i] = qlsp[i]*LSP_SCALE;

   id = lsp_quant(nlsp, high_lsp_cdbk, 64, order);
   speex_bits_pack(bits, id, 6);

   for (i=0;i<order;i++)
      nlsp[i]*=2;

   id = lsp_weight_quant(nlsp, quant_weight, high_lsp_cdbk2, 64, order);
   speex_bits_pack(bits, id, 6);

   for (i=0;i<order;i++)
      qlsp[i] = nlsp[i]*(0.0019531/LSP_OVERSCALE);

   for (i=0;i<order;i++)
      qlsp[i]=lsp[i]-qlsp[i];
}

void lsp_unquant_high(float *lsp, int order, SpeexBits *bits)
{

   int i, id;
   for (i=0;i<order;i++)
      lsp[i]=.3125*i+.75;


   id=speex_bits_unpack_unsigned(bits, 6);
   for (i=0;i<order;i++)
      lsp[i] += 0.0039062*high_lsp_cdbk[id*order+i];


   id=speex_bits_unpack_unsigned(bits, 6);
   for (i=0;i<order;i++)
      lsp[i] += 0.0019531*high_lsp_cdbk2[id*order+i];
}


#ifdef EPIC_48K

extern signed char cdbk_lsp_vlbr[5120];
extern signed char cdbk_lsp2_vlbr[160];

void lsp_quant_48k(float *lsp, float *qlsp, int order, SpeexBits *bits)
{
   int i;
   float tmp1, tmp2;
   int id;
   spx_word16_t nlsp[10];

   for (i=0;i<order;i++)
      qlsp[i]=lsp[i];

   quant_weight[0] = 10/(qlsp[1]-qlsp[0]);
   quant_weight[order-1] = 10/(qlsp[order-1]-qlsp[order-2]);
   for (i=1;i<order-1;i++)
   {
#if 1
      tmp1 = 10/((.15+qlsp[i]-qlsp[i-1])*(.15+qlsp[i]-qlsp[i-1]));
      tmp2 = 10/((.15+qlsp[i+1]-qlsp[i])*(.15+qlsp[i+1]-qlsp[i]));
#else
      tmp1 = 10/(qlsp[i]-qlsp[i-1]);
      tmp2 = 10/(qlsp[i+1]-qlsp[i]);
#endif
      quant_weight[i] = tmp1 > tmp2 ? tmp1 : tmp2;
   }

   for (i=0;i<order;i++)
      qlsp[i]-=(.25*i+.3125);
   for (i=0;i<order;i++)
      nlsp[i]=qlsp[i]*LSP_SCALE;
   
   id = lsp_quant(qlsp, cdbk_lsp_vlbr, 512, order);
   speex_bits_pack(bits, id, 9);

   for (i=0;i<order;i++)
      qlsp[i]*=4;
   
   id = lsp_weight_quant(qlsp, quant_weight, cdbk_lsp2_vlbr, 16, 10);
   speex_bits_pack(bits, id, 4);

   for (i=0;i<order;i++)
      qlsp[i]=nlsp[i]*(0.00097655/LSP_OVERSCALE);

   for (i=0;i<order;i++)
      qlsp[i]=lsp[i]-qlsp[i];
}

void lsp_unquant_48k(float *lsp, int order, SpeexBits *bits)
{
   int i, id;
   for (i=0;i<order;i++)
      lsp[i]=.25*i+.3125;


   id=speex_bits_unpack_unsigned(bits, 9);
   for (i=0;i<10;i++)
      lsp[i] += 0.0039062*cdbk_lsp_vlbr[id*10+i];

   id=speex_bits_unpack_unsigned(bits, 4);
   for (i=0;i<10;i++)
      lsp[i] += 0.00097655*cdbk_lsp2_vlbr[id*10+i];
   
}

#endif
