/* Copyright (C) 2002 Jean-Marc Valin 
   File: quant_lsp.c
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

#include "quant_lsp.h"
#include <math.h>
#include <stdio.h>

extern int lsp_nb_vqid[64];
static float quant_weight[MAX_LSP_SIZE];

/* Note: x is modified*/
static int lsp_quant(float *x, float *cdbk, int nbVec, int nbDim)
{
   int i,j;
   float dist, tmp;
   float best_dist=0;
   int best_id=0;
   float *ptr=cdbk;
   for (i=0;i<nbVec;i++)
   {
      dist=0;
      for (j=0;j<nbDim;j++)
      {
         tmp=(x[j]-*ptr++);
         dist+=tmp*tmp;
      }
      if (dist<best_dist || i==0)
      {
         best_dist=dist;
         best_id=i;
      }
   }

   for (j=0;j<nbDim;j++)
      x[j] -= cdbk[best_id*nbDim+j];
    
   return best_id;
}

/* Note: x is modified*/
static int lsp_weight_quant(float *x, float *weight, float *cdbk, int nbVec, int nbDim)
{
   int i,j;
   float dist, tmp;
   float best_dist=0;
   int best_id=0;
   float *ptr=cdbk;
   for (i=0;i<nbVec;i++)
   {
      dist=0;
      for (j=0;j<nbDim;j++)
      {
         tmp=(x[j]-*ptr++);
         dist+=weight[j]*tmp*tmp;
      }
      if (dist<best_dist || i==0)
      {
         best_dist=dist;
         best_id=i;
      }
   }
   
   for (j=0;j<nbDim;j++)
      x[j] -= cdbk[best_id*nbDim+j];
   return best_id;
}


void lsp_quant_nb(float *lsp, float *qlsp, int order, SpeexBits *bits)
{
   int i;
   float tmp1, tmp2;
   int id;

   for (i=0;i<order;i++)
      qlsp[i]=lsp[i];

   quant_weight[0] = 1/(qlsp[1]-qlsp[0]);
   quant_weight[order-1] = 1/(qlsp[order-1]-qlsp[order-2]);
   for (i=1;i<order-1;i++)
   {
#if 1
      tmp1 = 1/((.15+qlsp[i]-qlsp[i-1])*(.15+qlsp[i]-qlsp[i-1]));
      tmp2 = 1/((.15+qlsp[i+1]-qlsp[i])*(.15+qlsp[i+1]-qlsp[i]));
#else
      tmp1 = 1/(qlsp[i]-qlsp[i-1]);
      tmp2 = 1/(qlsp[i+1]-qlsp[i]);
#endif
      quant_weight[i] = tmp1 > tmp2 ? tmp1 : tmp2;
   }
   id = lsp_quant(qlsp, cdbk_nb, NB_CDBK_SIZE, order);
   speex_bits_pack(bits, id, 6);

   id = lsp_weight_quant(qlsp, quant_weight, cdbk_nb_low1, NB_CDBK_SIZE_LOW1, 5);
   speex_bits_pack(bits, id, 6);

   id = lsp_weight_quant(qlsp, quant_weight, cdbk_nb_low2, NB_CDBK_SIZE_LOW2, 5);
   speex_bits_pack(bits, id, 6);

   id = lsp_weight_quant(qlsp+5, quant_weight+5, cdbk_nb_high1, NB_CDBK_SIZE_HIGH1, 5);
   speex_bits_pack(bits, id, 6);

   id = lsp_weight_quant(qlsp+5, quant_weight+5, cdbk_nb_high2, NB_CDBK_SIZE_HIGH2, 5);
   speex_bits_pack(bits, id, 6);

   for (i=0;i<order;i++)
      qlsp[i]=lsp[i]-qlsp[i];
}

void lsp_unquant_nb(float *lsp, int order, SpeexBits *bits)
{
   int i, id;
   for (i=0;i<order;i++)
      lsp[i]=0;


   id=speex_bits_unpack_unsigned(bits, 6);
   for (i=0;i<10;i++)
      lsp[i] += cdbk_nb[id*10+i];

   id=speex_bits_unpack_unsigned(bits, 6);
   for (i=0;i<5;i++)
      lsp[i] += cdbk_nb_low1[id*5+i];

   id=speex_bits_unpack_unsigned(bits, 6);
   for (i=0;i<5;i++)
      lsp[i] += cdbk_nb_low2[id*5+i];

   id=speex_bits_unpack_unsigned(bits, 6);
   for (i=0;i<5;i++)
      lsp[i+5] += cdbk_nb_high1[id*5+i];
   
   id=speex_bits_unpack_unsigned(bits, 6);
   for (i=0;i<5;i++)
      lsp[i+5] += cdbk_nb_high2[id*5+i];
}


extern float lsp_cdbk_wb[];
extern float lsp_cdbk_wb11[];
extern float lsp_cdbk_wb12[];
extern float lsp_cdbk_wb21[];
extern float lsp_cdbk_wb22[];
extern float lsp_cdbk_wb31[];
extern float lsp_cdbk_wb32[];
extern float lsp_cdbk_wb41[];
extern float lsp_cdbk_wb42[];

void lsp_quant_wb(float *lsp, float *qlsp, int order, SpeexBits *bits)
{
   int i;
   float tmp1, tmp2;
   int id;
   for (i=0;i<order;i++)
      qlsp[i]=lsp[i];

   quant_weight[0] = 1/(qlsp[1]-qlsp[0]);
   quant_weight[order-1] = 1/(qlsp[order-1]-qlsp[order-2]);
   for (i=1;i<order-1;i++)
   {
      tmp1 = 1/(qlsp[i]-qlsp[i-1]);
      tmp2 = 1/(qlsp[i+1]-qlsp[i]);
      quant_weight[i] = tmp1 > tmp2 ? tmp1 : tmp2;
   }
   id = lsp_quant(qlsp, lsp_cdbk_wb, 64, order);
   speex_bits_pack(bits, id, 6);

   id = lsp_weight_quant(qlsp, quant_weight, lsp_cdbk_wb11, 64, 4);
   speex_bits_pack(bits, id, 6);

   id = lsp_weight_quant(qlsp, quant_weight, lsp_cdbk_wb12, 64, 4);
   speex_bits_pack(bits, id, 6);

   id = lsp_weight_quant(qlsp+4, quant_weight, lsp_cdbk_wb21, 64, 4);
   speex_bits_pack(bits, id, 6);

   id = lsp_weight_quant(qlsp+4, quant_weight, lsp_cdbk_wb22, 64, 4);
   speex_bits_pack(bits, id, 6);

   id = lsp_weight_quant(qlsp+8, quant_weight, lsp_cdbk_wb31, 64, 4);
   speex_bits_pack(bits, id, 6);

   id = lsp_weight_quant(qlsp+8, quant_weight, lsp_cdbk_wb32, 16, 4);
   speex_bits_pack(bits, id, 4);

   id = lsp_weight_quant(qlsp+12, quant_weight, lsp_cdbk_wb41, 64, 4);
   speex_bits_pack(bits, id, 6);

   id = lsp_weight_quant(qlsp+12, quant_weight, lsp_cdbk_wb42, 16, 4);
   speex_bits_pack(bits, id, 4);

   for (i=0;i<order;i++)
      qlsp[i]=lsp[i]-qlsp[i];

}


void lsp_unquant_wb(float *lsp, int order, SpeexBits *bits)
{

   int i, id;
   for (i=0;i<order;i++)
      lsp[i]=0;


   id=speex_bits_unpack_unsigned(bits, 6);
   for (i=0;i<16;i++)
      lsp[i] += lsp_cdbk_wb[id*16+i];


   id=speex_bits_unpack_unsigned(bits, 6);
   for (i=0;i<4;i++)
      lsp[i] += lsp_cdbk_wb11[id*4+i];

   id=speex_bits_unpack_unsigned(bits, 6);
   for (i=0;i<4;i++)
      lsp[i] += lsp_cdbk_wb12[id*4+i];


   id=speex_bits_unpack_unsigned(bits, 6);
   for (i=0;i<4;i++)
      lsp[i+4] += lsp_cdbk_wb21[id*4+i];

   id=speex_bits_unpack_unsigned(bits, 6);
   for (i=0;i<4;i++)
      lsp[i+4] += lsp_cdbk_wb22[id*4+i];


   id=speex_bits_unpack_unsigned(bits, 6);
   for (i=0;i<4;i++)
      lsp[i+8] += lsp_cdbk_wb31[id*4+i];

   id=speex_bits_unpack_unsigned(bits, 4);
   for (i=0;i<4;i++)
      lsp[i+8] += lsp_cdbk_wb32[id*4+i];


   id=speex_bits_unpack_unsigned(bits, 6);
   for (i=0;i<4;i++)
      lsp[i+12] += lsp_cdbk_wb41[id*4+i];

   id=speex_bits_unpack_unsigned(bits, 4);
   for (i=0;i<4;i++)
      lsp[i+12] += lsp_cdbk_wb42[id*4+i];

}

extern float high_lsp_cdbk[];
extern float high_lsp_cdbk2[];


void lsp_quant_high(float *lsp, float *qlsp, int order, SpeexBits *bits)
{
   int i;
   float tmp1, tmp2;
   int id;
   for (i=0;i<order;i++)
      qlsp[i]=lsp[i];

   quant_weight[0] = 1/(qlsp[1]-qlsp[0]);
   quant_weight[order-1] = 1/(qlsp[order-1]-qlsp[order-2]);
   for (i=1;i<order-1;i++)
   {
      tmp1 = 1/(qlsp[i]-qlsp[i-1]);
      tmp2 = 1/(qlsp[i+1]-qlsp[i]);
      quant_weight[i] = tmp1 > tmp2 ? tmp1 : tmp2;
   }
   id = lsp_quant(qlsp, high_lsp_cdbk, 64, order);
   speex_bits_pack(bits, id, 6);

   id = lsp_weight_quant(qlsp, quant_weight, high_lsp_cdbk2, 64, order);
   speex_bits_pack(bits, id, 6);

   for (i=0;i<order;i++)
      qlsp[i]=lsp[i]-qlsp[i];
}

void lsp_unquant_high(float *lsp, int order, SpeexBits *bits)
{

   int i, id;
   for (i=0;i<order;i++)
      lsp[i]=0;


   id=speex_bits_unpack_unsigned(bits, 6);
   for (i=0;i<order;i++)
      lsp[i] += high_lsp_cdbk[id*order+i];


   id=speex_bits_unpack_unsigned(bits, 6);
   for (i=0;i<order;i++)
      lsp[i] += high_lsp_cdbk2[id*order+i];
}
