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


void lsp_quant_nb(float *lsp, float *qlsp, int order, FrameBits *bits)
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
   id = lsp_quant(qlsp, cdbk_nb, NB_CDBK_SIZE, order);
   printf ("qid = %d\n", id);
   frame_bits_pack(bits, id, 6);

   id = lsp_weight_quant(qlsp, quant_weight, cdbk_nb_low1, NB_CDBK_SIZE_LOW1, 5);
   printf ("qid = %d\n", id);
   frame_bits_pack(bits, id, 6);

   id = lsp_weight_quant(qlsp, quant_weight, cdbk_nb_low2, NB_CDBK_SIZE_LOW2, 5);
   printf ("qid = %d\n", id);
   frame_bits_pack(bits, id, 6);

   id = lsp_weight_quant(qlsp+5, quant_weight+5, cdbk_nb_high1, NB_CDBK_SIZE_HIGH1, 5);
   printf ("qid = %d\n", id);
   frame_bits_pack(bits, id, 6);

   id = lsp_weight_quant(qlsp+5, quant_weight+5, cdbk_nb_high2, NB_CDBK_SIZE_HIGH2, 5);
   printf ("qid = %d\n", id);
   frame_bits_pack(bits, id, 6);

   for (i=0;i<order;i++)
      qlsp[i]=lsp[i]-qlsp[i];
}

void lsp_unquant_nb(float *lsp, int order, FrameBits *bits)
{
   int i, id;
   for (i=0;i<order;i++)
      lsp[i]=0;


   id=frame_bits_unpack_unsigned(bits, 6);
   printf ("did = %d\n", id);
   for (i=0;i<10;i++)
      lsp[i] += cdbk_nb[id*10+i];

   id=frame_bits_unpack_unsigned(bits, 6);
   printf ("did = %d\n", id);
   for (i=0;i<5;i++)
      lsp[i] += cdbk_nb_low1[id*5+i];

   id=frame_bits_unpack_unsigned(bits, 6);
   printf ("did = %d\n", id);
   for (i=0;i<5;i++)
      lsp[i] += cdbk_nb_low2[id*5+i];

   id=frame_bits_unpack_unsigned(bits, 6);
   printf ("did = %d\n", id);
   for (i=0;i<5;i++)
      lsp[i+5] += cdbk_nb_high1[id*5+i];
   
   id=frame_bits_unpack_unsigned(bits, 6);
   printf ("did = %d\n", id);
   for (i=0;i<5;i++)
      lsp[i+5] += cdbk_nb_high2[id*5+i];
}
