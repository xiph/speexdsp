/* Copyright (C) 2002 Jean-Marc Valin 
   File: quant_lsp.c
   LSP vector quantization
*/

#include "quant_lsp.h"
#include <math.h>
#include <stdio.h>

static float quant_weight[MAX_LSP_SIZE];

int lsp_quant(float *x, float *cdbk, int nbVec, int nbDim)
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

int lsp_weight_quant(float *x, float *weight, float *cdbk, int nbVec, int nbDim)
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


unsigned int lsp_quant_nb(float *lsp, int order)
{
   int i;
   float tmp1, tmp2;
   unsigned int id;
   quant_weight[0] = 1/(lsp[1]-lsp[0]);
   quant_weight[order-1] = 1/(lsp[order-1]-lsp[order-2]);
   for (i=1;i<order-1;i++)
   {
      tmp1 = 1/(lsp[i]-lsp[i-1]);
      tmp2 = 1/(lsp[i+1]-lsp[i]);
      quant_weight[i] = tmp1 > tmp2 ? tmp1 : tmp2;
   }
   id = lsp_quant(lsp, cdbk_nb, NB_CDBK_SIZE, order);
   
   id *= NB_CDBK_SIZE_LOW1;
   id += lsp_weight_quant(lsp, quant_weight, cdbk_nb_low1, NB_CDBK_SIZE_LOW1, 5);
   
   id *= NB_CDBK_SIZE_LOW2;
   id += lsp_weight_quant(lsp, quant_weight, cdbk_nb_low2, NB_CDBK_SIZE_LOW2, 5);
   
   id *= NB_CDBK_SIZE_HIGH1;
   id += lsp_weight_quant(lsp+5, quant_weight+5, cdbk_nb_high1, NB_CDBK_SIZE_HIGH1, 5);
   
   id *= NB_CDBK_SIZE_HIGH2;
   id += lsp_weight_quant(lsp+5, quant_weight+5, cdbk_nb_high2, NB_CDBK_SIZE_HIGH2, 5);
   
   return id;
}

void lsp_unquant_nb(float *lsp, int order, unsigned int id)
{
   int i, tmp;
   for (i=0;i<order;i++)
      lsp[i]=0;


   tmp = id % NB_CDBK_SIZE_HIGH2;
   for (i=0;i<5;i++)
      lsp[i+5] += cdbk_nb_high2[tmp*5+i];
   id /= NB_CDBK_SIZE_HIGH2;
   
   tmp = id % NB_CDBK_SIZE_HIGH1;
   for (i=0;i<5;i++)
      lsp[i+5] += cdbk_nb_high1[tmp*5+i];
   id /= NB_CDBK_SIZE_HIGH1;
   
   tmp = id % NB_CDBK_SIZE_LOW2;
   for (i=0;i<5;i++)
      lsp[i] += cdbk_nb_low2[tmp*5+i];
   id /= NB_CDBK_SIZE_LOW2;

   tmp = id % NB_CDBK_SIZE_LOW1;
   for (i=0;i<5;i++)
      lsp[i] += cdbk_nb_low1[tmp*5+i];
   id /= NB_CDBK_SIZE_LOW1;

   tmp=id;
   for (i=0;i<10;i++)
      lsp[i] += cdbk_nb[tmp*10+i];
   
}
