/* Copyright (C) 2002 Jean-Marc Valin 
   File: cb_search.c

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



#include <stdlib.h>
#include "cb_search.h"
#include "filters.h"
#include <math.h>
#ifdef DEBUG
#include <stdio.h>
#endif
#include "stack_alloc.h"
#include "vq.h"

#ifndef min
# define min(a,b) ((a) < (b) ? (a) : (b))
#endif
#ifndef max
# define max(a,b) ((a) > (b) ? (a) : (b))
#endif




void split_cb_search_shape_sign(
float target[],			/* target vector */
float ak[],			/* LPCs for this subframe */
float awk1[],			/* Weighted LPCs for this subframe */
float awk2[],			/* Weighted LPCs for this subframe */
void *par,                      /* Codebook/search parameters*/
int   p,                        /* number of LPC coeffs */
int   nsf,                      /* number of samples in subframe */
float *exc,
float *r,
SpeexBits *bits,
float *stack,
int   complexity
)
{
   int i,j,k,m,n,q;
   float *resp;
   float *t, *e, *E;
   /*FIXME: Should make this dynamic*/
   float *tmp, *_ot[20], *_nt[20];
   float *ndist, *odist;
   int *itmp, *_nind[20], *_oind[20];
   float **ot, **nt;
   int **nind, **oind;
   int *ind;
   float *shape_cb;
   int shape_cb_size, subvect_size, nb_subvect;
   split_cb_params *params;
   int N=2;
   int *best_index;
   float *best_dist;
   int have_sign;

   ot=_ot;
   nt=_nt;
   oind=_oind;
   nind=_nind;
   N=complexity;
   if (N<1)
      N=1;
   if (N>10)
      N=10;

   params = (split_cb_params *) par;
   subvect_size = params->subvect_size;
   nb_subvect = params->nb_subvect;
   shape_cb_size = 1<<params->shape_bits;
   shape_cb = params->shape_cb;
   have_sign = params->have_sign;
   resp = PUSH(stack, shape_cb_size*subvect_size);
   t = PUSH(stack, nsf);
   e = PUSH(stack, nsf);
   E = PUSH(stack, shape_cb_size);
   /*FIXME: This breaks if sizeof(int) != sizeof(float) */
   ind = (int*)PUSH(stack, nb_subvect);

   tmp = PUSH(stack, 2*N*nsf);
   for (i=0;i<N;i++)
   {
      ot[i]=tmp;
      tmp += nsf;
      nt[i]=tmp;
      tmp += nsf;
   }

   /*FIXME: This breaks if sizeof(int) != sizeof(float) */
   best_index = (int*)PUSH(stack, N);
   best_dist = PUSH(stack, N);
   ndist = PUSH(stack, N);
   odist = PUSH(stack, N);
   
   /*FIXME: This breaks if sizeof(int) != sizeof(float) */
   itmp = (int*)PUSH(stack, 2*N*nb_subvect);
   for (i=0;i<N;i++)
   {
      nind[i]=itmp;
      itmp+=nb_subvect;
      oind[i]=itmp;
      itmp+=nb_subvect;
      for (j=0;j<nb_subvect;j++)
         nind[i][j]=oind[i][j]=-1;
   }

   for (j=0;j<N;j++)
      for (i=0;i<nsf;i++)
         ot[j][i]=target[i];

   for (i=0;i<nsf;i++)
      t[i]=target[i];

   /* Pre-compute codewords response and energy */
   for (i=0;i<shape_cb_size;i++)
   {
      float *res;
      float *shape;

      res = resp+i*subvect_size;
      shape = shape_cb+i*subvect_size;
      /* Compute codeword response */

      for(j=0;j<subvect_size;j++)
      {
         res[j]=0;
         for (k=0;k<=j;k++)
            res[j] += shape[k]*r[j-k];
      }
      E[i]=0;
      for(j=0;j<subvect_size;j++)
         E[i]+=res[j]*res[j];
   }

   for (j=0;j<N;j++)
      odist[j]=0;
   /*For all subvectors*/
   for (i=0;i<nb_subvect;i++)
   {
      /*"erase" nbest list*/
      for (j=0;j<N;j++)
         ndist[j]=-1;

      /*For all n-bests of previous subvector*/
      for (j=0;j<N;j++)
      {
         float *x=ot[j]+subvect_size*i;
         /*Find new n-best based on previous n-best j*/
         if (have_sign)
            vq_nbest_sign(x, resp, subvect_size, shape_cb_size, E, N, best_index, best_dist);
         else
            vq_nbest(x, resp, subvect_size, shape_cb_size, E, N, best_index, best_dist);

         /*For all new n-bests*/
         for (k=0;k<N;k++)
         {
            float *ct;
            float err=0;
            ct = ot[j];
            /*update target*/

            /*previous target*/
            for (m=i*subvect_size;m<(i+1)*subvect_size;m++)
               t[m]=ct[m];

            /* New code: update only enough of the target to calculate error*/
            {
               int rind;
               float *res;
               float sign=1;
               rind = best_index[k];
               if (rind>shape_cb_size)
               {
                  sign=-1;
                  rind-=shape_cb_size;
               }
               res = resp+rind*subvect_size;
               if (sign>0)
                  for (m=0;m<subvect_size;m++)
                     t[subvect_size*i+m] -= res[m];
               else
                  for (m=0;m<subvect_size;m++)
                     t[subvect_size*i+m] += res[m];
            }
            
            /*compute error (distance)*/
            err=odist[j];
            for (m=i*subvect_size;m<(i+1)*subvect_size;m++)
               err += t[m]*t[m];
            /*update n-best list*/
            if (err<ndist[N-1] || ndist[N-1]<-.5)
            {

               /*previous target (we don't care what happened before*/
               for (m=(i+1)*subvect_size;m<nsf;m++)
                  t[m]=ct[m];
               /* New code: update the rest of the target only if it's worth it */
               for (m=0;m<subvect_size;m++)
               {
                  float g;
                  int rind;
                  float sign=1;
                  rind = best_index[k];
                  if (rind>shape_cb_size)
                  {
                     sign=-1;
                     rind-=shape_cb_size;
                  }

                  g=sign*shape_cb[rind*subvect_size+m];
                  q=subvect_size-m;
                  for (n=subvect_size*(i+1);n<nsf;n++,q++)
                     t[n] -= g*r[q];
               }


               for (m=0;m<N;m++)
               {
                  if (err < ndist[m] || ndist[m]<-.5)
                  {
                     for (n=N-1;n>m;n--)
                     {
                        for (q=0;q<nsf;q++)
                           nt[n][q]=nt[n-1][q];
                        for (q=0;q<nb_subvect;q++)
                           nind[n][q]=nind[n-1][q];
                        ndist[n]=ndist[n-1];
                     }
                     for (q=0;q<nsf;q++)
                        nt[m][q]=t[q];
                     for (q=0;q<nb_subvect;q++)
                        nind[m][q]=oind[j][q];
                     nind[m][i]=best_index[k];
                     ndist[m]=err;
                     break;
                  }
               }
            }
         }
         if (i==0)
           break;
      }

      /*update old-new data*/
      /* just swap pointers instead of a long copy */
      {
         float **tmp;
         tmp=ot;
         ot=nt;
         nt=tmp;
      }
      for (j=0;j<N;j++)
         for (m=0;m<nb_subvect;m++)
            oind[j][m]=nind[j][m];
      for (j=0;j<N;j++)
         odist[j]=ndist[j];
   }

   /*save indices*/
   for (i=0;i<nb_subvect;i++)
   {
      ind[i]=nind[0][i];
      speex_bits_pack(bits,ind[i],params->shape_bits+have_sign);
   }
   
   /* Put everything back together */
   for (i=0;i<nb_subvect;i++)
   {
      int rind;
      float sign=1;
      rind = ind[i];
      if (rind>shape_cb_size)
      {
         sign=-1;
         rind-=shape_cb_size;
      }

      for (j=0;j<subvect_size;j++)
         e[subvect_size*i+j]=sign*shape_cb[rind*subvect_size+j];
   }   
   /* Update excitation */
   for (j=0;j<nsf;j++)
      exc[j]+=e[j];
   
   /* Update target */
   syn_percep_zero(e, ak, awk1, awk2, r, nsf,p, stack);
   for (j=0;j<nsf;j++)
      target[j]-=r[j];

}


void split_cb_shape_sign_unquant(
float *exc,
void *par,                      /* non-overlapping codebook */
int   nsf,                      /* number of samples in subframe */
SpeexBits *bits,
float *stack
)
{
   int i,j;
   int *ind, *signs;
   float *shape_cb;
   int shape_cb_size, subvect_size, nb_subvect;
   split_cb_params *params;
   int have_sign;

   params = (split_cb_params *) par;
   subvect_size = params->subvect_size;
   nb_subvect = params->nb_subvect;
   shape_cb_size = 1<<params->shape_bits;
   shape_cb = params->shape_cb;
   have_sign = params->have_sign;

   /*FIXME: This breaks if sizeof(int) != sizeof(float) */
   ind = (int*)PUSH(stack, nb_subvect);
   signs = (int*)PUSH(stack, nb_subvect);

   /* Decode codewords and gains */
   for (i=0;i<nb_subvect;i++)
   {
      if (have_sign)
         signs[i] = speex_bits_unpack_unsigned(bits, 1);
      else
         signs[i] = 0;
      ind[i] = speex_bits_unpack_unsigned(bits, params->shape_bits);
   }
   /* Compute decoded excitation */
   for (i=0;i<nb_subvect;i++)
   {
      float s=1;
      if (signs[i])
         s=-1;
      for (j=0;j<subvect_size;j++)
         exc[subvect_size*i+j]+=s*shape_cb[ind[i]*subvect_size+j];
   }
}

void noise_codebook_quant(
float target[],			/* target vector */
float ak[],			/* LPCs for this subframe */
float awk1[],			/* Weighted LPCs for this subframe */
float awk2[],			/* Weighted LPCs for this subframe */
void *par,                      /* Codebook/search parameters*/
int   p,                        /* number of LPC coeffs */
int   nsf,                      /* number of samples in subframe */
float *exc,
float *r,
SpeexBits *bits,
float *stack,
int   complexity
)
{
   int i;
   float *tmp=PUSH(stack, nsf);
   residue_percep_zero(target, ak, awk1, awk2, tmp, nsf, p, stack);

   for (i=0;i<nsf;i++)
      exc[i]+=tmp[i];
   for (i=0;i<nsf;i++)
      target[i]=0;

}


void noise_codebook_unquant(
float *exc,
void *par,                      /* non-overlapping codebook */
int   nsf,                      /* number of samples in subframe */
SpeexBits *bits,
float *stack
)
{
   int i;

   for (i=0;i<nsf;i++)
      exc[i]+=3*((((float)rand())/RAND_MAX)-.5);
}
