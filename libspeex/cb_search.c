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


#include "cb_search.h"
#include "filters.h"
#include "stack_alloc.h"
#include "vq.h"
#include "misc.h"
#include <stdio.h>

void split_cb_search_shape_sign(
spx_sig_t target[],			/* target vector */
spx_coef_t ak[],			/* LPCs for this subframe */
spx_coef_t awk1[],			/* Weighted LPCs for this subframe */
spx_coef_t awk2[],			/* Weighted LPCs for this subframe */
void *par,                      /* Codebook/search parameters*/
int   p,                        /* number of LPC coeffs */
int   nsf,                      /* number of samples in subframe */
spx_sig_t *exc,
spx_sig_t *r,
SpeexBits *bits,
char *stack,
int   complexity
)
{
   int i,j,k,m,n,q;
   spx_word16_t *resp;
   spx_word16_t *t;
   spx_sig_t *e, *r2;
   spx_word32_t *E;
   spx_word16_t *tmp;
   spx_word32_t *ndist, *odist;
   int *itmp;
   spx_word16_t **ot, **nt;
   int **nind, **oind;
   int *ind;
   signed char *shape_cb;
   int shape_cb_size, subvect_size, nb_subvect;
   split_cb_params *params;
   int N=2;
   int *best_index;
   spx_word32_t *best_dist;
   int have_sign;

   N=complexity;
   if (N>10)
      N=10;

   ot=PUSH(stack, N, spx_word16_t*);
   nt=PUSH(stack, N, spx_word16_t*);
   oind=PUSH(stack, N, int*);
   nind=PUSH(stack, N, int*);

   params = (split_cb_params *) par;
   subvect_size = params->subvect_size;
   nb_subvect = params->nb_subvect;
   shape_cb_size = 1<<params->shape_bits;
   shape_cb = params->shape_cb;
   have_sign = params->have_sign;
   resp = PUSH(stack, shape_cb_size*subvect_size, spx_word16_t);
   t = PUSH(stack, nsf, spx_word16_t);
   e = PUSH(stack, nsf, spx_sig_t);
   r2 = PUSH(stack, nsf, spx_sig_t);
   E = PUSH(stack, shape_cb_size, spx_word32_t);
   ind = PUSH(stack, nb_subvect, int);

   tmp = PUSH(stack, 2*N*nsf, spx_word16_t);
   for (i=0;i<N;i++)
   {
      ot[i]=tmp;
      tmp += nsf;
      nt[i]=tmp;
      tmp += nsf;
   }

   best_index = PUSH(stack, N, int);
   best_dist = PUSH(stack, N, spx_word32_t);
   ndist = PUSH(stack, N, spx_word32_t);
   odist = PUSH(stack, N, spx_word32_t);
   
   itmp = PUSH(stack, 2*N*nb_subvect, int);
   for (i=0;i<N;i++)
   {
      nind[i]=itmp;
      itmp+=nb_subvect;
      oind[i]=itmp;
      itmp+=nb_subvect;
      for (j=0;j<nb_subvect;j++)
         nind[i][j]=oind[i][j]=-1;
   }

   for (i=0;i<nsf;i++)
      t[i]=SHR(target[i],6);

   for (j=0;j<N;j++)
      for (i=0;i<nsf;i++)
         ot[j][i]=t[i];

   /*for (i=0;i<nsf;i++)
     printf ("%d\n", (int)t[i]);*/

   /* Pre-compute codewords response and energy */
   for (i=0;i<shape_cb_size;i++)
   {
      spx_word16_t *res;
      signed char *shape;

      res = resp+i*subvect_size;
      shape = shape_cb+i*subvect_size;

      /* Compute codeword response using convolution with impulse response */
      for(j=0;j<subvect_size;j++)
      {
         spx_word32_t resj=0;
         for (k=0;k<=j;k++)
            resj += shape[k]*r[j-k];
         resj *= 0.03125;
         
         res[j] = SHR(resj,6);
         /*printf ("%d\n", (int)res[j]);*/
      }
      
      /* Compute codeword energy */
      E[i]=0;
      for(j=0;j<subvect_size;j++)
         E[i]+=MULT16_16(res[j],res[j]);
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
         spx_word16_t *x=ot[j]+subvect_size*i;
         /*Find new n-best based on previous n-best j*/
         if (have_sign)
            vq_nbest_sign(x, resp, subvect_size, shape_cb_size, E, N, best_index, best_dist);
         else
            vq_nbest(x, resp, subvect_size, shape_cb_size, E, N, best_index, best_dist);

         /*For all new n-bests*/
         for (k=0;k<N;k++)
         {
            spx_word16_t *ct;
            spx_word32_t err=0;
            ct = ot[j];
            /*update target*/

            /*previous target*/
            for (m=i*subvect_size;m<(i+1)*subvect_size;m++)
               t[m]=ct[m];

            /* New code: update only enough of the target to calculate error*/
            {
               int rind;
               spx_word16_t *res;
               spx_word16_t sign=1;
               rind = best_index[k];
               if (rind>=shape_cb_size)
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
                  spx_word16_t g;
                  int rind;
                  spx_word16_t sign=1;
                  rind = best_index[k];
                  if (rind>=shape_cb_size)
                  {
                     sign=-1;
                     rind-=shape_cb_size;
                  }

                  q=subvect_size-m;
#ifdef FIXED_POINT
                  g=sign*shape_cb[rind*subvect_size+m];
                  for (n=subvect_size*(i+1);n<nsf;n++,q++)
                     t[n] -= SHR(MULT16_16(g,r[q]),11);
#else
                  g=sign*0.03125*shape_cb[rind*subvect_size+m];
                  for (n=subvect_size*(i+1);n<nsf;n++,q++)
                     t[n] -= g*r[q];
#endif
               }


               for (m=0;m<N;m++)
               {
                  if (err < ndist[m] || ndist[m]<-.5)
                  {
                     for (n=N-1;n>m;n--)
                     {
                        for (q=(i+1)*subvect_size;q<nsf;q++)
                           nt[n][q]=nt[n-1][q];
                        for (q=0;q<nb_subvect;q++)
                           nind[n][q]=nind[n-1][q];
                        ndist[n]=ndist[n-1];
                     }
                     for (q=(i+1)*subvect_size;q<nsf;q++)
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
         spx_word16_t **tmp2;
         tmp2=ot;
         ot=nt;
         nt=tmp2;
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
      spx_word16_t sign=1;
      rind = ind[i];
      if (rind>=shape_cb_size)
      {
         sign=-1;
         rind-=shape_cb_size;
      }
#ifdef FIXED_POINT
      if (sign==1)
      {
         for (j=0;j<subvect_size;j++)
            e[subvect_size*i+j]=SHL((spx_word32_t)shape_cb[rind*subvect_size+j],SIG_SHIFT-5);
      } else {
         for (j=0;j<subvect_size;j++)
            e[subvect_size*i+j]=-SHL((spx_word32_t)shape_cb[rind*subvect_size+j],SIG_SHIFT-5);
      }
#else
      for (j=0;j<subvect_size;j++)
         e[subvect_size*i+j]=sign*0.03125*SIG_SCALING*shape_cb[rind*subvect_size+j];
#endif
   }   
   /* Update excitation */
   for (j=0;j<nsf;j++)
      exc[j]+=e[j];
   
   /* Update target */
   syn_percep_zero(e, ak, awk1, awk2, r2, nsf,p, stack);
   for (j=0;j<nsf;j++)
      target[j]-=r2[j];

}


void split_cb_shape_sign_unquant(
spx_sig_t *exc,
void *par,                      /* non-overlapping codebook */
int   nsf,                      /* number of samples in subframe */
SpeexBits *bits,
char *stack
)
{
   int i,j;
   int *ind, *signs;
   signed char *shape_cb;
   int shape_cb_size, subvect_size, nb_subvect;
   split_cb_params *params;
   int have_sign;

   params = (split_cb_params *) par;
   subvect_size = params->subvect_size;
   nb_subvect = params->nb_subvect;
   shape_cb_size = 1<<params->shape_bits;
   shape_cb = params->shape_cb;
   have_sign = params->have_sign;

   ind = PUSH(stack, nb_subvect, int);
   signs = PUSH(stack, nb_subvect, int);

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
      spx_word16_t s=1;
      if (signs[i])
         s=-1;
#ifdef FIXED_POINT
      if (s==1)
      {
         for (j=0;j<subvect_size;j++)
            exc[subvect_size*i+j]=SHL((spx_word32_t)shape_cb[ind[i]*subvect_size+j],SIG_SHIFT-5);
      } else {
         for (j=0;j<subvect_size;j++)
            exc[subvect_size*i+j]=-SHL((spx_word32_t)shape_cb[ind[i]*subvect_size+j],SIG_SHIFT-5);
      }
#else
      for (j=0;j<subvect_size;j++)
         exc[subvect_size*i+j]+=s*0.03125*shape_cb[ind[i]*subvect_size+j];      
#endif
   }
}

void noise_codebook_quant(
spx_sig_t target[],			/* target vector */
spx_coef_t ak[],			/* LPCs for this subframe */
spx_coef_t awk1[],			/* Weighted LPCs for this subframe */
spx_coef_t awk2[],			/* Weighted LPCs for this subframe */
void *par,                      /* Codebook/search parameters*/
int   p,                        /* number of LPC coeffs */
int   nsf,                      /* number of samples in subframe */
spx_sig_t *exc,
spx_sig_t *r,
SpeexBits *bits,
char *stack,
int   complexity
)
{
   int i;
   spx_sig_t *tmp=PUSH(stack, nsf, spx_sig_t);
   residue_percep_zero(target, ak, awk1, awk2, tmp, nsf, p, stack);

   for (i=0;i<nsf;i++)
      exc[i]+=tmp[i];
   for (i=0;i<nsf;i++)
      target[i]=0;

}


void noise_codebook_unquant(
spx_sig_t *exc,
void *par,                      /* non-overlapping codebook */
int   nsf,                      /* number of samples in subframe */
SpeexBits *bits,
char *stack
)
{
   speex_rand_vec(1, exc, nsf);
}
