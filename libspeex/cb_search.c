/* Copyright (C) 2002 Jean-Marc Valin 
   File: cb_search.c

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



#include <stdlib.h>
#include <cb_search.h>
#include "filters.h"
#include <math.h>
#include <stdio.h>
#include "stack_alloc.h"
#include "vq.h"

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))



void split_cb_search_nogain(
float target[],			/* target vector */
float ak[],			/* LPCs for this subframe */
float awk1[],			/* Weighted LPCs for this subframe */
float awk2[],			/* Weighted LPCs for this subframe */
void *par,                      /* Codebook/search parameters*/
int   p,                        /* number of LPC coeffs */
int   nsf,                      /* number of samples in subframe */
float *exc,
SpeexBits *bits,
float *stack,
int   complexity
)
{
   int i,j,k,m,n,q;
   float *resp;
   float *t, *r, *e, *E;
   /*FIXME: Should make this dynamic*/
   float *tmp, *ot[20], *nt[20];
   float *ndist, *odist;
   int *itmp, *nind[20], *oind[20];

   int *ind;
   float *shape_cb;
   int shape_cb_size, subvect_size, nb_subvect;
   split_cb_params *params;
   int N=2;
   int *best_index;
   float *best_dist;

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
   resp = PUSH(stack, shape_cb_size*subvect_size);
   t = PUSH(stack, nsf);
   r = PUSH(stack, nsf);
   e = PUSH(stack, nsf);
   E = PUSH(stack, shape_cb_size);
   ind = (int*)PUSH(stack, nb_subvect);

   tmp = PUSH(stack, 2*N*nsf);
   for (i=0;i<N;i++)
   {
      ot[i]=tmp;
      tmp += nsf;
      nt[i]=tmp;
      tmp += nsf;
   }

   best_index = (int*)PUSH(stack, N);
   best_dist = PUSH(stack, N);
   ndist = PUSH(stack, N);
   odist = PUSH(stack, N);
   
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

   e[0]=1;
   for (i=1;i<nsf;i++)
      e[i]=0;
   residue_zero(e, awk1, r, nsf, p);
   syn_filt_zero(r, ak, r, nsf, p);
   syn_filt_zero(r, awk2, r, nsf,p);
   
   /* Pre-compute codewords response and energy */
   for (i=0;i<shape_cb_size;i++)
   {
      float *res = resp+i*subvect_size;

      /* Compute codeword response */
      int k;
      for(j=0;j<subvect_size;j++)
         res[j]=0;
      for(j=0;j<subvect_size;j++)
      {
         for (k=j;k<subvect_size;k++)
            res[k]+=shape_cb[i*subvect_size+j]*r[k-j];
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
         vq_nbest(x, resp, subvect_size, shape_cb_size, E, N, best_index, best_dist);

         /*For all new n-bests*/
         for (k=0;k<N;k++)
         {
            float err=0;
            /*previous target*/
            for (m=0;m<nsf;m++)
               t[m]=ot[j][m];
            /*update target*/
            for (m=0;m<subvect_size;m++)
            {
               float g=shape_cb[best_index[k]*subvect_size+m];
               for (n=subvect_size*i+m,q=0;n<nsf;n++,q++)
                  t[n] -= g*r[q];
            }
            
            /*compute error (distance)*/
            err=odist[j];
            for (m=i*subvect_size;m<(i+1)*subvect_size;m++)
               err += t[m]*t[m];
            /*update n-best list*/
            if (err<ndist[N-1] || ndist[N-1]<-.5)
            {
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

      /*opdate old-new data*/
      for (j=0;j<N;j++)
         for (m=0;m<nsf;m++)
            ot[j][m]=nt[j][m];
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
      speex_bits_pack(bits,ind[i],params->shape_bits);
   }
   
   /* Put everything back together */
   for (i=0;i<nb_subvect;i++)
      for (j=0;j<subvect_size;j++)
         e[subvect_size*i+j]=shape_cb[ind[i]*subvect_size+j];

   /* Update excitation */
   for (j=0;j<nsf;j++)
      exc[j]+=e[j];
   
   /* Update target */
   residue_zero(e, awk1, r, nsf, p);
   syn_filt_zero(r, ak, r, nsf, p);
   syn_filt_zero(r, awk2, r, nsf,p);
   for (j=0;j<nsf;j++)
      target[j]-=r[j];

   POP(stack);
   POP(stack);
   POP(stack);
   POP(stack);
   POP(stack);
   POP(stack);
   POP(stack);
   POP(stack);
   POP(stack);
   POP(stack);
   POP(stack);
   POP(stack);
}





void split_cb_search_shape_sign(
float target[],			/* target vector */
float ak[],			/* LPCs for this subframe */
float awk1[],			/* Weighted LPCs for this subframe */
float awk2[],			/* Weighted LPCs for this subframe */
void *par,                      /* Codebook/search parameters*/
int   p,                        /* number of LPC coeffs */
int   nsf,                      /* number of samples in subframe */
float *exc,
SpeexBits *bits,
float *stack,
int   complexity
)
{
   int i,j;
   float *resp;
   float *t, *r, *e, *E;
   int *ind, *signs;
   float *shape_cb;
   int shape_cb_size, subvect_size, nb_subvect;
   split_cb_params *params;

   params = (split_cb_params *) par;
   subvect_size = params->subvect_size;
   nb_subvect = params->nb_subvect;
   shape_cb_size = 1<<params->shape_bits;
   shape_cb = params->shape_cb;
   resp = PUSH(stack, shape_cb_size*subvect_size);
   t = PUSH(stack, nsf);
   r = PUSH(stack, nsf);
   e = PUSH(stack, nsf);
   E = PUSH(stack, shape_cb_size);
   ind = (int*)PUSH(stack, nb_subvect);
   signs = (int*)PUSH(stack, nb_subvect);

   for (i=0;i<nsf;i++)
      t[i]=target[i];

   e[0]=1;
   for (i=1;i<nsf;i++)
      e[i]=0;
   residue_zero(e, awk1, r, nsf, p);
   syn_filt_zero(r, ak, r, nsf, p);
   syn_filt_zero(r, awk2, r, nsf,p);
   
   /* Pre-compute codewords response and energy */
   for (i=0;i<shape_cb_size;i++)
   {
      float *res = resp+i*subvect_size;

      /* Compute codeword response */
      int k;
      for(j=0;j<subvect_size;j++)
         res[j]=0;
      for(j=0;j<subvect_size;j++)
      {
         for (k=j;k<subvect_size;k++)
            res[k]+=shape_cb[i*subvect_size+j]*r[k-j];
      }
      E[i]=0;
      for(j=0;j<subvect_size;j++)
         E[i]+=res[j]*res[j];
   }

   for (i=0;i<nb_subvect;i++)
   {
      int best_index[2]={0,0}, k, m;
      float g, dist, best_dist[2]={-1,-1}, best_sign[2]={0,0};
      float *a, *x;
      float energy=0;
      x=t+subvect_size*i;

      for (k=0;k<subvect_size;k++)
         energy+=x[k]*x[k];
      /* Find best codeword for current sub-vector */
      for (j=0;j<shape_cb_size;j++)
      {
         int sign;
         dist=0;
         a=resp+j*subvect_size;
         dist=0;
         for (k=0;k<subvect_size;k++)
            dist -= 2*a[k]*x[k];
         if (dist > 0)
         {
            sign=1;
            dist =- dist;
         } else
            sign=0;
         dist += energy+E[j];
         if (dist<best_dist[0] || best_dist[0]<0)
         {
            best_dist[1]=best_dist[0];
            best_index[1]=best_index[0];
            best_sign[1]=best_sign[0];
            best_dist[0]=dist;
            best_index[0]=j;
            best_sign[0]=sign;
         } else if (dist<best_dist[1] || best_dist[1]<0)
         {
            best_dist[1]=dist;
            best_index[1]=j;
            best_sign[1]=sign;
         }
      }
      if (i<nb_subvect-1)
      {
         int nbest;
         float *tt, err[2];
         float best_score[2];
         tt=PUSH(stack,nsf);
         for (nbest=0;nbest<2;nbest++)
         {
            float s=1;
            if (best_sign[nbest])
               s=-1;
            for (j=0;j<nsf;j++)
               tt[j]=t[j];
            for (j=0;j<subvect_size;j++)
            {
               g=s*shape_cb[best_index[nbest]*subvect_size+j];
               for (k=subvect_size*i+j,m=0;k<nsf;k++,m++)
                  tt[k] -= g*r[m];
            }
            
            {
               int best_index2=0, best_sign2=0, sign2;
               float  best_dist2=0;
               x=t+subvect_size*(i+1);
               for (j=0;j<shape_cb_size;j++)
               {
                  a=resp+j*subvect_size;
                  dist = 0;
                  for (k=0;k<subvect_size;k++)
                     dist -= 2*a[k]*x[k];
                  if (dist > 0)
                  {
                     sign2=1;
                     dist =- dist;
                  } else
                     sign2=0;
                  dist += energy+E[j];
                  if (dist<best_dist2 || j==0)
                  {
                     best_dist2=dist;
                     best_index2=j;
                     best_sign2=sign2;
                  }
               }
               s=1;
               if (best_sign2)
                  s=-1;
               /*int i2=vq_index(&tt[subvect_size*(i+1)], resp, subvect_size, shape_cb_size);*/
               
               for (j=0;j<subvect_size;j++)
               {
                  g=s*shape_cb[best_index2*subvect_size+j];
                  for (k=subvect_size*(i+1)+j,m=0;k<nsf;k++,m++)
                     tt[k] -= g*r[m];
               }
            }

            err[nbest]=0;
            for (j=subvect_size*i;j<subvect_size*(i+2);j++)
               err[nbest]-=tt[j]*tt[j];
            
            best_score[nbest]=err[nbest];
         }

         if (best_score[1]>best_score[0])
         {
            best_sign[0]=best_sign[1];
            best_index[0]=best_index[1];
            best_score[0]=best_score[1];
         }
         POP(stack);

      }

      ind[i]=best_index[0];
      signs[i] = best_sign[0];

      /*printf ("best index: %d/%d\n", best_index, shape_cb_size);*/
      speex_bits_pack(bits,signs[i],1);
      speex_bits_pack(bits,ind[i],params->shape_bits);

      /* Update target for next subvector */
      for (j=0;j<subvect_size;j++)
      {
         g=shape_cb[ind[i]*subvect_size+j];
         if (signs[i])
            g=-g;
         for (k=subvect_size*i+j,m=0;k<nsf;k++,m++)
            t[k] -= g*r[m];
      }
   }
   
   /* Put everything back together */
   for (i=0;i<nb_subvect;i++)
   {
      float s=1;
      if (signs[i])
         s=-1;
      for (j=0;j<subvect_size;j++)
         e[subvect_size*i+j]=s*shape_cb[ind[i]*subvect_size+j];
   }
   /* Update excitation */
   for (j=0;j<nsf;j++)
      exc[j]+=e[j];
   
   /* Update target */
   residue_zero(e, awk1, r, nsf, p);
   syn_filt_zero(r, ak, r, nsf, p);
   syn_filt_zero(r, awk2, r, nsf,p);
   for (j=0;j<nsf;j++)
      target[j]-=r[j];

   
   POP(stack);
   POP(stack);
   POP(stack);
   POP(stack);
   POP(stack);
   POP(stack);
   POP(stack);
}


void split_cb_nogain_unquant(
float *exc,
void *par,                      /* non-overlapping codebook */
int   nsf,                      /* number of samples in subframe */
SpeexBits *bits,
float *stack
)
{
   int i,j;
   int *ind;
   float *shape_cb;
   int shape_cb_size, subvect_size, nb_subvect;
   split_cb_params *params;

   params = (split_cb_params *) par;
   subvect_size = params->subvect_size;
   nb_subvect = params->nb_subvect;
   shape_cb_size = 1<<params->shape_bits;
   shape_cb = params->shape_cb;
   
   ind = (int*)PUSH(stack, nb_subvect);

   /* Decode codewords and gains */
   for (i=0;i<nb_subvect;i++)
      ind[i] = speex_bits_unpack_unsigned(bits, params->shape_bits);

   /* Compute decoded excitation */
   for (i=0;i<nb_subvect;i++)
      for (j=0;j<subvect_size;j++)
         exc[subvect_size*i+j]+=shape_cb[ind[i]*subvect_size+j];

   POP(stack);
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

   params = (split_cb_params *) par;
   subvect_size = params->subvect_size;
   nb_subvect = params->nb_subvect;
   shape_cb_size = 1<<params->shape_bits;
   shape_cb = params->shape_cb;
   
   ind = (int*)PUSH(stack, nb_subvect);
   signs = (int*)PUSH(stack, nb_subvect);

   /* Decode codewords and gains */
   for (i=0;i<nb_subvect;i++)
   {
      signs[i] = speex_bits_unpack_unsigned(bits, 1);
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
   POP(stack);
   POP(stack);
}
