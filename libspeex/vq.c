/* Copyright (C) 2002 Jean-Marc Valin
   File: vq.c
   Vector quantization

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

#include "vq.h"

/*Finds the index of the entry in a codebook that best matches the input*/
int vq_index(float *in, float *codebook, int len, int entries)
{
   int i,j;
   float min_dist=0;
   int best_index=0;
   for (i=0;i<entries;i++)
   {
      float dist=0;
      for (j=0;j<len;j++)
      {
         float tmp = in[j]-*codebook++;
         dist += tmp*tmp;
      }
      if (i==0 || dist<min_dist)
      {
         min_dist=dist;
         best_index=i;
      }
   }
   return best_index;
}


/*Finds the indices of the n-best entries in a codebook*/
void vq_nbest(float *in, float *codebook, int len, int entries, float *E, int N, int *nbest, float *best_dist)
{
   int i,j,k;
   for (i=0;i<entries;i++)
   {
      float dist=.5*E[i];
      for (j=0;j<len;j++)
         dist -= in[j]**codebook++;
      if (i<N || dist<best_dist[N-1])
      {

         for (j=0;j<N;j++)
         {
            if (j >= i || dist < best_dist[j])
            {
               for (k=N-1;k>j;k--)
               {
                  best_dist[k]=best_dist[k-1];
                  nbest[k] = nbest[k-1];
               }
               best_dist[j]=dist;
               nbest[j]=i;
               break;
            }
         }
      }
   }
}

#ifdef _USE_SSE

void vq_nbest8(float *in, float *codebook, int len, int entries, float *E, int N, int *nbest, float *best_dist)
{
   int i,j,k;
   for (i=0;i<entries;i++)
   {
      float dist=-.5*E[i];
      /*for (j=0;j<8;j++)
        dist += in[j]**codebook++;*/
      __asm__ __volatile__ (
      "
      movups (%0), %%xmm0
      movups (%1), %%xmm1
      mulps %%xmm1, %%xmm0
      movups 16(%0), %%xmm2
      movups 16(%1), %%xmm3
      mulps %%xmm3, %%xmm2
      movss (%2), %%xmm4
      addps %%xmm2, %%xmm0
      movhlps %%xmm0, %%xmm1
      addps %%xmm1, %%xmm0
      movaps %%xmm0, %%xmm1
      shufps $0x01, %%xmm1, %%xmm1
      addss %%xmm1, %%xmm0
      addss %%xmm0, %%xmm4
      movss %%xmm4, (%2)

      " : : "r" (in), "r" (codebook), "r" (&dist) : "memory");
      codebook+=8;
      if (i<N || dist>best_dist[N-1])
      {

         for (j=0;j<N;j++)
         {
            if (j >= i || dist > best_dist[j])
            {
               for (k=N-1;k>j;k--)
               {
                  best_dist[k]=best_dist[k-1];
                  nbest[k] = nbest[k-1];
               }
               best_dist[j]=dist;
               nbest[j]=i;
               break;
            }
         }
      }
   }
}


void vq_nbest5(float *in, float *codebook, int len, int entries, float *E, int N, int *nbest, float *best_dist)
{
   int i,j,k;
   float scores[256];
   for (i=0;i<entries;i++)
   {
      float dist=-.5*E[i];
      /*for (j=0;j<8;j++)
        dist += in[j]**codebook++;*/
      __asm__ __volatile__ (
      "
      movups (%0), %%xmm0
      movups (%1), %%xmm1
      mulps %%xmm1, %%xmm0
      movss 16(%0), %%xmm2
      movss 16(%1), %%xmm3
      mulss %%xmm3, %%xmm2
      movss (%2), %%xmm4
      addss %%xmm2, %%xmm0
      movhlps %%xmm0, %%xmm1
      addps %%xmm1, %%xmm0
      movaps %%xmm0, %%xmm1
      shufps $0x01, %%xmm1, %%xmm1
      addss %%xmm1, %%xmm0
      addss %%xmm0, %%xmm4
      movss %%xmm4, (%2)

      " : : "r" (in), "r" (codebook), "r" (&dist) : "memory");
      codebook+=5;
      scores[i]=dist;
   }
   for (i=0;i<entries;i++)
   {
      float dist=scores[i];
      __asm__ __volatile__ ("tata2:");
      if (i<N || dist>best_dist[N-1])
      {

         for (j=0;j<N;j++)
         {
            if (j >= i || dist > best_dist[j])
            {
               for (k=N-1;k>j;k--)
               {
                  best_dist[k]=best_dist[k-1];
                  nbest[k] = nbest[k-1];
               }
               best_dist[j]=dist;
               nbest[j]=i;
               break;
            }
         }
      }
   }
}


#endif
