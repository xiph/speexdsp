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


