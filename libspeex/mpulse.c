/* Copyright (C) 2002 Jean-Marc Valin 
   File: mpulse.c

   Multi-pulse code

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

#include "mpulse.h"
#include "stack_alloc.h"
#include <stdio.h>
#include "filters.h"

void mpulse_search(
float target[],			/* target vector */
float ak[],			/* LPCs for this subframe */
float awk1[],			/* Weighted LPCs for this subframe */
float awk2[],			/* Weighted LPCs for this subframe */
void *par,                      /* Codebook/search parameters*/
int   p,                        /* number of LPC coeffs */
int   nsf,                      /* number of samples in subframe */
float *exc,
FrameBits *bits,
float *stack
)
{
   int i,j, nb_pulse=12;
   float *resp, *t, *e;
   
   resp=PUSH(stack, nsf);
   t=PUSH(stack, nsf);
   e=PUSH(stack, nsf);

   e[0]=1;
   for (i=1;i<nsf;i++)
      e[i]=0;

   residue_zero(e, awk1, resp, nsf, p);
   syn_filt_zero(resp, ak, resp, nsf, p);
   syn_filt_zero(resp, awk2, resp, nsf, p);
   
   for (i=0;i<nsf;i++)
      e[i]=0;

   for (i=0;i<nsf;i++)
      t[i]=target[i];

   /*For all pulses*/
   for (i=0;i<nb_pulse;i++)
   {
      float best_score=0, best_gain=0;
      int best_ind=0;
      /*For all positions*/
      for (j=0;j<nsf;j++)
      {
         float corr, energy, score;
         corr=xcorr(t+j,resp,nsf-j);
         /*This can be pre-calculated*/
         energy=xcorr(resp,resp,nsf-j);
         score=corr*corr/energy;
         if (score>best_score)
         {
            best_score=score;
            best_ind=j;
            best_gain=corr/energy;
         }
      }
    
      printf ("best pulse: %d %f\n", best_ind, best_gain);
      /*Remove pulse contribution from target*/
      for (j=best_ind;j<nsf;j++)
         t[j] -= best_gain * resp[j-best_ind];
      e[best_ind]+=best_gain;
   }
   
   for (i=0;i<nsf;i++)
      printf ("%f ", e[i]);
   printf ("\n");
   for (i=0;i<nsf;i++)
      exc[i]+=e[i];
   for (i=0;i<nsf;i++)
      target[i]+=t[i];

   POP(stack);
   POP(stack);
   POP(stack);
}
