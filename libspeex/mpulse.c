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
#include <math.h>

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
   float *resp, *t, *e, *pulses;
   float te=0,ee=0;
   float g;
   int tracks[4]={0,0,0,0};
   float *sign;
   int pulses_per_track=nb_pulse/4;
   resp=PUSH(stack, nsf);
   t=PUSH(stack, nsf);
   e=PUSH(stack, nsf);
   pulses=PUSH(stack, nsf);
   sign=PUSH(stack, nsf);
   
   syn_filt_zero(target, awk1, e, nsf, p);
   residue_zero(e, ak, e, nsf, p);
   residue_zero(e, awk2, e, nsf, p);
   for (i=0;i<nsf;i++)
   {
      pulses[i]=0;
      if (e[i]>0)
         sign[i]=1;
      else 
         sign[i]=-1;
      te+=target[i]*target[i];
      ee+=e[i]*e[i];
   }
   /*g=1.9154e-05*te+2.3396e-05*ee;*/

   g=2.2/sqrt(nb_pulse)*exp(0.18163*log(te+1)+0.17293*log(ee+1));
   
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
      float best_score=1e30, best_gain=0;
      int best_ind=0;
      /*For all positions*/
      for (j=0;j<nsf;j++)
      {
         int k;
         float dist=0;
         /*Constrain search in alternating tracks*/
         /*if ((i&3) != (j&3))
           continue;*/
         if (tracks[j&3]==pulses_per_track)
            continue;
#if 0
         for (k=0;k<j;k++)
            dist+=t[k]*t[k];
         for (k=0;k<nsf-j;k++)
            dist+=(t[k+j]-sign[j]*g*resp[k])*(t[k+j]-sign[j]*g*resp[k]);
         if (dist<best_score || j==0)
         {
            best_score=dist;
            best_gain=sign[j]*g;
            best_ind=j;
         }         
#else
         for (k=0;k<j;k++)
            dist+=t[k]*t[k];
         for (k=0;k<nsf-j;k++)
            dist+=(t[k+j]-g*resp[k])*(t[k+j]-g*resp[k]);
         if (dist<best_score || j==0)
         {
            best_score=dist;
            best_gain=g;
            best_ind=j;
         }
         /*printf ("dist: %f\n", dist);*/
         dist=0;
         for (k=0;k<j;k++)
            dist+=t[k]*t[k];
         for (k=0;k<nsf-j;k++)
            dist+=(t[k+j]+g*resp[k])*(t[k+j]+g*resp[k]);
         if (dist<best_score || j==0)
         {
            best_score=dist;
            best_gain=-g;
            best_ind=j;
         }
#endif
      }
      tracks[best_ind&3]++;
      printf ("best pulse: %d %d %f %f %f %f\n", i, best_ind, best_gain, te, ee, g);
      /*Remove pulse contribution from target*/
      for (j=best_ind;j<nsf;j++)
         t[j] -= best_gain * resp[j-best_ind];
      e[best_ind]+=best_gain;
      if (best_gain>0)
         pulses[best_ind]+=1;
      else
         pulses[best_ind]-=1;
   }
   
   /*Global gain re-estimation*/
   if (1) {
      float f;
      int quant_gain;
      residue_zero(e, awk1, resp, nsf, p);
      syn_filt_zero(resp, ak, resp, nsf, p);
      syn_filt_zero(resp, awk2, resp, nsf, p);

      f=((.1+(xcorr(resp,target,nsf)))/(.1+xcorr(resp,resp,nsf)));
      /*for (i=0;i<nsf;i++)
        e[i]*=f;*/
      g *= f;
      if (g<0)
         g=0;
      
      quant_gain=(int)floor(.5+8*(log(1+fabs(g))-1));
      if (quant_gain<0)
         quant_gain=0;
      if (quant_gain>127)
         quant_gain=127;
      g=exp((quant_gain/8.0)+1);
      
      for (i=0;i<nsf;i++)
         e[i]=g*pulses[i];
      printf ("global gain = %f\n", g);
      for (i=0;i<nsf;i++)
         t[i]=target[i]-f*resp[i];

   }
   for (i=0;i<nsf;i++)
      printf ("%f ", e[i]);
   printf ("\n");
   for (i=0;i<nsf;i++)
      exc[i]+=e[i];
   for (i=0;i<nsf;i++)
      target[i]=t[i];

   POP(stack);
   POP(stack);
   POP(stack);
   POP(stack);
}
