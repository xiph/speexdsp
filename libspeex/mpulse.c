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
#include <stdlib.h>
#include "filters.h"
#include <math.h>

#define MAX_PULSE 30
#define MAX_POS  100

int porder(int *p, int *s, int *o, int len)
{
   int i,j, bit1, nb_uniq=0;
   int *st, *en;
   int rep[MAX_POS];
   int uniq[MAX_PULSE];
   int n;
   /*Stupid bubble sort but for small N, we don't care!*/
   for (i=0;i<MAX_POS;i++)
      rep[i]=0;
   for (i=0;i<len;i++)
   {
      for (j=i+1;j<len;j++)
      {
         if (p[i]>p[j])
         {
            int tmp;
            tmp=p[j];
            p[j]=p[i];
            p[i]=tmp;
            tmp=s[j];
            s[j]=s[i];
            s[i]=tmp;
         }
      }
   }
   printf ("quant_pulse\n");
   for (i=0;i<len;i++)
      printf ("%d ", p[i]);
   printf ("\n");
   for (i=0;i<len;i++)
      printf ("%d ", s[i]);
   printf ("\n");

   for (i=0;i<len;i++)
   {
      rep[p[i]]++;
      if (i==0 || p[i]!=p[i-1])
      {
         uniq[nb_uniq]=p[i];
         s[nb_uniq]=s[i];
         nb_uniq++;
      }
   }
   st=uniq;
   en=&uniq[nb_uniq-1];


   bit1=s[0];
   n=0;
   for (i=0;i<nb_uniq;i++)
   {
      int next;
      if (i==nb_uniq-1)
      {
         next=*st;
         for (j=0;j<rep[next];j++)
         {
            o[n]=next;
            n++;
         }
      } else {
         if (s[i+1])
         {
            next=*en;
            for (j=0;j<rep[next];j++)
            {
               o[n]=next;
               n++;
            }
            en--;
         } else {
            next=*st;
            for (j=0;j<rep[next];j++)
            {
               o[n]=next;
               n++;
            }
            st++;
         }
      }
   }
   for (i=0;i<len;i++)
      printf ("%d ", o[i]);
   printf ("\n");

   return s[0];
}

void rorder(int *p, int *s, int *o, int bit, int len)
{
   int i,j,nb_uniq=0;
   int *st, *en;
   int rep[MAX_POS];
   int uniq[MAX_PULSE];
   int ss[MAX_PULSE];
   int n;
   /*Stupid bubble sort but for small N, we don't care!*/
   for (i=0;i<len;i++)
      o[i]=p[i];
   for (i=0;i<len;i++)
   {
      for (j=i+1;j<len;j++)
      {
         if (o[i]>o[j])
         {
            int tmp;
            tmp=o[j];
            o[j]=o[i];
            o[i]=tmp;
         }
      }
   }

   for (i=0;i<len;i++)  
   {
      rep[p[i]]++;
      if (i==0 || o[i]!=o[i-1])
      {
         uniq[nb_uniq]=o[i];
         s[nb_uniq]=s[i];
         nb_uniq++;
      }
   }
   st=uniq;
   en=&uniq[nb_uniq-1];

   ss[0]=bit;
   n=1;

   printf ("unquant_pulse\n");
   for (i=0;i<len;i++)
      printf ("%d ", o[i]);
   printf ("\n");

   for (i=1;i<len;i++)
   {
      if (i>1&&p[i-1]==p[i-2])
         continue;
      if (p[i-1]==*st)
      {
         ss[n++]=0;
         st++;
      } else if (p[i-1]==*en)
      {
         ss[n++]=1;
         en--;
      } else 
      {
         fprintf (stderr, "ERROR in decoding signs\n");
         exit(1);
      }
   }
   
   n=0;
   for (i=0;i<len;i++)
   {
      s[i]=ss[n];
      if (i<len&&o[i]!=o[i+1])
         n++;
   }

   for (i=0;i<len;i++)
      printf ("%d ", s[i]);
   printf ("\n");

}


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
   int i,j, nb_pulse;
   float *resp, *t, *e, *pulses;
   float te=0,ee=0;
   float g;
   int nb_tracks, track_ind_bits;
   int *tracks, *signs, *tr, *nb;
   mpulse_params *params;
   int pulses_per_track;
   params = (mpulse_params *) par;

   nb_pulse=params->nb_pulse;
   nb_tracks=params->nb_tracks;
   pulses_per_track=nb_pulse/nb_tracks;
   track_ind_bits=params->track_ind_bits;

   tracks = (int*)PUSH(stack,nb_pulse);
   signs = (int*)PUSH(stack,nb_pulse);
   tr = (int*)PUSH(stack,pulses_per_track);
   nb = (int*)PUSH(stack,nb_tracks);

   resp=PUSH(stack, nsf);
   t=PUSH(stack, nsf);
   e=PUSH(stack, nsf);
   pulses=PUSH(stack, nsf);
   
   syn_filt_zero(target, awk1, e, nsf, p);
   residue_zero(e, ak, e, nsf, p);
   residue_zero(e, awk2, e, nsf, p);
   for (i=0;i<nsf;i++)
   {
      pulses[i]=0;
      te+=target[i]*target[i];
      ee+=e[i]*e[i];
   }
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

   for (i=0;i<nb_tracks;i++)
      nb[i]=0;

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
         if ((i%nb_tracks) != (j%nb_tracks))
           continue;
         /*Try for positive sign*/
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
         /*Try again for negative sign*/
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
      }
      printf ("best pulse: %d %d %f %f %f %f\n", i, best_ind, best_gain, te, ee, g);
      /*Remove pulse contribution from target*/
      for (j=best_ind;j<nsf;j++)
         t[j] -= best_gain * resp[j-best_ind];
      e[best_ind]+=best_gain;
      if (best_gain>0)
         pulses[best_ind]+=1;
      else
         pulses[best_ind]-=1;
      {
         int t=best_ind%nb_tracks;
         tracks[t*pulses_per_track+nb[t]] = best_ind/nb_tracks;
         signs[t*pulses_per_track+nb[t]]  = best_gain >= 0 ? 0 : 1;
         nb[t]++;
      }
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
      frame_bits_pack(bits,quant_gain,7);
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
   
   for (i=0;i<nb_tracks;i++)
   {
      int bit1, ind=0;
      bit1=porder(tracks+i*pulses_per_track, signs+i*pulses_per_track,tr,pulses_per_track);
      frame_bits_pack(bits,bit1,1);
      for (j=0;j<pulses_per_track;j++)
      {
         ind*=nsf/nb_tracks;
         ind+=tr[j];
         /*printf ("%d ", ind);*/
      }
      
      frame_bits_pack(bits,ind,track_ind_bits);

      /*printf ("track %d %d:", i, ind);
      for (j=0;j<pulses_per_track;j++)
        printf ("%d ", tr[j]);
        printf ("\n");*/
   }
   POP(stack);
   POP(stack);
   POP(stack);
   POP(stack);
   POP(stack);
   POP(stack);
   POP(stack);
   POP(stack);
}


void mpulse_unquant(
float *exc,
void *par,                      /* non-overlapping codebook */
int   nsf,                      /* number of samples in subframe */
FrameBits *bits,
float *stack
)
{
   int i,j, bit1, nb_pulse, quant_gain;
   float g;
   int nb_tracks, track_ind_bits;
   int *track, *signs, *tr;
   mpulse_params *params;
   int pulses_per_track;
   params = (mpulse_params *) par;

   nb_pulse=params->nb_pulse;
   nb_tracks=params->nb_tracks;
   pulses_per_track=nb_pulse/nb_tracks;
   track_ind_bits=params->track_ind_bits;

   track = (int*)PUSH(stack,pulses_per_track);
   signs = (int*)PUSH(stack,pulses_per_track);
   tr = (int*)PUSH(stack,pulses_per_track);
   
   quant_gain=frame_bits_unpack_unsigned(bits, 7);
   g=exp((quant_gain/8.0)+1);
   /*Removes glitches when energy is near-zero*/
   if (g<3)
      g=0;
   for (i=0;i<nb_tracks;i++)
   {
      int ind;
      int max_val=nsf/nb_tracks;
      bit1=frame_bits_unpack_unsigned(bits, 1);
      ind = frame_bits_unpack_unsigned(bits,track_ind_bits);
      /*printf ("unquant ind = %d\n", ind);*/
      for (j=0;j<pulses_per_track;j++)
      {
         track[pulses_per_track-1-j]=ind%max_val;
         ind /= max_val;
      }
      rorder(track, signs, tr, bit1, pulses_per_track);
      for (j=0;j<pulses_per_track;j++)
      {
         exc[tr[j]*nb_tracks+i] += signs[j] ? -g : g;
      }
   }

   POP(stack);
   POP(stack);
   POP(stack);
}
