/*-----------------------------------------------------------------------*\

    FILE........: GAINSHAPE.C
    TYPE........: C Module
    AUTHOR......: David Rowe
    COMPANY.....: Voicetronix
    DATE CREATED: 19/2/02

    General gain-shape codebook search.

\*-----------------------------------------------------------------------*/

/* Modified by Jean-Marc Valin 2002

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

#define EXC_CB_SIZE 128
#define min(a,b) ((a) < (b) ? (a) : (b))
extern float exc_gains_table[];
extern float exc_table[];

/*---------------------------------------------------------------------------*\
                                                                             
 void overlap_cb_search()							      
									      
 Searches a gain/shape codebook consisting of overlapping entries for the    
 closest vector to the target.  Gives identical results to search() above   
 buts uses fast end correction algorithm for the synthesis of response       
 vectors.								      
                                                                             
\*---------------------------------------------------------------------------*/

float overlap_cb_search(
float target[],			/* target vector */
float ak[],			/* LPCs for this subframe */
float awk1[],			/* Weighted LPCs for this subframe */
float awk2[],			/* Weighted LPCs for this subframe */
float codebook[],		/* overlapping codebook */
int   entries,			/* number of overlapping entries to search */
float *gain,			/* gain of optimum entry */
int   *index,			/* index of optimum entry */
int   p,                        /* number of LPC coeffs */
int   nsf                       /* number of samples in subframe */
)
{
  float *resp;		        /* zero state response to current entry */
  float *h;		        /* impulse response of synthesis filter */
  float *impulse;		/* excitation vector containing one impulse */
  float d,e,g,score;		/* codebook searching variables */
  float bscore;			/* score of "best" vector so far */
  int i,k;			/* loop variables */

  /* Initialise */
  
  resp = (float*)malloc(sizeof(float)*nsf);
  h = (float*)malloc(sizeof(float)*nsf);
  impulse = (float*)malloc(sizeof(float)*nsf);

  for(i=0; i<nsf; i++)
    impulse[i] = 0.0;
   
  *gain = 0.0;
  *index = 0;
  bscore = 0.0;
  impulse[0] = 1.0;

  /* Calculate impulse response of  A(z/g2) / ( A(z)*(z/g1) ) */
  residue_zero(impulse, awk1, h, nsf, p);
  syn_filt_zero(h, ak, h, nsf, p);
  syn_filt_zero(h, awk2, h, nsf,p);
  
  /* Calculate codebook zero-response */
  residue_zero(&codebook[entries-1],awk1,resp,nsf,p);
  syn_filt_zero(resp,ak,resp,nsf,p);
  syn_filt_zero(resp,awk2,resp,nsf,p);
    
  /* Search codebook backwards using end correction for synthesis */
  
  for(k=entries-1; k>=0; k--) {

    d = 0.0; e = 0.0;
    for(i=0; i<nsf; i++) {
      d += target[i]*resp[i];
      e += resp[i]*resp[i];
    }
    g = d/(e+.1);
    score = g*d;
    /*printf ("score: %f %f %f %f\n", target[0],d,e,score);*/
    if (score >= bscore) {
      bscore = score;
      *gain = g;
      *index = k;
    }
    
    /* Synthesise next entry */
    
    if (k) {
      for(i=nsf-1; i>=1; i--)
        resp[i] = resp[i-1] + codebook[k-1]*h[i];
      resp[0] = codebook[k-1]*h[0];
    }
  }

  free(resp);
  free(h);
  free(impulse);
  return bscore;
}


split_cb_params split_cb_nb = {
   8,               /*subvect_size*/
   5,               /*nb_subvect*/
   exc_table,       /*shape_cb*/
   7,               /*shape_bits*/
   exc_gains_table, /*gain_cb*/
   8                /*gain_bits*/
};


void split_cb_search(
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
   int i,j;
   float *resp, *E, *Ee;
   float *t, *r, *e;
   float *gains;
   int *ind;
   float *shape_cb, *gain_cb;
   int shape_cb_size, gain_cb_size, subvect_size, nb_subvect;
   split_cb_params *params;

   params = (split_cb_params *) par;
   subvect_size = params->subvect_size;
   nb_subvect = params->nb_subvect;
   shape_cb_size = 1<<params->shape_bits;
   shape_cb = params->shape_cb;
   gain_cb_size = 1<<params->gain_bits;
   gain_cb = params->gain_cb;
   resp = PUSH(stack, shape_cb_size*subvect_size);
   E = PUSH(stack, shape_cb_size);
   Ee = PUSH(stack, shape_cb_size);
   t = PUSH(stack, nsf);
   r = PUSH(stack, nsf);
   e = PUSH(stack, nsf);
   gains = PUSH(stack, nb_subvect);
   ind = (int*)PUSH(stack, nb_subvect);
   

   for (i=0;i<nsf;i++)
      t[i]=target[i];
   for (i=0;i<shape_cb_size;i++)
   {
      float *res = resp+i*subvect_size;
      residue_zero(shape_cb+i*subvect_size, awk1, res, subvect_size, p);
      syn_filt_zero(res, ak, res, subvect_size, p);
      syn_filt_zero(res, awk2, res, subvect_size,p);
      E[i]=0;
      for(j=0;j<subvect_size;j++)
         E[i]+=res[j]*res[j];
      Ee[i]=0;
      for(j=0;j<subvect_size;j++)
         Ee[i]+=shape_cb[i*subvect_size+j]*shape_cb[i*subvect_size+j];
      
   }
   for (i=0;i<nb_subvect;i++)
   {
      int best_index=0;
      float g, corr, best_gain=0, score, best_score=-1;
      for (j=0;j<shape_cb_size;j++)
      {
         corr=xcorr(resp+j*subvect_size,t+subvect_size*i,subvect_size);
         score=corr*corr/(.001+E[j]);
         g = corr/(.001+E[j]);
         if (score>best_score)
         {
            best_index=j;
            best_score=score;
            best_gain=corr/(.001+E[j]);
         }
      }
      frame_bits_pack(bits,best_index,params->shape_bits);
      if (best_gain>0)
         frame_bits_pack(bits,0,1);
      else
          frame_bits_pack(bits,1,1);        
      ind[i]=best_index;
      gains[i]=best_gain*Ee[ind[i]];

      for (j=0;j<nsf;j++)
         e[j]=0;
      for (j=0;j<subvect_size;j++)
         e[subvect_size*i+j]=best_gain*shape_cb[best_index*subvect_size+j];
      residue_zero(e, awk1, r, nsf, p);
      syn_filt_zero(r, ak, r, nsf, p);
      syn_filt_zero(r, awk2, r, nsf,p);
      for (j=0;j<nsf;j++)
         t[j]-=r[j];
   }

   {
      int best_vq_index=0, max_index;
      float max_gain=0, log_max, min_dist=0, *sign;

      sign = PUSH(stack, nb_subvect);
      for (i=0;i<nb_subvect;i++)
      {
         if (gains[i]<0)
         {
            gains[i]=-gains[i];
            sign[i]=-1;
         } else {
            sign[i]=1;
         }
      }
      for (i=0;i<nb_subvect;i++)
         if (gains[i]>max_gain)
            max_gain=gains[i];
      log_max=log(max_gain+1);
      max_index = (int)(floor(.5+log_max-3));
      if (max_index>7)
         max_index=7;
      if (max_index<0)
         max_index=0;
      max_gain=1/exp(max_index+3.0);
      for (i=0;i<nb_subvect;i++)
        gains[i]*=max_gain;
      frame_bits_pack(bits,max_index,3);

      /*Vector quantize gains[i]*/
      best_vq_index = vq_index(gains, gain_cb, nb_subvect, gain_cb_size);
      frame_bits_pack(bits,best_vq_index,params->gain_bits);

      printf ("best_gains_vq_index %d %f %d\n", best_vq_index, min_dist, max_index);

#if 1 /* If 0, the gains are not quantized */
      for (i=0;i<nb_subvect;i++)
         gains[i]= sign[i]*gain_cb[best_vq_index*nb_subvect+i]/max_gain/(Ee[ind[i]]+.001);
#else 
      for (i=0;i<nb_subvect;i++)
         gains[i]= sign[i]*gains[i]/max_gain/(Ee[ind[i]]+.001);
#endif  
    
      for (i=0;i<nb_subvect;i++)
         for (j=0;j<subvect_size;j++)
            exc[subvect_size*i+j]+=gains[i]*shape_cb[ind[i]*subvect_size+j];

      POP(stack);
   }

   /*TODO: Perform joint optimization of gains*/
   
   for (i=0;i<nsf;i++)
      target[i]=t[i];

   POP(stack);
   POP(stack);
   POP(stack);
   POP(stack);
   POP(stack);
   POP(stack);
   POP(stack);
   POP(stack);
}

void split_cb_unquant(
float *exc,
float codebook[][8],		/* non-overlapping codebook */
int   nsf,                      /* number of samples in subframe */
FrameBits *bits
)
{
   int i,j;
   int ind[5];
   float gains[5];
   float sign[5];
   int max_gain_ind, vq_gain_ind;
   float max_gain, Ee[5];
   for (i=0;i<5;i++)
   {
      ind[i] = frame_bits_unpack_unsigned(bits, 7);
      if (frame_bits_unpack_unsigned(bits, 1))
         sign[i]=-1;
      else
         sign[i]=1;
      Ee[i]=.001;
      for (j=0;j<8;j++)
         Ee[i]+=codebook[ind[i]][j]*codebook[ind[i]][j];
   }
   max_gain_ind = frame_bits_unpack_unsigned(bits, 3);
   vq_gain_ind = frame_bits_unpack_unsigned(bits, 8);
   printf ("unquant gains ind: %d %d\n", max_gain_ind, vq_gain_ind);

   max_gain=exp(max_gain_ind+3.0);
   for (i=0;i<5;i++)
      gains[i] = sign[i]*exc_gains_table[vq_gain_ind*5+i]*max_gain/Ee[i];
   
   printf ("unquant gains: ");
   for (i=0;i<5;i++)
      printf ("%f ", gains[i]);
   printf ("\n");

   for (i=0;i<5;i++)
      for (j=0;j<8;j++)
         exc[8*i+j]+=gains[i]*codebook[ind[i]][j];
   
}
