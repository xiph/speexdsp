/* Original copyright */
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
#include "matrix.h"

static float scal_gains4[16] = {
   0.27713,
   0.49282,
   0.69570,
   0.90786,
   1.14235,
   1.42798,
   1.80756,
   2.42801
};

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

  /* Calculate impulse response of  A(z/g1) / ( A(z)*(z/g2) ) */
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
   int i,j, id;
   float *resp, *E, q;
   float *t, *r, *e, *tresp;
   float *gains;
   int *ind;
   float *shape_cb;
   int shape_cb_size, subvect_size, nb_subvect;
   float exc_energy=0;
   split_cb_params *params;

   params = (split_cb_params *) par;
   subvect_size = params->subvect_size;
   nb_subvect = params->nb_subvect;
   shape_cb_size = 1<<params->shape_bits;
   shape_cb = params->shape_cb;
   resp = PUSH(stack, shape_cb_size*subvect_size);
   tresp = PUSH(stack, shape_cb_size*nsf);
   E = PUSH(stack, shape_cb_size);
   t = PUSH(stack, nsf);
   r = PUSH(stack, nsf);
   e = PUSH(stack, nsf);
   gains = PUSH(stack, nb_subvect);
   ind = (int*)PUSH(stack, nb_subvect);

   /* Compute energy of the "real excitation" */
   syn_filt_zero(target, awk1, e, nsf, p);
   residue_zero(e, ak, e, nsf, p);
   residue_zero(e, awk2, e, nsf, p);
   for (i=0;i<nsf;i++)
      exc_energy += e[i]*e[i];
   exc_energy=sqrt(exc_energy/nb_subvect);

   /* Quantize global ("average") gain */
   q=log(exc_energy+.1);
   q=floor(.5+2*(q-2));
   if (q<0)
      q=0;
   if (q>15)
      q=15;
   id = (int)q;
   frame_bits_pack(bits, id, 4);
   exc_energy=exp(.5*q+2);


   for (i=0;i<nsf;i++)
      t[i]=target[i];

   /* Pre-compute codewords response and energy */
   for (i=0;i<shape_cb_size;i++)
   {
      float *res = resp+i*subvect_size;

      /* Compute codeword response */
      residue_zero(shape_cb+i*subvect_size, awk1, res, subvect_size, p);
      syn_filt_zero(res, ak, res, subvect_size, p);
      syn_filt_zero(res, awk2, res, subvect_size,p);

      /* Compute energy of codeword response */
      E[i]=0;
      for(j=0;j<subvect_size;j++)
         E[i]+=res[j]*res[j];
      
   }

   for (i=0;i<nb_subvect;i++)
   {
      int best_index=0;
      float g, corr, best_gain=0, score, best_score=-1;
      /* Find best codeword for current sub-vector */
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
      
      /* Quantize gain */
      {
         int s=0, best_id;
         best_gain /= .01+exc_energy;
         if (best_gain<0)
         {
            best_gain=-best_gain;
            s=1;
         }

         /* Find gain index (it's a scalar but we use the VQ code anyway)*/
         best_id = vq_index(&best_gain, scal_gains4, 1, 8);

         best_gain=scal_gains4[best_id];
         /*printf ("gain_quant: %f %d %f\n", best_gain, best_id, scal_gains4[best_id]);*/
         if (s)
            best_gain=-best_gain;
         best_gain *= exc_energy;
         frame_bits_pack(bits,s,1);
         frame_bits_pack(bits,best_id,3);
      }
      ind[i]=best_index;
      gains[i]=best_gain;

      for (j=0;j<nsf;j++)
         e[j]=0;
      for (j=0;j<subvect_size;j++)
         e[subvect_size*i+j]=best_gain*shape_cb[best_index*subvect_size+j];
      residue_zero(e, awk1, r, nsf, p);
      syn_filt_zero(r, ak, r, nsf, p);
      syn_filt_zero(r, awk2, r, nsf,p);
      for (j=0;j<nsf;j++)
         tresp[i*nsf+j]=r[j];
      for (j=0;j<nsf;j++)
         t[j]-=r[j];
   }
   
   /* Put everything back together */
   for (i=0;i<nb_subvect;i++)
      for (j=0;j<subvect_size;j++)
         e[subvect_size*i+j]=gains[i]*shape_cb[ind[i]*subvect_size+j];

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




void split_cb_unquant(
float *exc,
void *par,                      /* non-overlapping codebook */
int   nsf,                      /* number of samples in subframe */
FrameBits *bits,
float *stack
)
{
   int i,j;
   int *ind;
   float *gains;
   float *sign;
   float *shape_cb, exc_energy;
   int shape_cb_size, subvect_size, nb_subvect;
   split_cb_params *params;

   params = (split_cb_params *) par;
   subvect_size = params->subvect_size;
   nb_subvect = params->nb_subvect;
   shape_cb_size = 1<<params->shape_bits;
   shape_cb = params->shape_cb;
   
   ind = (int*)PUSH(stack, nb_subvect);
   gains = PUSH(stack, nb_subvect);
   sign = PUSH(stack, nb_subvect);

   /* Decode global (average) gain */
   {
      int id;
      id = frame_bits_unpack_unsigned(bits, 4);
      exc_energy=exp(.5*id+2);
   }

   /* Decode codewords and gains */
   for (i=0;i<nb_subvect;i++)
   {
      int gain_id;
      ind[i] = frame_bits_unpack_unsigned(bits, params->shape_bits);
      if (frame_bits_unpack_unsigned(bits, 1))
         sign[i]=-1;
      else
         sign[i]=1;
      
      gain_id = frame_bits_unpack_unsigned(bits, 3);
      gains[i]=scal_gains4[gain_id];
      gains[i] *= sign[i];
      gains[i] *= exc_energy;


   }

   /* Compute decoded excitation */
   for (i=0;i<nb_subvect;i++)
      for (j=0;j<subvect_size;j++)
         exc[subvect_size*i+j]+=gains[i]*shape_cb[ind[i]*subvect_size+j];

   POP(stack);
   POP(stack);
   POP(stack);
}
