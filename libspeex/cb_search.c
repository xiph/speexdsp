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

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))


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
int   nsf,                      /* number of samples in subframe */
float *stack
)
{
  float *resp;		        /* zero state response to current entry */
  float *h;		        /* impulse response of synthesis filter */
  float *impulse;		/* excitation vector containing one impulse */
  float d,e,g,score;		/* codebook searching variables */
  float bscore;			/* score of "best" vector so far */
  int i,k;			/* loop variables */

  /* Initialise */
  
  /*resp = (float*)malloc(sizeof(float)*nsf);
  h = (float*)malloc(sizeof(float)*nsf);
  impulse = (float*)malloc(sizeof(float)*nsf);
  */
  resp=PUSH(stack, nsf);
  h=PUSH(stack, nsf);
  impulse=PUSH(stack, nsf);

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
    g = d/(e+.0001);
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

  /*free(resp);
  free(h);
  free(impulse);*/
  POP(stack);
  POP(stack);
  POP(stack);
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
SpeexBits *bits,
float *stack
)
{
   int i,j, id;
   float *resp, *E, q;
   float *t, *r, *e;
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
   speex_bits_pack(bits, id, 4);
   exc_energy=exp(.5*q+2);


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
      /* Compute energy of codeword response */
      E[i]=0;
      for(j=0;j<subvect_size;j++)
         E[i]+=res[j]*res[j];
      E[i]=1/(.001+E[i]);
   }

   for (i=0;i<nb_subvect;i++)
   {
      int best_index=0, k, m;
      float g, corr, best_gain=0, score, best_score=-1;
      /* Find best codeword for current sub-vector */
      for (j=0;j<shape_cb_size;j++)
      {
         corr=xcorr(resp+j*subvect_size,t+subvect_size*i,subvect_size);
         score=corr*corr*E[j];
         g = corr*E[j];
         if (score>best_score)
         {
            best_index=j;
            best_score=score;
            best_gain=g;
         }
      }
      speex_bits_pack(bits,best_index,params->shape_bits);
      
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
         speex_bits_pack(bits,s,1);
         speex_bits_pack(bits,best_id,3);
      }
      ind[i]=best_index;
      gains[i]=best_gain;
      /* Update target for next subvector */
      for (j=0;j<subvect_size;j++)
      {
         g=best_gain*shape_cb[best_index*subvect_size+j];
         for (k=subvect_size*i+j,m=0;k<nsf;k++,m++)
            t[k] -= g*r[m];
      }
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
}

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
float *stack
)
{
   int i,j;
   float *resp;
   float *t, *r, *e;
   int *ind;
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
   ind = (int*)PUSH(stack, nb_subvect);

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
   }

   for (i=0;i<nb_subvect;i++)
   {
      int best_index=0, k, m;
      float g, dist, best_dist=-1;
      float *a, *b;

      /* Find best codeword for current sub-vector */
      for (j=0;j<shape_cb_size;j++)
      {
         dist=0;
         a=resp+j*subvect_size;
         b=t+subvect_size*i;
         for (k=0;k<subvect_size;k++)
            dist += (a[k]-b[k])*(a[k]-b[k]);
         if (dist<best_dist || j==0)
         {
            best_dist=dist;
            best_index=j;
         }
      }
      /*printf ("best index: %d/%d\n", best_index, shape_cb_size);*/
      speex_bits_pack(bits,best_index,params->shape_bits);

      ind[i]=best_index;
      /* Update target for next subvector */
      for (j=0;j<subvect_size;j++)
      {
         g=shape_cb[best_index*subvect_size+j];
         for (k=subvect_size*i+j,m=0;k<nsf;k++,m++)
            t[k] -= g*r[m];
      }

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
}


void split_cb_search_nogain2(
float target[],			/* target vector */
float ak[],			/* LPCs for this subframe */
float awk1[],			/* Weighted LPCs for this subframe */
float awk2[],			/* Weighted LPCs for this subframe */
void *par,                      /* Codebook/search parameters*/
int   p,                        /* number of LPC coeffs */
int   nsf,                      /* number of samples in subframe */
float *exc,
SpeexBits *bits,
float *stack
)
{
   int i,j;
   float *resp;
   float *t, *r, *e, *E;
   int *ind;
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
      float g, dist, best_dist[2]={-1,-1};
      float *a, *x;
      float energy=0;
      x=t+subvect_size*i;

      for (k=0;k<subvect_size;k++)
         energy+=x[k]*x[k];
      /* Find best codeword for current sub-vector */
      for (j=0;j<shape_cb_size;j++)
      {
         dist=0;
         a=resp+j*subvect_size;
         dist=energy+E[j];
         for (k=0;k<subvect_size;k++)
            dist -= 2*a[k]*x[k];
         if (dist<best_dist[0] || best_dist[0]<0)
         {
            best_dist[1]=best_dist[0];
            best_index[1]=best_index[0];
            best_dist[0]=dist;
            best_index[0]=j;
         } else if (dist<best_dist[1] || best_dist[1]<0)
         {
            best_dist[1]=dist;
            best_index[1]=j;
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
            for (j=0;j<nsf;j++)
               tt[j]=t[j];
            for (j=0;j<subvect_size;j++)
            {
               g=shape_cb[best_index[nbest]*subvect_size+j];
               for (k=subvect_size*i+j,m=0;k<nsf;k++,m++)
                  tt[k] -= g*r[m];
            }
            
            {
               int i2=vq_index(&tt[subvect_size*(i+1)], resp, subvect_size, shape_cb_size);
               for (j=0;j<subvect_size;j++)
               {
                  g=shape_cb[i2*subvect_size+j];
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
            best_index[0]=best_index[1];
            best_score[0]=best_score[1];
         }
         POP(stack);

      }

      ind[i]=best_index[0];

      /*printf ("best index: %d/%d\n", best_index, shape_cb_size);*/
      speex_bits_pack(bits,ind[i],params->shape_bits);

      /* Update target for next subvector */
      for (j=0;j<subvect_size;j++)
      {
         g=shape_cb[ind[i]*subvect_size+j];
         for (k=subvect_size*i+j,m=0;k<nsf;k++,m++)
            t[k] -= g*r[m];
      }
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
float *stack
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
}


void split_cb_search2(
float target[],			/* target vector */
float ak[],			/* LPCs for this subframe */
float awk1[],			/* Weighted LPCs for this subframe */
float awk2[],			/* Weighted LPCs for this subframe */
void *par,                      /* Codebook/search parameters*/
int   p,                        /* number of LPC coeffs */
int   nsf,                      /* number of samples in subframe */
float *exc,
SpeexBits *bits,
float *stack
)
{
   int i,j, id;
   float *resp, *E, q;
   float *t, *r, *e;
   float *gains;
   int *ind, *gain_ind;
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
   E = PUSH(stack, shape_cb_size);
   t = PUSH(stack, nsf);
   r = PUSH(stack, nsf);
   e = PUSH(stack, nsf);
   gains = PUSH(stack, nb_subvect);
   ind = (int*)PUSH(stack, nb_subvect);
   gain_ind = (int*)PUSH(stack, nb_subvect);

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
   speex_bits_pack(bits, id, 4);
   exc_energy=exp(.5*q+2);


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
      /* Compute energy of codeword response */
      E[i]=0;
      for(j=0;j<subvect_size;j++)
         E[i]+=res[j]*res[j];
      E[i]=1/(.001+E[i]);
   }

   for (i=0;i<nb_subvect;i++)
   {
      int best_index[2]={0,0}, k, m, best_gain_ind[2]={0,0};
      float g, corr, best_gain[2]={0,0}, score, best_score[2]={-1,-1};
      /* Find best codeword for current sub-vector */
      for (j=0;j<shape_cb_size;j++)
      {
         corr=xcorr(resp+j*subvect_size,t+subvect_size*i,subvect_size);
         score=corr*corr*E[j];
         g = corr*E[j];
         if (score>best_score[0])
         {
            best_index[1]=best_index[0];
            best_score[1]=best_score[0];
            best_gain[1]=best_gain[0];

            best_index[0]=j;
            best_score[0]=score;
            best_gain[0]=g;
         } else if (score>best_score[1]) {
            best_index[1]=j;
            best_score[1]=score;
            best_gain[1]=g;            
         }
      }
      
      /* Quantize gain */
      for (k=0;k<2;k++) {
         int s=0, best_id;
         best_gain[k] /= .01+exc_energy;
         if (best_gain[k]<0)
         {
            best_gain[k]=-best_gain[k];
            s=1;
         }

         /* Find gain index (it's a scalar but we use the VQ code anyway)*/
         best_id = vq_index(&best_gain[k], scal_gains4, 1, 8);

         best_gain_ind[k]=best_id;
         best_gain[k]=scal_gains4[best_id];
         /*printf ("gain_quant: %f %d %f\n", best_gain, best_id, scal_gains4[best_id]);*/
         if (s)
            best_gain[k]=-best_gain[k];
         best_gain[k] *= exc_energy;
      }



      if (i<nb_subvect-1) {
         int best_index2=0;
         float best_score2=-1, best_gain2=0;
         int nbest;
         float err[2]={0,0};
         float *tt=PUSH(stack,nsf);
         for (nbest=0;nbest<2;nbest++)
         {
            for (j=0;j<nsf;j++)
               tt[j]=t[j];
            for (j=0;j<subvect_size;j++)
            {
               g=best_gain[nbest]*shape_cb[best_index[nbest]*subvect_size+j];
               for (k=subvect_size*i+j,m=0;k<nsf;k++,m++)
                  tt[k] -= g*r[m];
            }
            

            for (j=0;j<shape_cb_size;j++)
            {
               corr=xcorr(resp+j*subvect_size,tt+subvect_size*(i+1),subvect_size);
               score=corr*corr*E[j];
               g = corr*E[j];
               if (score>best_score2)
               {
                  best_index2=j;
                  best_score2=score;
                  best_gain2=g;
               }
            }

            {
               int s=0, best_id;
               best_gain2 /= .01+exc_energy;
               if (best_gain2<0)
               {
                  best_gain2=-best_gain2;
                  s=1;
               }
               best_id = vq_index(&best_gain2, scal_gains4, 1, 8);
               best_gain2=scal_gains4[best_id];
               if (s)
                  best_gain2=-best_gain2;
               best_gain2 *= exc_energy;
            }

            for (j=0;j<subvect_size;j++)
            {
               g=best_gain2*shape_cb[best_index2*subvect_size+j];
               for (k=subvect_size*(i+1)+j,m=0;k<nsf;k++,m++)
                  tt[k] -= g*r[m];
            }
            for (j=subvect_size*i;j<subvect_size*(i+2);j++)
               err[nbest]-=tt[j]*tt[j];
            
            best_score[nbest]=err[nbest];
         }

         if (best_score[1]>best_score[0])
         {
            best_index[0]=best_index[1];
            best_score[0]=best_score[1];
            best_gain[0]=best_gain[1];
            best_gain_ind[0]=best_gain_ind[1];
         }
         POP(stack);
      }


      

      ind[i]=best_index[0];
      gain_ind[i]=best_gain_ind[0];
      gains[i]=best_gain[0];
      /* Update target for next subvector */
      for (j=0;j<subvect_size;j++)
      {
         g=best_gain[0]*shape_cb[best_index[0]*subvect_size+j];
         for (k=subvect_size*i+j,m=0;k<nsf;k++,m++)
            t[k] -= g*r[m];
      }
   }
   for (i=0;i<nb_subvect;i++)
   {
      speex_bits_pack(bits, ind[i], params->shape_bits);
      if (gains[i]<0)
         speex_bits_pack(bits, 1, 1);
      else
         speex_bits_pack(bits, 0, 1);
      speex_bits_pack(bits, gain_ind[i], 3);
      /*printf ("encode split: %d %d %f\n", i, ind[i], gains[i]);*/

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
SpeexBits *bits,
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
      id = speex_bits_unpack_unsigned(bits, 4);
      exc_energy=exp(.5*id+2);
   }

   /* Decode codewords and gains */
   for (i=0;i<nb_subvect;i++)
   {
      int gain_id;
      ind[i] = speex_bits_unpack_unsigned(bits, params->shape_bits);
      if (speex_bits_unpack_unsigned(bits, 1))
         sign[i]=-1;
      else
         sign[i]=1;
      
      gain_id = speex_bits_unpack_unsigned(bits, 3);
      gains[i]=scal_gains4[gain_id];
      gains[i] *= sign[i];
      gains[i] *= exc_energy;

      /*printf ("decode split: %d %d %f\n", i, ind[i], gains[i]);*/
   }

   /* Compute decoded excitation */
   for (i=0;i<nb_subvect;i++)
      for (j=0;j<subvect_size;j++)
         exc[subvect_size*i+j]+=gains[i]*shape_cb[ind[i]*subvect_size+j];

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
