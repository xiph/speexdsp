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

extern float exc_gains_wb2_table[];
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

      if (gain_cb) /*If no gain codebook, do not quantize (for testing/debugging) */
      {
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
         if (nb_subvect<=5)
         {
         best_vq_index = vq_index(gains, gain_cb, nb_subvect, gain_cb_size);
         frame_bits_pack(bits,best_vq_index,params->gain_bits);
         printf ("best_gains_vq_index %d %f %d\n", best_vq_index, min_dist, max_index);
         for (i=0;i<nb_subvect;i++)
            gains[i]= sign[i]*gain_cb[best_vq_index*nb_subvect+i]/max_gain/(Ee[ind[i]]+.001);
         } else
         {
            float tmp[5];
            int best_vq_index2;
         best_vq_index = vq_index(gains, gain_cb, nb_subvect/2, gain_cb_size);
         for (i=0;i<5;i++)
            tmp[i]=gains[i]-gain_cb[best_vq_index*nb_subvect/2+i];
         best_vq_index2 = vq_index(tmp, exc_gains_wb2_table, nb_subvect/2, 256);

         frame_bits_pack(bits,best_vq_index,params->gain_bits);
         printf ("best_gains_vq_index %d %f %d\n", best_vq_index, min_dist, max_index);
         for (i=0;i<nb_subvect/2;i++)
            gains[i]= sign[i]*(gain_cb[best_vq_index*nb_subvect/2+i]+exc_gains_wb2_table[best_vq_index2*nb_subvect/2+i])/max_gain/(Ee[ind[i]]+.001);


         best_vq_index = vq_index(gains+5, gain_cb, nb_subvect/2, gain_cb_size);
         frame_bits_pack(bits,best_vq_index,params->gain_bits);
         for (i=0;i<5;i++)
            tmp[i]=gains[i+5]-gain_cb[best_vq_index*nb_subvect/2+i];
         best_vq_index2 = vq_index(tmp, exc_gains_wb2_table, nb_subvect/2, 256);

         printf ("best_gains_vq_index %d %f %d\n", best_vq_index, min_dist, max_index);
         for (i=0;i<nb_subvect/2;i++)
            gains[i+5]= sign[i+5]*(gain_cb[best_vq_index*nb_subvect/2+i]+exc_gains_wb2_table[best_vq_index2*nb_subvect/2+i])/max_gain/(Ee[ind[i+5]]+.001);
         }

    

         POP(stack);
      } else {
         printf ("exc: ");
         for (i=0;i<nb_subvect;i++)
            printf ("%f ", gains[i]);
         printf ("\n");
         for (i=0;i<nb_subvect;i++)
            gains[i]= gains[i]/(Ee[ind[i]]+.001);
      }

      for (i=0;i<nb_subvect;i++)
         for (j=0;j<subvect_size;j++)
            exc[subvect_size*i+j]+=gains[i]*shape_cb[ind[i]*subvect_size+j];

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




void split_cb_search_wb(
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
   float *t, *r, *e, *tresp;
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
   tresp = PUSH(stack, shape_cb_size*nsf);
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
#if 0
   {
      int half,k;
   for (half=0;half<2;half++)
   {
      int nb_half=nb_subvect/2;
      int max_subvect=0;
      float max_energy=0;
      syn_filt_zero(t+half*nsf/2, awk1, r, nsf/2, p);
      residue_zero(r, ak, r, nsf/2, p);
      residue_zero(r, awk2, r, nsf/2,p);
      for (i=0;i<nb_half;i++)
      {
         float energy=0;
         for (k=0;k<subvect_size;k++)
            energy+=r[subvect_size*i+k]*r[subvect_size*i+k];
         if (energy>max_energy)
         {
            max_subvect=i;
            max_energy=energy;
         }
      }
      printf ("max_energy: %d %f\n", max_subvect, max_energy);
      
      for (i=0;i<nb_half;i++)
      {
         int nb_times=1;
         if (i==max_subvect)
            nb_times++;

         for (k=0;k<nb_times;k++)
         {
         int best_index=0;
         float g, corr, best_gain=0, score, best_score=-1;
         for (j=0;j<shape_cb_size;j++)
         {
            corr=xcorr(resp+j*subvect_size,t+subvect_size*(i+half*nb_half),subvect_size);
            score=corr*corr/(.001+E[j]);
            g = corr/(.001+E[j]);
            if (score>best_score)
            {
               best_index=j;
               best_score=score;
               best_gain=corr/(.001+E[j]);
            }
         }
         for (j=0;j<nsf;j++)
            e[j]=0;
         for (j=0;j<subvect_size;j++)
            e[subvect_size*(i+half*nb_half)+j]=best_gain*shape_cb[best_index*subvect_size+j];
         residue_zero(e, awk1, r, nsf, p);
         syn_filt_zero(r, ak, r, nsf, p);
         syn_filt_zero(r, awk2, r, nsf,p);
         for (j=0;j<nsf;j++)
            t[j]-=r[j];
         for (j=0;j<nsf;j++)
            exc[j]+=e[j];
         }
      }
   }
   }
#else
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
      /*for (j=0;j<nsf;j++)
        exc[j]+=e[j];*/
   }
   {
      float A[10][10];
      float b[10];
      float c[10];
      for (i=0;i<10;i++)
         for (j=0;j<10;j++)
            A[i][j]=xcorr(tresp+i*nsf, tresp+j*nsf, nsf);
      for (i=0;i<10;i++)
         b[i]=xcorr(target,tresp+i*nsf,nsf);
      for (i=0;i<10;i++)
         A[i][i]+=.01;


      solve(&A[0][0],b,c, 10);
      for (i=0;i<10;i++)
        gains[i]*=c[i];
      
      for (i=0;i<10;i++)
         gains[i]*=Ee[ind[i]];



   {
      int best_vq_index=0, max_index;
      float max_gain=0, log_max, min_dist=0, *sign;

      if (gain_cb) /*If no gain codebook, do not quantize (for testing/debugging) */
      {
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
         if (nb_subvect<=5)
         {
         best_vq_index = vq_index(gains, gain_cb, nb_subvect, gain_cb_size);
         frame_bits_pack(bits,best_vq_index,params->gain_bits);
         printf ("best_gains_vq_index %d %f %d\n", best_vq_index, min_dist, max_index);
         for (i=0;i<nb_subvect;i++)
            gains[i]= sign[i]*gain_cb[best_vq_index*nb_subvect+i]/max_gain/(Ee[ind[i]]+.001);
         } else
         {
            float tmp[5];
            int best_vq_index2;
         best_vq_index = vq_index(gains, gain_cb, nb_subvect/2, gain_cb_size);
         for (i=0;i<5;i++)
            tmp[i]=gains[i]-gain_cb[best_vq_index*nb_subvect/2+i];
         best_vq_index2 = vq_index(tmp, exc_gains_wb2_table, nb_subvect/2, 256);

         frame_bits_pack(bits,best_vq_index,params->gain_bits);
         printf ("best_gains_vq_index %d %f %d\n", best_vq_index, min_dist, max_index);
         for (i=0;i<nb_subvect/2;i++)
            gains[i]= sign[i]*(gain_cb[best_vq_index*nb_subvect/2+i]+exc_gains_wb2_table[best_vq_index2*nb_subvect/2+i])/max_gain/(Ee[ind[i]]+.001);


         best_vq_index = vq_index(gains+5, gain_cb, nb_subvect/2, gain_cb_size);
         frame_bits_pack(bits,best_vq_index,params->gain_bits);
         for (i=0;i<5;i++)
            tmp[i]=gains[i+5]-gain_cb[best_vq_index*nb_subvect/2+i];
         best_vq_index2 = vq_index(tmp, exc_gains_wb2_table, nb_subvect/2, 256);

         printf ("best_gains_vq_index %d %f %d\n", best_vq_index, min_dist, max_index);
         for (i=0;i<nb_subvect/2;i++)
            gains[i+5]= sign[i+5]*(gain_cb[best_vq_index*nb_subvect/2+i]+exc_gains_wb2_table[best_vq_index2*nb_subvect/2+i])/max_gain/(Ee[ind[i+5]]+.001);
         }
         
         
      } else {
         
         for (i=0;i<10;i++)
            gains[i]/=Ee[ind[i]]+.001;
         
      }
   }
      for (i=0;i<10;i++)
         for (j=0;j<subvect_size;j++)
            e[subvect_size*i+j]=gains[i]*shape_cb[ind[i]*subvect_size+j];

      for (j=0;j<nsf;j++)
         exc[j]+=e[j];
      residue_zero(e, awk1, r, nsf, p);
      syn_filt_zero(r, ak, r, nsf, p);
      syn_filt_zero(r, awk2, r, nsf,p);
      for (j=0;j<nsf;j++)
         target[j]-=r[j];

   }
#endif

   /*TODO: Perform joint optimization of gains*/
   
   /*for (i=0;i<nsf;i++)
      target[i]=t[i];
   */
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
   int max_gain_ind, vq_gain_ind;
   float max_gain, *Ee;
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
   
   ind = (int*)PUSH(stack, nb_subvect);
   gains = PUSH(stack, nb_subvect);
   sign = PUSH(stack, nb_subvect);
   Ee=PUSH(stack, nb_subvect);

   for (i=0;i<nb_subvect;i++)
   {
      ind[i] = frame_bits_unpack_unsigned(bits, params->shape_bits);
      if (frame_bits_unpack_unsigned(bits, 1))
         sign[i]=-1;
      else
         sign[i]=1;
      Ee[i]=.001;
      for (j=0;j<subvect_size;j++)
         Ee[i]+=shape_cb[ind[i]*subvect_size+j]*shape_cb[ind[i]*subvect_size+j];
   }
   max_gain_ind = frame_bits_unpack_unsigned(bits, 3);
   vq_gain_ind = frame_bits_unpack_unsigned(bits, params->gain_bits);
   printf ("unquant gains ind: %d %d\n", max_gain_ind, vq_gain_ind);

   max_gain=exp(max_gain_ind+3.0);
   for (i=0;i<nb_subvect;i++)
      gains[i] = sign[i]*gain_cb[vq_gain_ind*nb_subvect+i]*max_gain/Ee[i];
   
   printf ("unquant gains: ");
   for (i=0;i<nb_subvect;i++)
      printf ("%f ", gains[i]);
   printf ("\n");

   for (i=0;i<nb_subvect;i++)
      for (j=0;j<subvect_size;j++)
         exc[subvect_size*i+j]+=gains[i]*shape_cb[ind[i]*subvect_size+j];
   
   POP(stack);
   POP(stack);
   POP(stack);
   POP(stack);
}
