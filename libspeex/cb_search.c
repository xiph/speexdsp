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

#define min(a,b) ((a) < (b) ? (a) : (b))

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





float split_cb_search(
float target[],			/* target vector */
float ak[],			/* LPCs for this subframe */
float awk1[],			/* Weighted LPCs for this subframe */
float awk2[],			/* Weighted LPCs for this subframe */
float codebook[][8],		/* overlapping codebook */
int   entries,			/* number of entries to search */
float *gain,			/* gain of optimum entries */
int   *index,			/* index of optimum entries */
int   p,                        /* number of LPC coeffs */
int   nsf,                      /* number of samples in subframe */
float *exc
)
{
   int i,j;
   float resp[64][8], E[64];
   float t[40], r[40], e[40];
   float gains[5];
   for (i=0;i<40;i++)
      t[i]=target[i];
   for (i=0;i<64;i++)
   {
      residue_zero(codebook[i], awk1, resp[i], 8, p);
      syn_filt_zero(resp[i], ak, resp[i], 8, p);
      syn_filt_zero(resp[i], awk2, resp[i], 8,p);
      E[i]=0;
      for(j=0;j<8;j++)
         E[i]+=resp[i][j]*resp[i][j];
   }
   for (i=0;i<5;i++)
   {
      int best_index;
      float corr, best_gain, score, best_score=-1;
      for (j=0;j<64;j++)
      {
         corr=xcorr(resp[j],t+8*i,8);
         score=corr*corr/(.001+E[j]);
         if (score>best_score)
         {
            best_index=j;
            best_score=score;
            best_gain=corr/(.001+E[j]);
         }
      }

      if (1) { /* Simulating scalar quantization of the gain*/
         float sign=1;
         printf("before: %f\n", best_gain);
         if (best_gain<0)
            sign=-1;
         best_gain = abs(best_gain)+.1;
         best_gain = log(best_gain);
         if (best_gain>8)
            best_gain=8;
         if (best_gain<0)
            best_gain=0;
         /*best_gain=.25*rint(4*best_gain);*/
         best_gain=.25*floor(4*best_gain+.5);
         best_gain=sign*exp(best_gain);
         printf("after: %f\n", best_gain);
      }
      gains[i]=best_gain;

      printf ("search: %d %f %f %f\n", best_index, best_gain, best_gain*best_gain*E[best_index], best_score);
      for (j=0;j<40;j++)
         e[j]=0;
      for (j=0;j<8;j++)
         e[8*i+j]=best_gain*codebook[best_index][j];
      residue_zero(e, awk1, r, 40, p);
      syn_filt_zero(r, ak, r, 40, p);
      syn_filt_zero(r, awk2, r, 40,p);
      for (j=0;j<40;j++)
         t[j]-=r[j];

      /*FIXME: Should move that out of the loop if we are to vector-quantize the gains*/
      for (j=0;j<40;j++)
         exc[j]+=e[j];
   }

   for (i=0;i<5;i++)
      printf ("%f ", gains[i]);
   printf ("cbgains: \n");
   /*TODO: Perform joint optimization of gains and quantization with prediction*/
   
   for (i=0;i<40;i++)
      target[i]=t[i];
}

