/*-----------------------------------------------------------------------*\

    FILE........: GAINSHAPE.C
    TYPE........: C Module
    AUTHOR......: David Rowe
    COMPANY.....: Voicetronix
    DATE CREATED: 19/2/02

    General gain-shape codebbok search.

\*-----------------------------------------------------------------------*/

/* Modified by Jean-Marc Valin 2002 */


#include <stdlib.h>
#include <cb_search.h>
#include "filters.h"

#define min(a,b) ((a) < (b) ? (a) : (b))

/*---------------------------------------------------------------------------*\
                                                                             
 void overlap_cb_search()							      
									      
 Searches a gain/shape codebook consisting of overlapping entries for the    
 closest vector to the target.  Gives identical results to search() above   
 buts uses fast end correction algorithm for the synthesis of response       
 vectors.								      
                                                                             
\*---------------------------------------------------------------------------*/

void overlap_cb_search(
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
      e += resp[i]*resp[i]+1;
    }
    g = d/e;
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

}

