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
float awk[],			/* Weighted LPCs for this subframe */
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
  int i,j,k;			/* loop variables */

  /* Initialise */
  
  resp = (float*)malloc(sizeof(float)*(p+nsf));
  h = (float*)malloc(sizeof(float)*(p+nsf));
  impulse = (float*)malloc(sizeof(float)*(p+nsf));

  for(i=0; i<p+nsf; i++)
    impulse[i] = 0.0;
  for(i=0; i<p; i++) {
    h[i] = 0.0;
    resp[i] = 0.0;
  }  
   
  *gain = 0.0;
  *index = 0;
  bscore = 0.0;
  impulse[p] = 1.0;

  /*synthesis_filter(impulse,ak,nsf,p,&h[p]);*/
  for (i=0;i<nsf;i++)
  {
     h[p+i] = impulse[p+i];
     for (k=1;k<=p;k++)
        h[p+i] += awk[k]*impulse[p+i-k] - ak[k]*h[p+i-k];
  }

  for (i=0;i<nsf;i++)
  {
     resp[p+i] = codebook[entries-1];
     for (k=1;k<=min(p,i);k++)
        resp[p+i] += awk[k]*codebook[entries-1+i-k];
     for (k=1;k<p+1;k++)
        resp[p+i] -= ak[k]*resp[p+i-k];
  }

  /*synthesis_filter(&codebook[entries-1],ak,nsf,p,&resp[p]);  */  
    
  /* Search codebook backwards using end correction for synthesis */
  
  for(k=entries-1; k>=0; k--) {

     /*FIXME: This is horribly slow, but I've got problems with the faster version.
      This loop should not be there at all once the bug is fixed*/
     for (i=0;i<nsf;i++)
     {
        resp[p+i] = codebook[k];
        for (j=1;j<=min(p,i);j++)
           resp[p+i] += awk[j]*codebook[k+i-j];
        for (j=1;j<p+1;j++)
           resp[p+i] -= ak[j]*resp[p+i-j];
     }

    d = 0.0; e = 0.0;
    for(i=0; i<nsf; i++) {
      d += target[i]*resp[p+i];
      e += resp[p+i]*resp[p+i];
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
        resp[i+p] = resp[i+p-1] + codebook[k-1]*h[p+i];
      resp[p] = codebook[k-1]*h[p];
    }
  }

  free(resp);
  free(h);
  free(impulse);

}

