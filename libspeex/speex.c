#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "speex.h"
#include "lpc.h"
#include "lsp.h"

#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif

#define sqr(x) ((x)*(x))

void encoder_init(EncState *st)
{
   int i;
   st->frameSize = 128;
   st->windowSize = 256;
   st->nbSubframes=4;
   st->subframeSize=32;
   st->lpcSize = 10;
   st->bufSize = 256;
   st->inBuf = malloc(st->bufSize*sizeof(float));
   for (i=0;i<st->bufSize;i++)
      st->inBuf[i]=0;
   st->window = malloc(st->windowSize*sizeof(float));
   /* Hanning window */
   for (i=0;i<st->windowSize;i++)
      st->window[i]=.5*(1-cos(2*M_PI*i/st->windowSize));
   st->buf2 = malloc(st->windowSize*sizeof(float));
   st->lpc = malloc((st->lpcSize+1)*sizeof(float));
   st->autocorr = malloc((st->lpcSize+1)*sizeof(float));
   /* Create the window for autocorrelation (lag-windowing) */
   st->lagWindow = malloc((st->lpcSize+1)*sizeof(float));
   for (i=0;i<st->lpcSize+1;i++)
      st->lagWindow[i]=exp(-.5*sqr(2*M_PI*.01*i));
   st->lsp = malloc(st->lpcSize*sizeof(float));
   st->rc = malloc(st->lpcSize*sizeof(float));
}

void encoder_destroy(EncState *st)
{
   free(st->inBuf);
   free(st->window);
   free(st->buf2);
   free(st->lpc);
   free(st->autocorr);
   free(st->lagWindow);
   free(st->lsp);
   free(st->rc);
}

void encode(EncState *st, float *in, int *outSize, void *bits)
{
   int i, roots;
   float error;

   /* Copy new data in input buffer */
   memmove(st->inBuf, st->inBuf+st->bufSize-st->frameSize, (st->bufSize-st->frameSize)*sizeof(float));
   for (i=0;i<st->frameSize;i++)
      st->inBuf[st->bufSize-st->frameSize+i] = in[i];

   /* Window for analysis */
   for (i=0;i<st->windowSize;i++)
      st->buf2[i] = st->inBuf[i] * st->window[i];
   /* Compute auto-correlation */
   autocorr(st->buf2, st->autocorr, st->lpcSize+1, st->windowSize);
   st->autocorr[0] += 1;        /* prevents NANs */
   st->autocorr[0] *= 1.0001;   /* 40 dB noise floor */
   /* Perform lag windowing here, equivalent to filtering in the power-spectrum domain */
   for (i=0;i<st->lpcSize+1;i++)
      st->autocorr[i] *= st->lagWindow[i];
   /* Levinson-Durbin */
   error = wld(st->lpc+1, st->autocorr, st->rc, st->lpcSize);
   st->lpc[0]=1;
   for (i=0;i<st->lpcSize+1;i++)
      printf("%f ", st->lpc[i]);
   printf ("\nprediction error = %f, R[0] = %f, gain = %f\n", error, st->autocorr[0], 
           st->autocorr[0]/error);
   /* LPC to LSPs (x-domain) transform */
   roots=lpc_to_lsp (st->lpc, st->lpcSize, st->lsp, 6, 0.02);
   for (i=0;i<roots;i++)
      printf("%f ", st->lsp[i]);
   printf ("\nfound %d roots\n", roots);

   /* Quantize LSPs */

   /* Back to "quantized" LPC */
   lsp_to_lpc(st->lsp, st->lpc, roots);
   for (i=0;i<st->lpcSize+1;i++)
      printf("%f ", st->lpc[i]);
   printf ("\n\n");

   /* Compute "weighted" residue */

   /* Find pitch */
}


void decoder_init(DecState *st)
{
}

void decoder_destroy(DecState *st)
{
}

void decode(DecState *st, float *bits, float *out)
{
}
