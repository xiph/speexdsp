#include "speex.h"
#include "lpc.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

void encoder_init(EncState *st)
{
   int i;
   st->frameSize = 128;
   st->windowSize = 256;
   st->subframeSize=64;
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
   st->lpc = malloc(st->lpcSize*sizeof(float));
   st->autocorr = malloc(st->lpcSize*sizeof(float));
   st->lsf = malloc(st->lpcSize*sizeof(float));
   st->rc = malloc(st->lpcSize*sizeof(float));
}

void encoder_destroy(EncState *st)
{
   free(st->inBuf);
   free(st->window);
   free(st->buf2);
   free(st->lpc);
   free(st->autocorr);
   free(st->lsf);
   free(st->rc);
}

void encode(EncState *st, float *in, int *outSize, void *bits)
{
   int i;
   float error;

   /* Copy new data in input buffer */
   memmove(st->inBuf, st->inBuf+st->bufSize-st->frameSize, (st->bufSize-st->frameSize)*sizeof(float));
   for (i=0;i<st->frameSize;i++)
      st->inBuf[st->bufSize-st->frameSize+i] = in[i];

   /* Window for analysis */
   for (i=0;i<st->windowSize;i++)
      st->buf2[i] = st->inBuf[i] * st->window[i];
   
   autocorr(st->buf2, st->autocorr, st->lpcSize-1, st->windowSize);
   st->autocorr[0] += 1;        /* prevents nan */
   st->autocorr[0] *= 1.0001;   /* 40 dB noise floor */
   /*Should do lag windowing here */
   error = wld(st->lpc, st->autocorr, st->rc, st->lpcSize-1);
   printf ("prediction error = %f, R[0] = %f, gain = %f\n", error, st->autocorr[0], 
           st->autocorr[0]/error);

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
