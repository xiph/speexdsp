#ifndef SPEEX_H
#define SPEEX_H


typedef struct EncState {
   int    frameSize;
   int    subframeSize;
   int    lpcSize;
   int    bufSize;
   float *inBuf;
   float *window;
   int    windowSize;
   float *buf2;
   float *autocorr;
   float *lpc;
   float *lsf;
   float *rc;
} EncState;

typedef struct DecState {
} DecState;

void encoder_init(EncState *st);
void encoder_destroy(EncState *st);
void encode(EncState *st, float *in, int *outSize, void *bits);

void decoder_init(DecState *st);
void decoder_destroy(DecState *st);
void decode(DecState *st, float *bits, float *out);



#endif
