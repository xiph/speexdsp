#ifndef SPEEX_H
#define SPEEX_H

/**Structure representing the full state of the encoder*/
typedef struct EncState {
   int    frameSize;      /* Size of frames */
   int    subframeSize;   /* Size of sub-frames */
   int    nbSubframes;    /* Number of sub-frames */
   int    windowSize;     /* Analysis (LPC) window length */
   int    lpcSize;        /* LPC order */
   int    bufSize;        /* Buffer size */
   float *inBuf;          /* Input buffer */
   float *frame;          /* Start of encoded frame */
   float *window;         /* Temporary (Hanning) window */
   float *buf2;           /* 2nd temporary buffer */
   float *autocorr;       /* auto-correlation */
   float *lagWindow;      /* Window applied to auto-correlation */
   float *lpc;            /* LPCs for current frame */
   float *lsp;            /* LSPs for current frame */
   float *old_lsp;        /* LSPs for previous frame */
   float *interp_lsp;     /* Interpolated LSPs */
   float *interp_lpc;     /* Interpolated LPCs */
   float *rc;             /* Reflection coefficients */
   int    first;          /* Is this the first frame? */
} EncState;

typedef struct DecState {
} DecState;

/**Initializes encoder state*/
void encoder_init(EncState *st);
/**De-allocates encoder state resources*/
void encoder_destroy(EncState *st);
/**Encodes one frame*/
void encode(EncState *st, float *in, int *outSize, void *bits);

void decoder_init(DecState *st);
void decoder_destroy(DecState *st);
void decode(DecState *st, float *bits, float *out);



#endif
