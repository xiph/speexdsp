#ifndef LPC_H
#define LPC_H

void autocorr(
        const float * x,   /*  in: [0...n-1] samples x   */
        float *ac,   /* out: [0...lag-1] ac values */
        int lag, int   n);

float                      /* returns minimum mean square error    */
wld(
        float       * lpc, /*      [0...p-1] LPC coefficients      */
        const float * ac,  /*  in: [0...p] autocorrelation values  */
        float       * ref, /* out: [0...p-1] reflection coef's     */
        int p
        );


#endif
