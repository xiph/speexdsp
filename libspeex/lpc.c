/*
  Copyright 1992, 1993, 1994 by Jutta Degener and Carsten Bormann,
  Technische Universitaet Berlin

  Any use of this software is permitted provided that this notice is not
  removed and that neither the authors nor the Technische Universitaet Berlin
  are deemed to have made any representations as to the suitability of this
  software for any purpose nor are held responsible for any defects of
  this software.  THERE IS ABSOLUTELY NO WARRANTY FOR THIS SOFTWARE.

  As a matter of courtesy, the authors request to be informed about uses
  this software has found, about bugs in this software, and about any
  improvements that may be of general interest.

  Berlin, 28.11.1994
  Jutta Degener
  Carsten Bormann
*/


/* LPC- and Reflection Coefficients
 *
 * The next two functions calculate linear prediction coefficients
 * and/or the related reflection coefficients from the first P_MAX+1
 * values of the autocorrelation function.
 */

/* Invented by N. Levinson in 1947, modified by J. Durbin in 1959.
 */

#include "lpc.h"

float                      /* returns minimum mean square error    */
wld(
    float       * lpc, /*      [0...p-1] LPC coefficients      */
    const float * ac,  /*  in: [0...p] autocorrelation values  */
    float       * ref, /* out: [0...p-1] reflection coef's     */
    int p
    )
{
   int i, j;  float r, error = ac[0];

   if (ac[0] == 0) {
      for (i = 0; i < p; i++) ref[i] = 0; return 0; }

   for (i = 0; i < p; i++) {

      /* Sum up this iteration's reflection coefficient.
       */
      r = -ac[i + 1];
      for (j = 0; j < i; j++) r -= lpc[j] * ac[i - j];
      ref[i] = r /= error;

      /*  Update LPC coefficients and total error.
       */
      lpc[i] = r;
      for (j = 0; j < i/2; j++) {
         float tmp  = lpc[j];
         lpc[j]     += r * lpc[i-1-j];
         lpc[i-1-j] += r * tmp;
      }
      if (i % 2) lpc[j] += lpc[j] * r;

      error *= 1.0 - r * r;
   }
   return error;
}


/* Compute the autocorrelation
 *                      ,--,
 *              ac(i) = >  x(n) * x(n-i)  for all n
 *                      `--'
 * for lags between 0 and lag-1, and x == 0 outside 0...n-1
 */
void autocorr(
              const float * x,   /*  in: [0...n-1] samples x   */
              float *ac,   /* out: [0...lag-1] ac values */
              int lag, int   n)
{
   float d; int i;
   while (lag--) {
      for (i = lag, d = 0; i < n; i++) d += x[i] * x[i-lag];
      ac[lag] = d;
   }
}
