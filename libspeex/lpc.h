/* Copyright (C) 2002 Jean-Marc Valin 
   File: lpc.h
   Functions for LPC (Linear Prediction Coefficients) analysis

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
