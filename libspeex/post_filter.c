/* Copyright (C) 2002 Jean-Marc Valin 
   File: post_filter.c
   Post-filtering routines: processing to enhance perceptual quality

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

#include <math.h>
#include "filters.h"
#include "stack_alloc.h"

/* Perceptual post-filter for narrowband */
void nb_post_filter(
float *exc, 
float *new_exc, 
float *ak, 
int p, 
int nsf,
int pitch,
float pitch_gain,
float *stack)
{
   int i;
   float exc_energy=0, new_exc_energy=0;
   float *awk;
   float gain;

   awk = PUSH(stack, p);
   
   bw_lpc (0.1, ak, awk, p);

   for (i=0;i<nsf;i++)
      exc_energy+=exc[i]*exc[i];

   for (i=0;i<nsf;i++)
   {
      new_exc[i] = exc[i] + .4*pitch_gain*exc[i-pitch];
   }
   
   syn_filt(new_exc, awk, new_exc, nsf, p);

   for (i=0;i<nsf;i++)
      new_exc_energy+=new_exc[i]*new_exc[i];

   gain = sqrt(exc_energy)/sqrt(.1+new_exc_energy);
   for (i=0;i<nsf;i++)
      new_exc[i] *= gain;

   POP(stack);
}
