/* Copyright (C) 2002 Jean-Marc Valin 
   File: vbr.c

   VBR-related routines

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

#include "vbr.h"

void vbr_init(VBRState *vbr)
{
   vbr->average_energy=0;
   vbr->last_energy=1;
   vbr->accum_sum=0;
   vbr->energy_alpha=.1;
   vbr->soft_pitch=0;
   vbr->last_pitch_coef=0;
   vbr->last_quality=0;
}


/*
  This function should analyse the signal and decide how critical the
  coding error will be perceptually. The following factors should be
  taken into account:

  -Attacks (positive energy derivative) should be coded with more bits

  -Stationary voiced segments should receive more bits

  -Segments with (very) low absolute energy should receive less bits (maybe
  only shaped noise?)

  -DTX for near-zero energy?

  -Stationary fricative segments should have less bits

  -Temporal masking: when energy slope is decreasing, decrease the bit-rate

  -Decrease bit-rate for males (low pitch)?

  -(wideband only) less bits in the high-band when signal is very 
  non-stationary (harder to notice high-frequency noise)???

*/
float vbr_analysis(VBRState *vbr, float *sig, int len, int pitch, float pitch_coef)
{
   int i;
   float ener=0, ener1=0, ener2=0;
   float qual=0;

   for (i=0;i<len>>1;i++)
      ener1 += sig[i]*sig[i];

   for (i=len>>1;i<len;i++)
      ener2 += sig[i]*sig[i];
   ener=ener1+ener2;

   vbr->average_energy = (1-vbr->energy_alpha)*vbr->average_energy + vbr->energy_alpha*ener;
   
   /* Checking for "pseudo temporal masking" */
   if (ener < .1*vbr->average_energy)
      qual -= .7;
   if (ener < .01*vbr->average_energy)
      qual -= .7;
   if (ener < .001*vbr->average_energy)
      qual -= .7;
   /* Checking for very low absolute energy */
   if (ener < 30000)
   {
      qual -= .7;
      if (ener < 10000)
         qual-=.7;
      if (ener < 3000)
         qual-=.7;
   } else {
      /* Checking for energy increases */
      if (ener > vbr->last_energy*4.0)
         qual += 1;
      if (ener > vbr->last_energy*1.8)
         qual += 1;
      if (ener > 3*vbr->average_energy)
         qual += 1;
      if (ener2 > 1.6*ener1)
         qual += .7;
      if (ener2 < .6*ener1)
         qual -= .5;

      if (ener < .3*vbr->last_energy)
         qual -= .6;
   }
   vbr->soft_pitch = .6*vbr->soft_pitch + .4*pitch_coef;
   qual += (pitch_coef-.4) + (vbr->soft_pitch-.4);

   if (qual < vbr->last_quality)
      qual = .5*qual + .5*vbr->last_quality;
   if (qual<-3)
      qual=-3;
   if (qual>3)
      qual=3;

   vbr->last_energy = ener;
   vbr->last_pitch_coef = pitch_coef;
   vbr->last_quality = qual;
   return qual;
}

void vbr_destroy(VBRState *vbr)
{
}
