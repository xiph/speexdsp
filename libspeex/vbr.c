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
#include <math.h>
#include <stdio.h>

#define sqr(x) ((x)*(x))

#define MIN_ENERGY 1000
#define NOISE_POW .3

void vbr_init(VBRState *vbr)
{
   int i;

   vbr->average_energy=0;
   vbr->last_energy=1;
   vbr->accum_sum=0;
   vbr->energy_alpha=.1;
   vbr->soft_pitch=0;
   vbr->last_pitch_coef=0;
   vbr->last_quality=0;

   vbr->noise_accum = .05*pow(MIN_ENERGY, NOISE_POW);
   vbr->noise_accum_count=.05;
   vbr->noise_level=vbr->noise_accum/vbr->noise_accum_count;
   vbr->consec_noise=0;


   for (i=0;i<VBR_MEMORY_SIZE;i++)
      vbr->last_log_energy[i] = log(MIN_ENERGY);
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
   int va;
   float log_energy;
   float non_st=0;
   float voicing;
   float pow_ener;

   for (i=0;i<len>>1;i++)
      ener1 += sig[i]*sig[i];

   for (i=len>>1;i<len;i++)
      ener2 += sig[i]*sig[i];
   ener=ener1+ener2;

   log_energy = log(ener+MIN_ENERGY);
   for (i=0;i<VBR_MEMORY_SIZE;i++)
      non_st += sqr(log_energy-vbr->last_log_energy[i]);
   non_st =  non_st/(30*VBR_MEMORY_SIZE);
   if (non_st>1)
      non_st=1;

   voicing = 3*(pitch_coef-.4)*fabs(pitch_coef-.4);
   vbr->average_energy = (1-vbr->energy_alpha)*vbr->average_energy + vbr->energy_alpha*ener;
   vbr->noise_level=vbr->noise_accum/vbr->noise_accum_count;
   pow_ener = pow(ener,NOISE_POW);
   if ((voicing<.3 && non_st < .2 && pow_ener < 1.2*vbr->noise_level)
       || (voicing<.2 && non_st < .1))
   {
      float tmp;
      va = 0;
      vbr->consec_noise++;
      if (pow_ener > 3*vbr->noise_level)
         tmp = 3*vbr->noise_level;
      else 
         tmp = pow_ener;
      if (vbr->consec_noise>=4)
      {
         vbr->noise_accum = .95*vbr->noise_accum + .05*tmp;
         vbr->noise_accum_count = .95*vbr->noise_accum_count + .05;
      }
   } else {
      va = 1;
      vbr->consec_noise=0;
   }

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
         qual += .7;
      if (ener > vbr->last_energy*1.8)
         qual += .7;
      if (ener > 3*vbr->average_energy)
         qual += .7;
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

   if (vbr->consec_noise>=1)
      qual-=1.2;
   if (vbr->consec_noise>=4)
      qual-=1.2;
   if (vbr->consec_noise>=8)
      qual-=1.2;

   vbr->last_pitch_coef = pitch_coef;
   vbr->last_quality = qual;

   for (i=VBR_MEMORY_SIZE-1;i>0;i--)
      vbr->last_log_energy[i] = vbr->last_log_energy[i-1];
   vbr->last_log_energy[0] = log_energy;

   /*printf ("VBR: %f %f %f %d %f\n", (float)(log_energy-log(vbr->average_energy+MIN_ENERGY)), non_st, voicing, va, vbr->noise_level);*/

   return qual;
}

void vbr_destroy(VBRState *vbr)
{
}
