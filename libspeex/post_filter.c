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
#include "post_filter.h"

/* Perceptual post-filter for narrowband */
void nb_post_filter(
float *exc,          /*decoded excitation*/
float *new_exc,      /*enhanced excitation*/
float *ak,           /*LPC filter coefs*/
int p,               /*LPC order*/
int nsf,             /*sub-frame size*/
int pitch,           /*pitch period*/
float *pitch_gain,   /*pitch gain (3-tap)*/
void *par,           /*post-filter parameters*/
float *mem,          /*filter memory #1*/
float *mem2,         /*filter memory #2*/
float *stack)
{
   int i;
   float exc_energy=0, new_exc_energy=0;
   float *awk1, *awk2;
   float gain;
   pf_params *params;
   float *tmp_exc;
   float formant_num, formant_den, voiced_fact;

   params = (pf_params*) par;

   awk1 = PUSH(stack, p);
   awk2 = PUSH(stack, p);
   tmp_exc = PUSH(stack, nsf);
  
   /*Compute excitation energy prior to enhancement*/
   for (i=0;i<nsf;i++)
      exc_energy+=exc[i]*exc[i];

   /*Apply pitch comb-filter (filter out noise between pitch harmonics)*/
   for (i=0;i<nsf;i++)
   {
      new_exc[i] = exc[i] + params->pitch_enh*(
                                               pitch_gain[0]*exc[i-pitch+1] +
                                               pitch_gain[1]*exc[i-pitch] +
                                               pitch_gain[2]*exc[i-pitch-1]
                                               );
   }
   
   /*Compute "voicing coefficient" (0 for unvoiced, 1 for voiced)*/
   voiced_fact = pitch_gain[0]+pitch_gain[1]+pitch_gain[2];
   if (voiced_fact<0)
      voiced_fact=0;
   if (voiced_fact>1)
      voiced_fact=1;

   /*Adapting post-filter to voicing*/
   formant_num = params->formant_enh_num;
   formant_den =   voiced_fact    * params->formant_enh_den 
                 + (1-voiced_fact)* params->formant_enh_num;

   /*Short-term post-filter using "flatified" versions of ak*/
   lpc_flat (formant_num, formant_den, ak, awk1, awk2, p);
   residue_mem(new_exc, awk1, tmp_exc, nsf, p, mem);
   syn_filt_mem(tmp_exc, awk2, new_exc, nsf, p, mem2);

   /*Gain after enhancement*/
   for (i=0;i<nsf;i++)
      new_exc_energy+=new_exc[i]*new_exc[i];

   /*Compute scaling factor and normalize energy*/
   gain = sqrt(exc_energy)/sqrt(.1+new_exc_energy);
   for (i=0;i<nsf;i++)
      new_exc[i] *= gain;
   
   POP(stack);
   POP(stack);
   POP(stack);
}
