/* Copyright (C) 2002 Jean-Marc Valin 
   File: vbr.h

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


#ifndef VBR_H
#define VBR_H

#define VBR_MEMORY_SIZE 5

typedef struct VBRState {
   float energy_alpha;
   float average_energy;
   float last_energy;
   float last_log_energy[VBR_MEMORY_SIZE];
   float accum_sum;
   float last_pitch_coef;
   float soft_pitch;
   float last_quality;
   float noise_level;
   float noise_accum;
   float noise_accum_count;
   int   consec_noise;
} VBRState;

void vbr_init(VBRState *vbr);

float vbr_analysis(VBRState *vbr, float *sig, int len, int pitch, float pitch_coef);

void vbr_destroy(VBRState *vbr);

#endif
