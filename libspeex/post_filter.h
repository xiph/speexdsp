/* Copyright (C) 2002 Jean-Marc Valin 
   File: post_filter.h
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

#ifndef POST_FILTER_H
#define POST_FILTER_H

/** Post-filter parameters*/
typedef struct pf_params {
   float formant_enh;
   float pitch_enh;
} pf_params;

/** Perceptual post-filter for narrowband */
void nb_post_filter(
float *exc, 
float *new_exc, 
float *ak, 
int p, 
int nsf,
int pitch,
float *pitch_gain,
void *params,
float *stack);


#endif
