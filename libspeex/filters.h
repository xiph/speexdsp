/* Copyright (C) 2002 Jean-Marc Valin 
   File: filters.h
   Various analysis/synthesis filters

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

#ifndef FILTERS_H
#define FILTERS_H

void syn_filt(float *x, float *a, float *y, int N, int ord);
void syn_filt_zero(float *x, float *a, float *y, int N, int ord);
void syn_filt_mem(float *x, float *a, float *y, int N, int ord, float *mem);

void residue(float *x, float *a, float *y, int N, int ord);
void residue_zero(float *x, float *a, float *y, int N, int ord);
void residue_mem(float *x, float *a, float *y, int N, int ord, float *mem);

#endif
