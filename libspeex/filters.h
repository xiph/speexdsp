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

void print_vec(float *vec, int len, char *name);

void filter_mem2(float *x, float *num, float *den, float *y, int N, int ord, float *mem);
void fir_mem2(float *x, float *num, float *y, int N, int ord, float *mem);
void iir_mem2(float *x, float *den, float *y, int N, int ord, float *mem);

/* Apply bandwidth expansion on LPC coef */
void bw_lpc(float gamma, float *lpc_in, float *lpc_out, int order);

void poly(float *re, float *im, float *p, int ord);

void enh_lpc(float *ak, int order, float *num, float *den, float k1, float k2, float *stack);

/*LPC polynomial "flatifier"*/
void lpc_flat(float g1, float g2, float *lpc_in, float *lpc_out1, float *lpc_out2, int order);


#if 0
/* Synthesis filter using the past of y[n] (negative indices) as memory */
void syn_filt(float *x, float *a, float *y, int N, int ord);
#endif

/* Synthesis filter using zero memory */
void syn_filt_zero(float *x, float *a, float *y, int N, int ord);

#if 0
/* Synthesis filter using memory */
void syn_filt_mem(float *x, float *a, float *y, int N, int ord, float *mem);

/* Analysis (FIR) filter using the past of x[n] (negative indices) as memory */
void residue(float *x, float *a, float *y, int N, int ord);
#endif

/* Analysis (FIR) filter using zero memory */
void residue_zero(float *x, float *a, float *y, int N, int ord);

#if 0
/* Analysis (FIR) filter using memory */
void residue_mem(float *x, float *a, float *y, int N, int ord, float *mem);
#endif



/* FIR filter */
void fir_mem(float *x, float *a, float *y, int N, int M, float *mem);

void syn_percep_zero(float *x, float *ak, float *awk1, float *awk2, float *y, int N, int ord);

void comb_filter(
float *exc,          /*decoded excitation*/
float *new_exc,      /*enhanced excitation*/
float *ak,           /*LPC filter coefs*/
int p,               /*LPC order*/
int nsf,             /*sub-frame size*/
int pitch,           /*pitch period*/
float *pitch_gain,   /*pitch gain (3-tap)*/
float  comb_gain     /*gain of comb filter*/
);

void pole_zero_mem(float *x, float *num, float *den, float *y, int N, int ord, float *mem, float *stack);

#endif
