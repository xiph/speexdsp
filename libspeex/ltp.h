/* Copyright (C) 2002 Jean-Marc Valin 
   File: ltp.h
   Lont-Term Prediction functions

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

#include "bits.h"

extern float gain_cdbk_nb[];

/* Finds open-loop pitch */
void open_loop_pitch(float *sw, int start, int end, int len, int *pitch, int *vuv);

void closed_loop_fractional_pitch(
float target[],                 /* Target vector */
float ak[],                     /* LPCs for this subframe */
float awk1[],                   /* Weighted LPCs #1 for this subframe */
float awk2[],                   /* Weighted LPCs #2 for this subframe */
float exc[],                    /* Overlapping codebook */
float *filt,                    /* Over-sampling filter */
int   filt_side,                /* Over-sampling factor */
int   fact,                     /* Over-sampling factor */
int   start,                    /* Smallest pitch value allowed */
int   end,                      /* Largest pitch value allowed */
float *gain,                    /* 3-tab gains of optimum entry */
int   *pitch,                   /* Index of optimum entry */
int   p,                        /* Number of LPC coeffs */
int   nsf,                      /* Number of samples in subframe */
float *stack
);

/** Finds the best quantized 3-tap pitch predictor by analysis by synthesis */
void pitch_search_3tap(
float target[],                 /* Target vector */
float ak[],                     /* LPCs for this subframe */
float awk1[],                   /* Weighted LPCs #1 for this subframe */
float awk2[],                   /* Weighted LPCs #2 for this subframe */
float exc[],                    /* Overlapping codebook */
int   start,                    /* Smallest pitch value allowed */
int   end,                      /* Largest pitch value allowed */
float *gain,                    /* 3-tab gains of optimum entry */
int   *pitch,                   /* Best pitch delay */
int   *gain_index,              /* Index of optimum gain */
int   p,                        /* Number of LPC coeffs */
int   nsf,                      /* Number of samples in subframe */
FrameBits *bits,
float *stack
);

/** Finds the best quantized 3-tap pitch predictor by analysis by synthesis */
float pitch_search_3tap_unquant(
float target[],                 /* Target vector */
float ak[],                     /* LPCs for this subframe */
float awk1[],                   /* Weighted LPCs #1 for this subframe */
float awk2[],                   /* Weighted LPCs #2 for this subframe */
float exc[],                    /* Overlapping codebook */
int   start,                    /* Smallest pitch value allowed */
int   end,                      /* Largest pitch value allowed */
float *gain,                    /* 3-tab gains of optimum entry */
int   *pitch,                   /* Index of optimum entry */
int   p,                        /* Number of LPC coeffs */
int   nsf,                      /* Number of samples in subframe */
float *stack
);

/*Unquantize adaptive codebook and update pitch contribution*/
void pitch_unquant_3tap(
float exc[],                    /* Excitation */
int   start,                    /* Smallest pitch value allowed */
int   end,                      /* Largest pitch value allowed */
int   nsf,                      /* Number of samples in subframe */
FrameBits *bits,
float *stack
);
