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


extern float gain_cdbk_nb[];

/** Computes the open-loop pitch prediction. Returns pitch period and pitch gain */
int open_loop_ltp(float *x, int len, int start, int end, float *gain);


/** Computes a 3-tap pitch predictor */
int three_tap_ltp(float *x, int len, int start, int end, float *gain);

/** Finds the best 3-tap pitch predictor from a codebook*/
int ltp_closed_loop(float *x, int len, int start, int end, float *gain);

/** In place 3-tap pitch predictor (FIR)*/
void predictor_three_tap(float *x, int len, int period, float *gain);


/** In place 3-tap inverse pitch predictor (IIR)*/
void inverse_three_tap(float *x, int len, int period, float *gain);
