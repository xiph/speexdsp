/* Copyright (C) 2002 Jean-Marc Valin 
   File: ltp.h
   Lont-Term Prediction functions
*/


/** Computes the open-loop pitch prediction. Returns pitch period and pitch gain */
int open_loop_ltp(float *x, int len, int start, int end, float *gain);


/** Computes a 3-tap pitch predictor */
int three_tap_ltp(float *x, int len, int start, int end, float *gain);


/** In place 3-tap pitch predictor (FIR)*/
void predictor_three_tap(float *x, int len, int period, float *gain);


/** In place 3-tap inverse pitch predictor (IIR)*/
void inverse_three_tap(float *x, int len, int period, float *gain);
