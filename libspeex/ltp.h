/* Copyright (C) 2002 Jean-Marc Valin 
   File: ltp.h
   Lont-Term Prediction functions
*/

int open_loop_ltp(float *x, int len, int start, int end, float *gain);

int three_tap_ltp(float *x, int len, int start, int end, float gain[3]);
