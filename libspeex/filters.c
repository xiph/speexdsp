/* Copyright (C) 2002 Jean-Marc Valin 
   File: filters.c
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

#include "filters.h"
#include <stdio.h>

#define min(a,b) ((a) < (b) ? (a) : (b))

void print_vec(float *vec, int len, char *name)
{
   int i;
   printf ("%s ", name);
   for (i=0;i<len;i++)
      printf (" %f", vec[i]);
   printf ("\n");
}

void bw_lpc(float gamma, float *lpc_in, float *lpc_out, int order)
{
   int i;
   float tmp=1;
   for (i=0;i<order+1;i++)
   {
      lpc_out[i] = tmp * lpc_in[i];
      tmp *= gamma;
   }
}

void syn_filt(float *x, float *a, float *y, int N, int ord)
{
   int i,j;
   for (i=0;i<N;i++)
   {
      y[i]=x[i];
      for (j=1;j<=ord;j++)
         y[i] -= a[j]*y[i-j];
   }
}

void syn_filt_zero(float *x, float *a, float *y, int N, int ord)
{
   int i,j;
   for (i=0;i<N;i++)
   {
      y[i]=x[i];
      for (j=1;j<=min(ord,i);j++)
         y[i] -= a[j]*y[i-j];
   }
}

void syn_filt_mem(float *x, float *a, float *y, int N, int ord, float *mem)
{
   int i,j;
   for (i=0;i<N;i++)
   {
      y[i]=x[i];
      for (j=1;j<=min(ord,i);j++)
         y[i] -= a[j]*y[i-j];
      for (j=i+1;j<=ord;j++)
         y[i] -= a[j]*mem[j-i-1];
   }
   for (i=0;i<ord;i++)
      mem[i]=y[N-i-1];
}


void residue(float *x, float *a, float *y, int N, int ord)
{
   int i,j;
   for (i=N-1;i>=0;i--)
   {
      y[i]=x[i];
      for (j=1;j<=ord;j++)
         y[i] += a[j]*x[i-j];
   }
}

void residue_zero(float *x, float *a, float *y, int N, int ord)
{
   int i,j;
   for (i=N-1;i>=0;i--)
   {
      y[i]=x[i];
      for (j=1;j<=min(ord,i);j++)
         y[i] += a[j]*x[i-j];
   }
}

void residue_mem(float *x, float *a, float *y, int N, int ord, float *mem)
{
   int i,j;
   for (i=N-1;i>=0;i--)
   {
      y[i]=x[i];
      for (j=1;j<=min(ord,i);j++)
         y[i] += a[j]*x[i-j];
      for (j=i+1;j<=ord;j++)
         y[i] += a[j]*mem[j-i-1];
   }
   for (i=0;i<ord;i++)
      mem[i]=x[N-i-1];
}

float xcorr(float *x, float *y, int len)
{
   int i;
   float sum=0;
   for (i=0;i<len;i++)
      sum += x[i]*y[i];
   return sum;
}

/*
void fir_mem(float *x, float *a, float *y, int N, int M, float *mem)
{
   int i,j;
   for (i=0;i<N;i++)
   {
      y[i]=0;
      for (j=0;j<=min(M-1,i);j++)
         y[i] += a[j]*x[i-j];
      for (j=i+1;j<=M-1;j++)
         y[i] += a[j]*mem[j-i-1];
   }
   for (i=0;i<M-1;i++)
      mem[i]=x[N-i-1];
}
*/

#define MAX_FILTER 100
#define MAX_SIGNAL 1000
void fir_mem(float *xx, float *aa, float *y, int N, int M, float *mem)
{
   int i,j;
   float a[MAX_FILTER];
   float x[MAX_SIGNAL];
   for (i=0;i<M;i++)
      a[M-i-1]=aa[i];
   for (i=0;i<M-1;i++)
      x[i]=mem[M-i-2];
   for (i=0;i<N;i++)
      x[i+M-1]=xx[i];
   for (i=0;i<N;i++)
   {
      y[i]=0;
      for (j=0;j<M;j++)
         y[i]+=a[j]*x[i+j];
   }
   for (i=0;i<M-1;i++)
     mem[i]=xx[N-i-1];
}
