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
#include "roots.h"
#include <math.h>

#define min(a,b) ((a) < (b) ? (a) : (b))

#define MAX_ORD 20

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

void poly(float *re, float *im, float *p, int ord)
{
   int i,j;

   float p_re[MAX_ORD], p_im[MAX_ORD];
   for(i=0;i<ord+1;i++)
      p_re[i]=p_im[i]=0;
   p_re[0]=1;
   for (i=0;i<ord;i++)
   {
      for (j=i;j>=0;j--)
      {
         /* complex version of: p[j+1] -= p[j]*root[i] */
         p_re[j+1] -= p_re[j]*re[i] - p_im[j]*im[i];
         p_im[j+1] -= p_re[j]*im[i] + p_im[j]*re[i];
      }
   }
   for (i=0;i<ord+1;i++)
      p[i]=p_re[i];
}

/*LPC polynomial "flatifier"*/
void lpc_flat(float gamma, float *lpc_in, float *lpc_out, int order)
{
   int i;
   float re[10], im[10], conv[10];
   float alpha;
   alpha = 1/(4-4*gamma);
   poly_roots(lpc_in, re, im, conv, 10, 20, 7);
   for (i=0;i<order;i++)
   {
      float fact,tmp;
      float radius = sqrt(re[i]*re[i]+im[i]*im[i]);
      tmp=1-radius;
      if (tmp>2-2*gamma)
         fact = tmp;
      else
         fact = alpha*tmp*tmp-gamma+1;
      fact = (1-fact)/(radius+.001);
      re[i]*=fact;
      im[i]*=fact;
   }
   poly(re, im, lpc_out, order);
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

void syn_percep_zero(float *xx, float *ak, float *awk1, float *awk2, float *y, int N, int ord)
{
   int i,j;
   float long_filt[21];
   float iir[20];
   float fir[10];
   float x[40];
   int ord2=ord<<1;
   for (i=0;i<=ord;i++)
      long_filt[i]=ak[i];
   for (i=ord+1;i<=ord2;i++)
      long_filt[i]=0;
   residue_zero(long_filt,awk2,long_filt,1+(ord<<1),ord);

   for (i=0;i<N;i++)
      x[i]=xx[i];
   for (i=0;i<ord2;i++)
      iir[i]=long_filt[ord2-i];
   for (i=0;i<ord;i++)
      fir[i]=awk1[ord-i];

   for (i=0;i<ord2;i++)
   {
      y[i]=x[i];
      for (j=1;j<=min(ord2,i);j++)
         y[i] -= long_filt[j]*y[i-j];
      for (j=1;j<=min(ord,i);j++)
         y[i] += awk1[j]*x[i-j];
   }
#if 0
   for (i=ord2;i<N;i++)
   {
      float *yptr=y+i-ord2;
      float *xptr=x+i-ord;
      y[i]=x[i];
      for (j=0;j<ord2;j++)
         y[i]-= iir[j]*yptr[j];
      for (j=0;j<ord;j++)
         y[i]+= fir[j]*xptr[j];
   }
#else
   for (i=ord2;i<N;i++)
   {
      float *f, *ptr;
      float sum1=x[i], sum2=0;
      f=iir;
      ptr=y+i-ord2;
      for (j=0;j<ord2;j+=2)
      {
         sum1-= f[0]*ptr[0];
         sum2-= f[1]*ptr[1];
         f+=2;
         ptr+=2;
      }
      f=fir;
      ptr=x+i-ord;
      for (j=0;j<ord;j+=2)
      {
         sum1+= f[0]*ptr[0];
         sum2+= f[1]*ptr[1];
         f+=2;
         ptr+=2;
      }
      y[i]=sum1+sum2;
   }
#endif
}


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
