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
void lpc_flat(float g1, float g2, float *lpc_in, float *lpc_out1, float *lpc_out2, int order)
{
   int i;
   float re[10], im[10], conv[10];
   float re1[10], im1[10];
   float alpha;
   poly_roots(lpc_in, re, im, conv, 10, 20, 4);
   for (i=0;i<order;i++)
   {
      re1[i]=re[i];
      im1[i]=im[i];
   }
   alpha = 1/(4-4*g1);
   for (i=0;i<order;i++)
   {
      float fact,tmp;
      float radius = sqrt(re[i]*re[i]+im[i]*im[i]);
      if (radius>1)
      {
         re[i]=im[i]=1;
         radius=1;
      }
      if (!(radius<1))
      {
         re[i]=im[i]=0;
         radius=0;
         /*fprintf(stderr, "radius = %f\n", radius);*/
      }
      tmp=1-radius;
      if (tmp>2-2*g1)
         fact = tmp;
      else
         fact = alpha*tmp*tmp-g1+1;
      fact = (1-fact)/(radius+.001);
      re[i]*=fact;
      im[i]*=fact;
   }
   poly(re, im, lpc_out1, order);
   for (i=0;i<order;i++)
   {
      re[i]=re1[i];
      im[i]=im1[i];
   }
   alpha = 1/(4-4*g2);
   for (i=0;i<order;i++)
   {
      float fact,tmp;
      float radius = sqrt(re[i]*re[i]+im[i]*im[i]);
      if (radius>1)
      {
         re[i]=im[i]=1;
         radius=1;
      }
      if (!(radius<1))
      {
         re[i]=im[i]=0;
         radius=0;
         /*fprintf(stderr, "radius = %f\n", radius);*/
      }
      tmp=1-radius;
      if (tmp>2-2*g2)
         fact = tmp;
      else
         fact = alpha*tmp*tmp-g2+1;
      fact = (1-fact)/(radius+.001);
      re[i]*=fact;
      im[i]*=fact;
   }
   poly(re, im, lpc_out2, order);

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

#ifndef _USE_SSE
float xcorr(float *x, float *y, int len)
{
   int i;
   float sum=0;
   for (i=0;i<len;i++)
      sum += x[i]*y[i];
   return sum;
}

#else
float xcorr(float *a, float *b, int len)
{
  float sum;
  __asm__ __volatile__ (
  "
  push %%eax
  push %%edi
  push %%ecx
  xorps %%xmm3, %%xmm3
  xorps %%xmm4, %%xmm4

  sub $8, %%ecx
  jb mul8_skip%=

mul8_loop%=:
  movups (%%eax), %%xmm0
  movups (%%edi), %%xmm1
  movups 16(%%eax), %%xmm5
  movups 16(%%edi), %%xmm6
  add $32, %%eax
  add $32, %%edi
  mulps %%xmm0, %%xmm1
  mulps %%xmm5, %%xmm6
  addps %%xmm1, %%xmm3
  addps %%xmm6, %%xmm4

  sub $8,  %%ecx

  jae mul8_loop%=

mul8_skip%=:

  addps %%xmm4, %%xmm3

  add $4, %%ecx
  jl mul4_skip%=

  movups (%%eax), %%xmm0
  movups (%%edi), %%xmm1
  add $16, %%eax
  add $16, %%edi
  mulps %%xmm0, %%xmm1
  addps %%xmm1, %%xmm3

  sub $4,  %%ecx

mul4_skip%=:


  add $4, %%ecx

  jmp cond1%=

mul1_loop%=:
  movss (%%eax), %%xmm0
  movss (%%edi), %%xmm1
  add $4, %%eax
  add $4, %%edi
  mulss %%xmm0, %%xmm1
  addss %%xmm1, %%xmm3

cond1%=:
  sub $1, %%ecx
  jae mul1_loop%=

  movhlps %%xmm3, %%xmm4
  addps %%xmm4, %%xmm3
  movaps %%xmm3, %%xmm4
  //FIXME: which one?
  shufps $0x55, %%xmm4, %%xmm4
  //shufps $33, %%xmm4, %%xmm4
  addss %%xmm4, %%xmm3
  movss %%xmm3, (%%edx)
  
  pop %%ecx
  pop %%edi
  pop %%eax
  "
  : : "a" (a), "D" (b), "c" (len), "d" (&sum) : "memory");
  return sum;
}
#endif

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
