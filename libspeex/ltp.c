/* Copyright (C) 2002 Jean-Marc Valin 
   File: ltp.c
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

#include <math.h>
#include <stdio.h>
#include "ltp.h"
#include "stack_alloc.h"
#include "filters.h"
#include "speex_bits.h"


static float inner_prod(float *x, float *y, int len)
{
   int i;
   float sum1=0,sum2=0,sum3=0,sum4=0;
   for (i=0;i<len;)
   {
      sum1 += x[i]*y[i];
      sum2 += x[i+1]*y[i+1];
      sum3 += x[i+2]*y[i+2];
      sum4 += x[i+3]*y[i+3];
      i+=4;
   }
   return sum1+sum2+sum3+sum4;
}

/*Original, non-optimized version*/
/*static float inner_prod(float *x, float *y, int len)
{
   int i;
   float sum=0;
   for (i=0;i<len;i++)
      sum += x[i]*y[i];
   return sum;
}
*/


void open_loop_nbest_pitch(float *sw, int start, int end, int len, int *pitch, float *gain, int N, float *stack)
{
   int i,j,k;
   float corr=0;
   float energy;
   float score=0, *best_score;
   float e0;

   best_score = PUSH(stack,N);
   for (i=0;i<N;i++)
        best_score[i]=-1;
   energy=inner_prod(sw-start, sw-start, len);
   e0=inner_prod(sw, sw, len);
   for (i=start;i<=end;i++)
   {
      if (0&&score < .3*best_score[N-1])
         score = .9*best_score[N-1];
      else
      {
         corr=inner_prod(sw, sw-i, len);
         score=corr*corr/(energy+1);
      }
      if (score>best_score[N-1])
      {
         float g1, g;
         g1 = corr/(energy+10);
         g = sqrt(g1*corr/(e0+10));
         if (g>g1)
            g=g1;
         if (g<0)
            g=0;
         for (j=0;j<N;j++)
         {
            if (score > best_score[j])
            {
               for (k=N-1;k>j;k--)
               {
                  best_score[k]=best_score[k-1];
                  pitch[k]=pitch[k-1];
                  gain[k] = gain[k-1];
               }
               best_score[j]=score;
               pitch[j]=i;
               gain[j]=g;
               break;
            }
         }
      }
      /* Update energy for next pitch*/
      energy+=sw[-i-1]*sw[-i-1] - sw[-i+len-1]*sw[-i+len-1];
   }

   POP(stack);
}




/** Finds the best quantized 3-tap pitch predictor by analysis by synthesis */
float pitch_gain_search_3tap(
float target[],                 /* Target vector */
float ak[],                     /* LPCs for this subframe */
float awk1[],                   /* Weighted LPCs #1 for this subframe */
float awk2[],                   /* Weighted LPCs #2 for this subframe */
float exc[],                    /* Excitation */
void *par,
int   pitch,                    /* Pitch value */
int   p,                        /* Number of LPC coeffs */
int   nsf,                      /* Number of samples in subframe */
SpeexBits *bits,
float *stack,
float *exc2,
int  *cdbk_index
)
{
   int i,j;
   float *tmp, *tmp2;
   float *x[3];
   float *e[3];
   float corr[3];
   float A[3][3];
   float gain[3];
   int   gain_cdbk_size;
   float *gain_cdbk;
   float err1,err2;
   ltp_params *params;
   params = (ltp_params*) par;
   gain_cdbk=params->gain_cdbk;
   gain_cdbk_size=1<<params->gain_bits;
   tmp = PUSH(stack, 3*nsf);
   tmp2 = PUSH(stack, 3*nsf);

   x[0]=tmp;
   x[1]=tmp+nsf;
   x[2]=tmp+2*nsf;

   e[0]=tmp2;
   e[1]=tmp2+nsf;
   e[2]=tmp2+2*nsf;
   
   for (i=0;i<3;i++)
   {
      int pp=pitch+1-i;
      for (j=0;j<nsf;j++)
      {
         if (j-pp<0)
            e[i][j]=exc2[j-pp];
         else if (j-pp-pitch<0)
            e[i][j]=exc2[j-pp-pitch];
         else
            e[i][j]=0;
      }

      if (p==10)
      {
         syn_percep_zero(e[i], ak, awk1, awk2, x[i], nsf, p);
      } else {
         residue_zero(e[i],awk1,x[i],nsf,p);
         syn_filt_zero(x[i],ak,x[i],nsf,p);
         syn_filt_zero(x[i],awk2,x[i],nsf,p);
      }
   }

   for (i=0;i<3;i++)
      corr[i]=inner_prod(x[i],target,nsf);
   
   for (i=0;i<3;i++)
      for (j=0;j<=i;j++)
         A[i][j]=A[j][i]=inner_prod(x[i],x[j],nsf);
   
   {
      int j;
      float C[9];
      float *ptr=gain_cdbk;
      int best_cdbk=0;
      float best_sum=0;
      C[0]=corr[2];
      C[1]=corr[1];
      C[2]=corr[0];
      C[3]=A[1][2];
      C[4]=A[0][1];
      C[5]=A[0][2];
      C[6]=A[2][2];
      C[7]=A[1][1];
      C[8]=A[0][0];
      
      for (i=0;i<gain_cdbk_size;i++)
      {
         float sum=0;
         ptr = gain_cdbk+12*i;
         for (j=0;j<9;j++)
            sum+=C[j]*ptr[j+3];
         if (0) {
            float tot=ptr[0]+ptr[1]+ptr[2];
            if (tot < 1.1)
               sum *= 1+.15*tot;
         }
         if (sum>best_sum || i==0)
         {
            best_sum=sum;
            best_cdbk=i;
         }
      }
      gain[0] = gain_cdbk[best_cdbk*12];
      gain[1] = gain_cdbk[best_cdbk*12+1];
      gain[2] = gain_cdbk[best_cdbk*12+2];

      *cdbk_index=best_cdbk;
   }
   /* Calculate gains by matrix inversion... (unquantized) */
   if (0) {
      float tmp;
      float B[3][3];
      A[0][0]+=1;
      A[1][1]+=1;
      A[2][2]+=1;
      
      for (i=0;i<3;i++)
         for (j=0;j<3;j++)
            B[i][j]=A[i][j];


      tmp=A[1][0]/A[0][0];
      for (i=0;i<3;i++)
         A[1][i] -= tmp*A[0][i];
      corr[1] -= tmp*corr[0];

      tmp=A[2][0]/A[0][0];
      for (i=0;i<3;i++)
         A[2][i] -= tmp*A[0][i];
      corr[2] -= tmp*corr[0];
      
      tmp=A[2][1]/A[1][1];
      A[2][2] -= tmp*A[1][2];
      corr[2] -= tmp*corr[1];

      corr[2] /= A[2][2];
      corr[1] = (corr[1] - A[1][2]*corr[2])/A[1][1];
      corr[0] = (corr[0] - A[0][2]*corr[2] - A[0][1]*corr[1])/A[0][0];
      /*printf ("\n%f %f %f\n", best_corr[0], best_corr[1], best_corr[2]);*/

   
      /* Put gains in right order */
      gain[0]=corr[2];gain[1]=corr[1];gain[2]=corr[0];

      {
         float gain_sum = gain[0]+gain[1]+gain[2];
         if (fabs(gain_sum)>2.5)
         {
            float fact = 2.5/gain_sum;
            for (i=0;i<3;i++)
               gain[i]*=fact;
         }
      }
      
   }
   
   for (i=0;i<nsf;i++)
      exc[i]=gain[0]*e[2][i]+gain[1]*e[1][i]+gain[2]*e[0][i];
#ifdef DEBUG
   printf ("3-tap pitch = %d, gains = [%f %f %f]\n",pitch, gain[0], gain[1], gain[2]);
#endif
   
   err1=0;
   err2=0;
   for (i=0;i<nsf;i++)
      err1+=target[i]*target[i];
   for (i=0;i<nsf;i++)
      err2+=(target[i]-gain[2]*x[0][i]-gain[1]*x[1][i]-gain[0]*x[2][i])
      * (target[i]-gain[2]*x[0][i]-gain[1]*x[1][i]-gain[0]*x[2][i]);
#ifdef DEBUG
   printf ("prediction gain = %f\n",err1/(err2+1));
#endif

   POP(stack);
   POP(stack);
   return err2;
}


/** Finds the best quantized 3-tap pitch predictor by analysis by synthesis */
int pitch_search_3tap(
float target[],                 /* Target vector */
float *sw,
float ak[],                     /* LPCs for this subframe */
float awk1[],                   /* Weighted LPCs #1 for this subframe */
float awk2[],                   /* Weighted LPCs #2 for this subframe */
float exc[],                    /* Excitation */
void *par,
int   start,                    /* Smallest pitch value allowed */
int   end,                      /* Largest pitch value allowed */
float pitch_coef,               /* Voicing (pitch) coefficient */
int   p,                        /* Number of LPC coeffs */
int   nsf,                      /* Number of samples in subframe */
SpeexBits *bits,
float *stack,
float *exc2,
int complexity
)
{
   int i,j;
   int cdbk_index, pitch=0, best_gain_index=0;
   float *best_exc;
   int best_pitch=0;
   float err, best_err=-1;
   int N=3;
   ltp_params *params;
   int *nbest;
   float *gains;

   N=complexity;
   if (N<1)
      N=1;
   if (N>10)
      N=10;

   nbest=(int*)PUSH(stack, N);
   gains = PUSH(stack, N);
   params = (ltp_params*) par;
   
   best_exc=PUSH(stack,nsf);
   
   if (N>end-start+1)
      N=end-start+1;
   open_loop_nbest_pitch(sw, start, end, nsf, nbest, gains, N, stack);
   for (i=0;i<N;i++)
   {
      pitch=nbest[i];
      for (j=0;j<nsf;j++)
         exc[j]=0;
      err=pitch_gain_search_3tap(target, ak, awk1, awk2, exc, par, pitch, p, nsf,
                                 bits, stack, exc2, &cdbk_index);
      if (err<best_err || best_err<0)
      {
         for (j=0;j<nsf;j++)
            best_exc[j]=exc[j];
         best_err=err;
         best_pitch=pitch;
         best_gain_index=cdbk_index;
      }
   }
   
   /*printf ("pitch: %d %d\n", best_pitch, best_gain_index);*/
   speex_bits_pack(bits, best_pitch-start, params->pitch_bits);
   speex_bits_pack(bits, best_gain_index, params->gain_bits);
   /*printf ("encode pitch: %d %d\n", best_pitch, best_gain_index);*/
   for (i=0;i<nsf;i++)
      exc[i]=best_exc[i];

   POP(stack);
   POP(stack);
   POP(stack);
   return pitch;
}


void pitch_unquant_3tap(
float exc[],                    /* Excitation */
int   start,                    /* Smallest pitch value allowed */
int   end,                      /* Largest pitch value allowed */
float pitch_coef,               /* Voicing (pitch) coefficient */
void *par,
int   nsf,                      /* Number of samples in subframe */
int *pitch_val,
float *gain_val,
SpeexBits *bits,
float *stack,
int lost)
{
   int i;
   int pitch;
   int gain_index;
   float gain[3];
   float *gain_cdbk;
   ltp_params *params;
   params = (ltp_params*) par;
   gain_cdbk=params->gain_cdbk;

   pitch = speex_bits_unpack_unsigned(bits, params->pitch_bits);
   pitch += start;
   gain_index = speex_bits_unpack_unsigned(bits, params->gain_bits);
   /*printf ("decode pitch: %d %d\n", pitch, gain_index);*/
   gain[0] = gain_cdbk[gain_index*12];
   gain[1] = gain_cdbk[gain_index*12+1];
   gain[2] = gain_cdbk[gain_index*12+2];
   if (lost)
   {
      float gain_sum;
      gain_sum = fabs(gain[0])+fabs(gain[1])+fabs(gain[2]);
      if (gain_sum>.95)
      {
         float fact = .95/gain_sum;
         for (i=0;i<3;i++)
            gain[i]*=fact;
      }
   }

   *pitch_val = pitch;
   /**gain_val = gain[0]+gain[1]+gain[2];*/
   gain_val[0]=gain[0];
   gain_val[1]=gain[1];
   gain_val[2]=gain[2];

#ifdef DEBUG
   printf ("unquantized pitch: %d %f %f %f\n", pitch, gain[0], gain[1], gain[2]);
#endif
   for(i=0;i<nsf;i++)
      exc[i]=0;

   {
      float *e[3];
      float *tmp2;
      tmp2=PUSH(stack, 3*nsf);
      e[0]=tmp2;
      e[1]=tmp2+nsf;
      e[2]=tmp2+2*nsf;
      
      for (i=0;i<3;i++)
      {
         int j;
         int pp=pitch+1-i;
         for (j=0;j<nsf;j++)
         {
            if (j-pp<0)
               e[i][j]=exc[j-pp];
            else
               e[i][j]=exc[j-pp-pitch];
         }
      }
      for (i=0;i<nsf;i++)
         exc[i]=gain[0]*e[2][i]+gain[1]*e[1][i]+gain[2]*e[0][i];
      
      POP(stack);
   }
}


/** Forced pitch delay and gain */
int forced_pitch_quant(
float target[],                 /* Target vector */
float *sw,
float ak[],                     /* LPCs for this subframe */
float awk1[],                   /* Weighted LPCs #1 for this subframe */
float awk2[],                   /* Weighted LPCs #2 for this subframe */
float exc[],                    /* Excitation */
void *par,
int   start,                    /* Smallest pitch value allowed */
int   end,                      /* Largest pitch value allowed */
float pitch_coef,               /* Voicing (pitch) coefficient */
int   p,                        /* Number of LPC coeffs */
int   nsf,                      /* Number of samples in subframe */
SpeexBits *bits,
float *stack,
float *exc2,
int complexity
)
{
   int i;
   if (pitch_coef>.9)
      pitch_coef=.9;
   for (i=0;i<nsf;i++)
   {
      exc[i]=exc[i-start]*pitch_coef;
   }
   return start;
}

/** Unquantize forced pitch delay and gain */
void forced_pitch_unquant(
float exc[],                    /* Excitation */
int   start,                    /* Smallest pitch value allowed */
int   end,                      /* Largest pitch value allowed */
float pitch_coef,               /* Voicing (pitch) coefficient */
void *par,
int   nsf,                      /* Number of samples in subframe */
int *pitch_val,
float *gain_val,
SpeexBits *bits,
float *stack,
int lost)
{
   int i;
   /*pitch_coef=.9;*/
   if (pitch_coef>.9)
      pitch_coef=.9;
   for (i=0;i<nsf;i++)
   {
      exc[i]=exc[i-start]*pitch_coef;
   }
   *pitch_val = start;
   gain_val[0]=gain_val[2]=0;
   gain_val[1] = pitch_coef;
}
