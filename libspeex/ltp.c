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
#include "cb_search.h"
#include "stack_alloc.h"
#include "filters.h"
#include "speex_bits.h"

#define abs(x) ((x)<0 ? -(x) : (x))

void open_loop_pitch(float *sw, int start, int end, int len, int *pitch)
{
   int i;
   float e0, corr, energy, best_gain=0, pred_gain, best_corr=0, best_energy=0;
   float score, best_score=-1;
   e0=xcorr(sw, sw, len);
   energy=xcorr(sw-start, sw-start, len);
   for (i=start;i<=end;i++)
   {
      corr=xcorr(sw, sw-i, len);
      score=corr*corr/(energy+1);
      if (score>best_score)
      {
         /*if ((abs(i-2**pitch)>4 && abs(i-3**pitch)>6) || score>0.9*best_score)*/
         {
         best_score=score;
         best_gain=corr/(energy+1);
         best_corr=corr;
         best_energy=energy;
         *pitch=i;
         }
      }

      /* Update energy for next pitch*/
      energy+=sw[-i-1]*sw[-i-1] - sw[-i+len]*sw[-i+len];
   }
   pred_gain=e0/(1+fabs(e0+best_gain*best_gain*best_energy-2*best_gain*best_corr));
#ifdef DEBUG
   printf ("pred = %f\n", pred_gain);
#endif
}

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
)
{
   int i, j, size, filt_size, base, frac=0, best_cor=0;
   float *oexc_mem, *oexc, *exc_ptr, *fexc, *f, frac_pitch=0, best_score=-1, best_gain=0;
   float sc;
   float corr[3];
   float A[3][3];
#if 1
   sc = overlap_cb_search(target, ak, awk1, awk2,
                     &exc[-end], end-start+1, gain, pitch, p,
                     nsf, stack);
                     *pitch=end-*pitch;
#ifdef DEBUG
   printf ("ol score: %d %f\n", *pitch, sc);
#endif
#endif
   base=*pitch;
   exc_ptr=exc-*pitch;
   size = fact*nsf + filt_side*2 + 16*fact;
   filt_size = 2*filt_side+1;
   oexc_mem = PUSH(stack, size);
   oexc=oexc_mem+filt_side;
   fexc = PUSH(stack, size/fact);
   f=filt+filt_side;

   for(i=0;i<size;i++)
      oexc_mem[i]=0;
   for (i=-8;i<nsf+8;i++)
   {
      for (j=-filt_side;j<=filt_side;j++)
         oexc[fact*(i+8)+j] += fact*exc_ptr[i]*f[j];
   }

   /*for (i=0;i<size;i++)
     printf ("%f ", oexc_mem[i]);
   printf ("eee\n");
   */
   for (j=0;j<fact;j++)
   {
      int correction;
      float score;
      for (i=0;i<size/fact;i++)
         fexc[i]=oexc[fact*i+j];
      score=overlap_cb_search(target, ak, awk1, awk2,
                        fexc, 16, gain, &correction, p,
                        nsf, stack);
      if (score>best_score)
      {
         best_cor = correction;
         *pitch = base+8-correction;
         frac = j;
         best_gain=*gain;
         frac_pitch = *pitch-(j/(float)fact);
         best_score=score;
      }
#ifdef DEBUG
      printf ("corr: %d %d %f\n", correction, *pitch, score);
#endif
   }
   /*for (i=0;i<nsf;i++)
     printf ("%f ", oexc[4*(i+8)]);
   printf ("aaa\n");
   for (i=0;i<nsf;i++)
     printf ("%f ", exc[i-base]);
     printf ("bbb\n");*/

   /*if (best_gain>1.2)
      best_gain=1.2;
   if (best_gain<-.2)
   best_gain=-.2;*/
   for (i=0;i<nsf;i++)
     exc[i]=best_gain*oexc[fact*(best_cor+i)+frac];
   
   {
      float *x[3];
      x[0] = PUSH(stack, nsf);
      x[1] = PUSH(stack, nsf);
      x[2] = PUSH(stack, nsf);
      
      for (j=0;j<3;j++)
         for (i=0;i<nsf;i++)
            x[j][i]=oexc[fact*(best_cor+i)+frac+fact*(j-1)];


   for (i=0;i<3;i++)
   {
      residue_zero(x[i],awk1,x[i],nsf,p);
      syn_filt_zero(x[i],ak,x[i],nsf,p);
      syn_filt_zero(x[i],awk2,x[i],nsf,p);
   }

   for (i=0;i<3;i++)
      corr[i]=xcorr(x[i],target,nsf);
   
   for (i=0;i<3;i++)
      for (j=0;j<=i;j++)
         A[i][j]=A[j][i]=xcorr(x[i],x[j],nsf);
   /*for (i=0;i<3;i++)
   {
      for (j=0;j<3;j++)
         printf ("%f ", A[i][j]);
      printf ("\n");
      }*/
   A[0][0]+=1;
   A[1][1]+=1;
   A[2][2]+=1;
   {
      float tmp=A[1][0]/A[0][0];
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

   }
   /* Put gains in right order */
   /*gain[0]=corr[2];gain[1]=corr[1];gain[2]=corr[0];*/
   gain[0]=corr[0];gain[1]=corr[1];gain[2]=corr[2];

   for (i=0;i<nsf;i++)
      exc[i]=gain[0]*oexc[fact*(best_cor+i)+frac-fact] + 
             gain[1]*oexc[fact*(best_cor+i)+frac] + 
             gain[2]*oexc[fact*(best_cor+i)+frac+fact];
   
   /*for (i=0;i<nsf;i++)
     exc[i]=best_gain*x[1][i];*/
   /*for (i=0;i<nsf;i++)
     exc[i]=best_gain*oexc[fact*(best_cor+i)+frac];*/
#ifdef DEBUG
   printf ("frac gains: %f %f %f\n", gain[0], gain[1], gain[2]);
#endif
      POP(stack);
      POP(stack);
      POP(stack);
   }
#ifdef DEBUG
   printf ("frac pitch = %f %f\n", frac_pitch, best_score);
#endif
   POP(stack);
   POP(stack);

}

void pitch_search_3tap_unquant(
float target[],                 /* Target vector */
float ak[],                     /* LPCs for this subframe */
float awk1[],                   /* Weighted LPCs #1 for this subframe */
float awk2[],                   /* Weighted LPCs #2 for this subframe */
float exc[],                    /* Excitation */
void *par,
int   start,                    /* Smallest pitch value allowed */
int   end,                      /* Largest pitch value allowed */
int   p,                        /* Number of LPC coeffs */
int   nsf,                      /* Number of samples in subframe */
FrameBits *bits,
float *stack
)
{
   int i,j;
   float *tmp;
   float *x[3];
   float corr[3];
   float A[3][3];
   int pitch;
   float gain[3];
   tmp = PUSH(stack, 3*nsf);
   x[0]=tmp;
   x[1]=tmp+nsf;
   x[2]=tmp+2*nsf;

   /* Perform closed-loop 1-tap search*/
   overlap_cb_search(target, ak, awk1, awk2,
                     &exc[-end], end-start+1, gain, &pitch, p,
                     nsf, stack);
   /* Real pitch value */
   pitch=end-pitch;
   
   
   for (i=0;i<3;i++)
   {
      residue_zero(&exc[-pitch-1+i],awk1,x[i],nsf,p);
      syn_filt_zero(x[i],ak,x[i],nsf,p);
      syn_filt_zero(x[i],awk2,x[i],nsf,p);
   }

   for (i=0;i<3;i++)
      corr[i]=xcorr(x[i],target,nsf);
   
   for (i=0;i<3;i++)
      for (j=0;j<=i;j++)
         A[i][j]=A[j][i]=xcorr(x[i],x[j],nsf);
   A[0][0]+=1;
   A[1][1]+=1;
   A[2][2]+=1;
   {
      float tmp=A[1][0]/A[0][0];
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

   }
   /* Put gains in right order */
   gain[0]=corr[2];gain[1]=corr[1];gain[2]=corr[0];

   for (i=nsf-1;i>=0;i--)
      exc[i]=gain[0]*exc[i-pitch+1]+gain[1]*exc[i-pitch]+gain[2]*exc[i-pitch-1];
#if 0
   if (0){
      float tmp1=0,tmp2=0;
      for (i=0;i<nsf;i++)
         tmp1+=target[i]*target[i];
      for (i=0;i<nsf;i++)
         tmp2+=(target[i]-gain[2]*x[0][i]-gain[1]*x[1][i]-gain[0]*x[2][i])
         * (target[i]-gain[2]*x[0][i]-gain[1]*x[1][i]-gain[0]*x[2][i]);
#ifdef DEBUG
      printf ("prediction gain = %f\n",tmp1/(tmp2+1));
#endif
      return tmp1/(tmp2+1);
   }
#endif
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
FrameBits *bits,
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
         else
            e[i][j]=exc2[j-pp-pitch];
      }
      residue_zero(e[i],awk1,x[i],nsf,p);
      syn_filt_zero(x[i],ak,x[i],nsf,p);
      syn_filt_zero(x[i],awk2,x[i],nsf,p);
   }

   for (i=0;i<3;i++)
      corr[i]=xcorr(x[i],target,nsf);
   
   for (i=0;i<3;i++)
      for (j=0;j<=i;j++)
         A[i][j]=A[j][i]=xcorr(x[i],x[j],nsf);
   
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
   
   for (i=0;i<nsf;i++)
      exc[i]=gain[0]*e[2][i]+gain[1]*e[1][i]+gain[2]*e[0][i];
#ifdef DEBUG
   printf ("3-tap pitch = %d, gains = [%f %f %f]\n",pitch, gain[0], gain[1], gain[2]);
#endif
   {
      float tmp1=0,tmp2=0;
      for (i=0;i<nsf;i++)
         tmp1+=target[i]*target[i];
      for (i=0;i<nsf;i++)
         tmp2+=(target[i]-gain[2]*x[0][i]-gain[1]*x[1][i]-gain[0]*x[2][i])
         * (target[i]-gain[2]*x[0][i]-gain[1]*x[1][i]-gain[0]*x[2][i]);
#ifdef DEBUG
      printf ("prediction gain = %f\n",tmp1/(tmp2+1));
#endif
      return tmp2;
   }
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
int   p,                        /* Number of LPC coeffs */
int   nsf,                      /* Number of samples in subframe */
FrameBits *bits,
float *stack,
float *exc2
)
{
   int i,j;
   int cdbk_index, pitch, ol_pitch, best_gain_index=0;
   float *best_exc;
   int best_pitch=0;
   float err, best_err=-1;
   int N=1;
   ltp_params *params;
   params = (ltp_params*) par;

   best_exc=PUSH(stack,nsf);

   open_loop_pitch(sw, start, nsf, nsf, &ol_pitch);

   if (ol_pitch-N<start)
      ol_pitch=start+N;
   if (ol_pitch+N>end)
      ol_pitch=end-N;

   for (pitch=ol_pitch-N; pitch<=ol_pitch+N; pitch++)
   {
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


   open_loop_pitch(sw, nsf, end, nsf, &ol_pitch);

   for (pitch=ol_pitch-N; pitch<=ol_pitch+N; pitch++)
   {
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
   return pitch;
}


void pitch_unquant_3tap(
float exc[],                    /* Excitation */
int   start,                    /* Smallest pitch value allowed */
int   end,                      /* Largest pitch value allowed */
void *par,
int   nsf,                      /* Number of samples in subframe */
FrameBits *bits,
float *stack
)
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
