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

#include "ltp.h"
#include "cb_search.h"
#include "stack_alloc.h"
#include "filters.h"

#define abs(x) ((x)<0 ? -(x) : (x))



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
int   nsf                       /* Number of samples in subframe */
)
{
   int i,j;
   float tmp[3*nsf];
   float *x[3];
   float corr[3];
   float A[3][3];

   x[0]=tmp;
   x[1]=tmp+nsf;
   x[2]=tmp+2*nsf;

   /* Perform closed-loop 1-tap search*/
   overlap_cb_search(target, ak, awk1, awk2,
                     &exc[-end], end-start+1, gain, pitch, p,
                     nsf);
   /* Real pitch value */
   *pitch=end-*pitch;
   
   
   for (i=0;i<3;i++)
   {
      residue_zero(&exc[-*pitch-1+i],awk1,x[i],nsf,p);
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

   --*pitch;

   {
      float tmp1=0,tmp2=0;
      for (i=0;i<nsf;i++)
         tmp1+=target[i]*target[i];
      for (i=0;i<nsf;i++)
         tmp2+=(target[i]-gain[2]*x[0][i]-gain[1]*x[1][i]-gain[0]*x[2][i])
         * (target[i]-gain[2]*x[0][i]-gain[1]*x[1][i]-gain[0]*x[2][i]);
      printf ("prediction gain = %f\n",tmp1/(tmp2+1));
      return tmp1/(tmp2+1);
   }
}


/** Finds the best quantized 3-tap pitch predictor by analysis by synthesis */
float pitch_search_3tap(
float target[],                 /* Target vector */
float ak[],                     /* LPCs for this subframe */
float awk1[],                   /* Weighted LPCs #1 for this subframe */
float awk2[],                   /* Weighted LPCs #2 for this subframe */
float exc[],                    /* Overlapping codebook */
int   start,                    /* Smallest pitch value allowed */
int   end,                      /* Largest pitch value allowed */
float *gain,                    /* 3-tab gains of optimum entry */
int   *pitch,                   /* Index of optimum entry */
int   *gain_index,              /* Index of optimum gain */
int   p,                        /* Number of LPC coeffs */
int   nsf                       /* Number of samples in subframe */
)
{
   int i,j;
   float tmp[3*nsf];
   float *x[3];
   float corr[3];
   float A[3][3];

   x[0]=tmp;
   x[1]=tmp+nsf;
   x[2]=tmp+2*nsf;

   /* Perform closed-loop 1-tap search*/
   overlap_cb_search(target, ak, awk1, awk2,
                     &exc[-end], end-start+1, gain, pitch, p,
                     nsf);
   /* Real pitch value */
   *pitch=end-*pitch;
   
   
   for (i=0;i<3;i++)
   {
      residue_zero(&exc[-*pitch-1+i],awk1,x[i],nsf,p);
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
      float *ptr=gain_cdbk_nb;
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
      
      for (i=0;i<127;i++)
      {
         float sum=0;
         ptr = gain_cdbk_nb+12*i;
         for (j=0;j<9;j++)
            sum+=C[j]*ptr[j+3];
         if (sum>best_sum || i==0)
         {
            best_sum=sum;
            best_cdbk=i;
         }
      }
      gain[0] = gain_cdbk_nb[best_cdbk*12];
      gain[1] = gain_cdbk_nb[best_cdbk*12+1];
      gain[2] = gain_cdbk_nb[best_cdbk*12+2];
      *gain_index=best_cdbk;

   }
   --*pitch;

   {
      float tmp1=0,tmp2=0;
      for (i=0;i<nsf;i++)
         tmp1+=target[i]*target[i];
      for (i=0;i<nsf;i++)
         tmp2+=(target[i]-gain[2]*x[0][i]-gain[1]*x[1][i]-gain[0]*x[2][i])
         * (target[i]-gain[2]*x[0][i]-gain[1]*x[1][i]-gain[0]*x[2][i]);
      printf ("prediction gain = %f\n",tmp1/(tmp2+1));
   }
}
