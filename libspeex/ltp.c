/* Copyright (C) 2002 Jean-Marc Valin 
   File: ltp.c
   Lont-Term Prediction functions
*/

#include "ltp.h"


#define abs(x) ((x)<0 ? -(x) : (x))

/** Computes the open-loop pitch prediction. Returns pitch period and pitch gain */
int open_loop_ltp(float *x, int len, int start, int end, float *gain)
   /*  x:     time-domain signal (note, x[-end] must be valid)
       len:   length of the x signal
       start: smallest pitch period possible
       end:   largest pitch period possible
       gain:  return value for the pitch predictor gain
    */
{
   int i, period, best_period=0;
   float score, best_score=-1, corr, energy;
   
   for (period = start; period <= end; period++)
   {
      corr=energy=0;
      for (i=0;i<len;i++)
      {
         corr += x[i]*x[i-period];
         energy += x[i-period]*x[i-period];
      }
      score = corr*corr/(energy+1);
      if (score > best_score)
      {
         best_score=score;
         best_period=period;
         *gain = corr/(energy+1);
      }
      
   }
   return best_period;
}

/** Computes a 3-tap pitch predictor */
int three_tap_ltp(float *x, int len, int start, int end, float *gain)

   /*  x:     time-domain signal (note, x[-end] must be valid)
       len:   length of the x signal
       start: smallest pitch period possible
       end:   largest pitch period possible
       gain:  return value for the pitch predictor gain
    */
{
   int i, period, best_period=0;
   float total, score[3]={0,0,0}, best_score=-1, corr[3]={0,0,0}, energy[3]={0,0,0};
   float best_energy[3], best_corr[3], best_gain=0;
   float A[3][3];
   /*Check all periods */
   for (period = start; period <= end; period++)
   {
      corr[0]=energy[0]=0;
      for (i=0;i<len;i++)
      {
         corr[0] += x[i]*x[i-period];
         energy[0] += x[i-period]*x[i-period];
      }
      score[0] = corr[0]*corr[0]/(energy[0]+1);
      /* Smooth the "correlation score" */
      total = score[0]+2*score[1]+score[2];
      if (total > best_score && period >= start+3 && period <= end-3)
      {
         best_score=total;
         best_period=period;
         for (i=0;i<3;i++)
         {
            best_corr[i]=corr[i];
            best_energy[i]=energy[i];
         }
         best_gain = corr[1]/(energy[1]+1);
      }
      score[2]=score[1];
      score[1]=score[0];
      corr[2]=corr[1];
      corr[1]=corr[0];
      energy[2]=energy[1];
      energy[1]=energy[0];
   }
   /* build the correlation matrix */
   A[0][0]=best_energy[0]+1;
   A[1][1]=best_energy[1]+1;
   A[2][2]=best_energy[2]+1;
   A[0][1]=0;
   A[0][2]=0;
   for (i=0;i<len;i++)
   {
      A[0][1] += x[i-best_period+1]*x[i-best_period];
      A[0][2] += x[i-best_period+2]*x[i-best_period];
   }
   A[1][0]=A[0][1];
   A[2][0]=A[0][2];
   A[1][2]=A[0][1] - x[-best_period+1]*x[-best_period] 
                   + x[len-best_period]*x[len-best_period+1];
   A[2][1]=A[1][2];
   /*for (i=0;i<3;i++)
      printf ("%f %f %f\n", A[i][0], A[i][1], A[i][2]);
      printf ("\n%f %f %f\n", best_corr[0], best_corr[1], best_corr[2]);*/
   /*Solve the linear system to find gains*/
   {
      float tmp=A[1][0]/A[0][0];
      for (i=0;i<3;i++)
         A[1][i] -= tmp*A[0][i];
      best_corr[1] -= tmp*best_corr[0];

      tmp=A[2][0]/A[0][0];
      for (i=0;i<3;i++)
         A[2][i] -= tmp*A[0][i];
      best_corr[2] -= tmp*best_corr[0];
      
      tmp=A[2][1]/A[1][1];
      A[2][2] -= tmp*A[1][2];
      best_corr[2] -= tmp*best_corr[1];

      best_corr[2] /= A[2][2];
      best_corr[1] = (best_corr[1] - A[1][2]*best_corr[2])/A[1][1];
      best_corr[0] = (best_corr[0] - A[0][2]*best_corr[2] - A[0][1]*best_corr[1])/A[0][0];
      /*printf ("\n%f %f %f\n", best_corr[0], best_corr[1], best_corr[2]);*/

   }
   /* Put gains in right order */
   gain[0]=best_corr[2];gain[1]=best_corr[1];gain[2]=best_corr[0];
   
   /* Check that 3-tap predictor "makes sense" */
   /*if (!((abs(gain[0]) + abs(gain[1]) + abs(gain[2])<2) 
         && abs(gain[0]+gain[1]+gain[2]) < 1.2 
         && abs(gain[0]+gain[1]+gain[2]) > -.2 ))
   {
      gain[0]=0;gain[1]=best_gain;gain[2]=0;
      if (best_gain > 1.2)
         gain[1]=1.2;
      else if (best_gain<-.2)
         gain[1]=-.2;
      else 
         gain[1]=best_gain;
         }*/
   return best_period-2;
}


/** Finds the best 3-tap pitch predictor from a codebook*/
int ltp_closed_loop(float *x, int len, int start, int end, float *gain)

   /*  x:     time-domain signal (note, x[-end] must be valid)
       len:   length of the x signal
       start: smallest pitch period possible
       end:   largest pitch period possible
       gain:  return value for the pitch predictor gain
    */
{
   int i, period, best_period=0;
   float total, score[3]={0,0,0}, best_score=-1, corr[3]={0,0,0}, energy[3]={0,0,0};
   float best_energy[3], best_corr[3], best_gain=0;
   float A[3][3];
   /*Check all periods */
   for (period = start; period <= end; period++)
   {
      corr[0]=energy[0]=0;
      for (i=0;i<len;i++)
      {
         corr[0] += x[i]*x[i-period];
         energy[0] += x[i-period]*x[i-period];
      }
      score[0] = corr[0]*corr[0]/(energy[0]+1);
      /* Smooth the "correlation score" */
      total = score[0]+2*score[1]+score[2];
      if (total > best_score && period >= start+3 && period <= end-3)
      {
         best_score=total;
         best_period=period;
         for (i=0;i<3;i++)
         {
            best_corr[i]=corr[i];
            best_energy[i]=energy[i];
         }
         best_gain = corr[1]/(energy[1]+1);
      }
      score[2]=score[1];
      score[1]=score[0];
      corr[2]=corr[1];
      corr[1]=corr[0];
      energy[2]=energy[1];
      energy[1]=energy[0];
   }
   /* build the correlation matrix */
   A[0][0]=best_energy[0]+1;
   A[1][1]=best_energy[1]+1;
   A[2][2]=best_energy[2]+1;
   A[0][1]=0;
   A[0][2]=0;
   for (i=0;i<len;i++)
   {
      A[0][1] += x[i-best_period+1]*x[i-best_period];
      A[0][2] += x[i-best_period+2]*x[i-best_period];
   }
   A[1][0]=A[0][1];
   A[2][0]=A[0][2];
   A[1][2]=A[0][1] - x[-best_period+1]*x[-best_period] 
                   + x[len-best_period]*x[len-best_period+1];
   A[2][1]=A[1][2];

   {
      int j;
      float C[9];
      float *ptr=gain_cdbk_nb;
      int best_cdbk=0;
      float best_sum=0;
      C[0]=best_corr[2];
      C[1]=best_corr[1];
      C[2]=best_corr[0];
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
   }
   return best_period-2;
}



/** In place 3-tap pitch predictor (FIR)*/
void predictor_three_tap(float *x, int len, int period, float *gain)
{
   int i;
   for (i=len-1;i>=0;i--)
   {
      x[i] -= gain[0]*x[i-period] + gain[1]*x[i-period-1] + gain[2]*x[i-period-2];
   }
}

/** In place 3-tap inverse pitch predictor (IIR)*/
void inverse_three_tap(float *x, int len, int period, float *gain)
{
   int i;
   for (i=0;i<len;i++)
   {
      x[i] += gain[0]*x[i-period] + gain[1]*x[i-period-1] + gain[2]*x[i-period-2];
   }
}
