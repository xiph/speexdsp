/* Copyright (C) 2002 Jean-Marc Valin 
   File: ltp.c
   Lont-Term Prediction functions
*/


/*static float lp4[] = {0.00139, 0.00853, 0.01465, -0.00000, -0.04766, -0.10610, -0.11333, 0.00000, 0.25616, 0.59397, 0.88498, 1.00000, 0.88498, 0.59397, 0.25616,  0.00000, -0.11333, -0.10610, -0.04766, -0.00000, 0.01465, 0.00853, 0.00139};
 */


/* Computes the open-loop pitch prediction. Returns pitch period and pitch gain */
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
   float best_energy[3], best_corr[3];
   float A[3][3];
   for (period = start; period <= end; period++)
   {
      corr[0]=energy[0]=0;
      for (i=0;i<len;i++)
      {
         corr[0] += x[i]*x[i-period];
         energy[0] += x[i-period]*x[i-period];
      }
      score[0] = corr[0]*corr[0]/(energy[0]+1);
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
         *gain = corr[1]/(energy[1]+1);
      }
      corr[2]=corr[1];
      corr[1]=corr[0];
      energy[2]=energy[1];
      energy[1]=energy[0];
   }
   A[0][0]=energy[0];
   A[1][1]=energy[1];
   A[2][2]=energy[2];
   A[0][1]=0;
   A[0][2]=0;
   for (i=0;i<len;i++)
   {
      A[0][1] += x[i-best_period]*x[i-best_period-1];
      A[0][2] += x[i-best_period]*x[i-best_period-2];
   }
   A[1][0]=A[0][1];
   A[2][0]=A[0][2];
   A[1][2]=A[0][1] - x[-best_period]*x[-best_period-1] 
                   + x[len-best_period]*x[len-best_period-1];
   A[2][1]=A[1][2];
   return best_period-1;
}


