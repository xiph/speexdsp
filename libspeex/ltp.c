/* Copyright (C) 2002 Jean-Marc Valin 
   File: ltp.c
   Lont-Term Prediction functions
*/

int open_loop_ltp(float *x, int len, int start, int end, float *gain)
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
