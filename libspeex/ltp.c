/* Copyright (C) 2002 Jean-Marc Valin 
   File: ltp.c
   Long-Term Prediction functions

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:
   
   - Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
   
   - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
   
   - Neither the name of the Xiph.org Foundation nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.
   
   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE FOUNDATION OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include "ltp.h"
#include "stack_alloc.h"
#include "filters.h"
#include <speex/speex_bits.h>
#include "math_approx.h"

#include <stdio.h>

#ifdef _USE_SSE
#include "ltp_sse.h"
#else
static spx_word32_t inner_prod(const spx_word16_t *x, const spx_word16_t *y, int len)
{
   int i;
   spx_word32_t sum=0;
   for (i=0;i<len;i+=4)
   {
      spx_word32_t part=0;
      part = MAC16_16(part,x[i],y[i]);
      part = MAC16_16(part,x[i+1],y[i+1]);
      part = MAC16_16(part,x[i+2],y[i+2]);
      part = MAC16_16(part,x[i+3],y[i+3]);
      sum = ADD32(sum,SHR(part,6));
   }
   return sum;
}

static void pitch_xcorr(const spx_word16_t *_x, const spx_word16_t *_y, spx_word32_t *corr, int len, int nb_pitch, char *stack)
{
   int i;
   for (i=0;i<nb_pitch;i++)
   {
      /* Compute correlation*/
      corr[nb_pitch-1-i]=inner_prod(_x, _y+i, len);
   }

}

#endif

void open_loop_nbest_pitch(spx_sig_t *sw, int start, int end, int len, int *pitch, spx_word16_t *gain, int N, char *stack)
{
   int i,j,k;
   spx_word32_t *best_score;
   spx_word32_t e0;
   spx_word32_t *corr, *energy;
   spx_word32_t *score;
   spx_word16_t *swn;

   best_score = PUSH(stack,N, spx_word32_t);
   corr = PUSH(stack,end-start+1, spx_word32_t);
   energy = PUSH(stack,end-start+2, spx_word32_t);
   score = PUSH(stack,end-start+1, spx_word32_t);

#ifdef FIXED_POINT
   swn = PUSH(stack, end+len, spx_word16_t);
   normalize16(sw-end, swn, 16384, end+len);
   swn += end;
#else
   swn = sw;
#endif

   for (i=0;i<N;i++)
   {
        best_score[i]=-1;
        gain[i]=0;
        pitch[i]=start;
   }


   energy[0]=inner_prod(swn-start, swn-start, len);
   e0=inner_prod(swn, swn, len);
   for (i=start;i<=end;i++)
   {
      /* Update energy for next pitch*/
      energy[i-start+1] = energy[i-start] + SHR(MULT16_16(swn[-i-1],swn[-i-1]),6) - SHR(MULT16_16(swn[-i+len-1],swn[-i+len-1]),6);
   }
   for (i=start;i<=end;i++)
   {
      corr[i-start]=0;
      score[i-start]=0;
   }

   pitch_xcorr(swn, swn-end, corr, len, end-start+1, stack);

#ifdef FIXED_POINT
   {
      spx_word16_t *corr16;
      spx_word16_t *ener16;
      corr16 = PUSH(stack, end-start+1, spx_word16_t);
      ener16 = PUSH(stack, end-start+1, spx_word16_t);
      normalize16(corr, corr16, 16384, end-start+1);
      normalize16(energy, ener16, 16384, end-start+1);

      for (i=start;i<=end;i++)
      {
         spx_word16_t g;
         spx_word32_t tmp;
         tmp = corr16[i-start];
         if (SHR(corr16[i-start],4)>ener16[i-start])
            tmp = SHL((spx_word32_t)ener16[i-start],14);
         else if (-SHR(corr16[i-start],4)>ener16[i-start])
            tmp = -SHL((spx_word32_t)ener16[i-start],14);
         else
            tmp = SHL(tmp,10);
         g = DIV32_16(tmp, 8+ener16[i-start]);
         score[i-start] = MULT16_16(corr16[i-start],g);
      }
   }
#else
   for (i=start;i<=end;i++)
   {
      float g = corr[i-start]/(1+energy[i-start]);
      if (g>16)
         g = 16;
      else if (g<-16)
         g = -16;
      score[i-start] = g*corr[i-start];
   }
#endif

   /* Extract best scores */
   for (i=start;i<=end;i++)
   {
      if (score[i-start]>best_score[N-1])
      {
         for (j=0;j<N;j++)
         {
            if (score[i-start] > best_score[j])
            {
               for (k=N-1;k>j;k--)
               {
                  best_score[k]=best_score[k-1];
                  pitch[k]=pitch[k-1];
               }
               best_score[j]=score[i-start];
               pitch[j]=i;
               break;
            }
         }
      }
   }

   /* Compute open-loop gain */
   for (j=0;j<N;j++)
   {
      spx_word16_t g;
      i=pitch[j];
      g = DIV32(corr[i-start], 10+SHR(MULT16_16(spx_sqrt(e0),spx_sqrt(energy[i-start])),6));
      /* FIXME: g = max(g,corr/energy) */
      if (g<0)
         g = 0;
      gain[j]=g;
   }
}


/** Finds the best quantized 3-tap pitch predictor by analysis by synthesis */
static spx_word64_t pitch_gain_search_3tap(
spx_sig_t target[],                 /* Target vector */
spx_coef_t ak[],                     /* LPCs for this subframe */
spx_coef_t awk1[],                   /* Weighted LPCs #1 for this subframe */
spx_coef_t awk2[],                   /* Weighted LPCs #2 for this subframe */
spx_sig_t exc[],                    /* Excitation */
const void *par,
int   pitch,                    /* Pitch value */
int   p,                        /* Number of LPC coeffs */
int   nsf,                      /* Number of samples in subframe */
SpeexBits *bits,
char *stack,
spx_sig_t *exc2,
spx_sig_t *r,
spx_sig_t *new_target,
int  *cdbk_index,
int cdbk_offset
)
{
   int i,j;
   spx_sig_t *tmp, *tmp2;
   spx_sig_t *x[3];
   spx_sig_t *e[3];
   spx_word32_t corr[3];
   spx_word32_t A[3][3];
   int   gain_cdbk_size;
   const signed char *gain_cdbk;
   spx_word16_t gain[3];
   spx_word64_t err;

   ltp_params *params;
   params = (ltp_params*) par;
   gain_cdbk_size = 1<<params->gain_bits;
   gain_cdbk = params->gain_cdbk + 3*gain_cdbk_size*cdbk_offset;
   tmp = PUSH(stack, 3*nsf, spx_sig_t);
   tmp2 = PUSH(stack, 3*nsf, spx_sig_t);

   x[0]=tmp;
   x[1]=tmp+nsf;
   x[2]=tmp+2*nsf;

   e[0]=tmp2;
   e[1]=tmp2+nsf;
   e[2]=tmp2+2*nsf;
   for (i=2;i>=0;i--)
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

      if (i==2)
         syn_percep_zero(e[i], ak, awk1, awk2, x[i], nsf, p, stack);
      else {
         for (j=0;j<nsf-1;j++)
            x[i][j+1]=x[i+1][j];
         x[i][0]=0;
         for (j=0;j<nsf;j++)
         {
            x[i][j]=ADD32(x[i][j],SHL(MULT16_32_Q15(r[j], e[i][0]),1));
         }
      }
   }

#ifdef FIXED_POINT
   {
      /* If using fixed-point, we need to normalize the signals first */
      spx_word16_t *y[3];
      spx_word16_t *t;

      spx_sig_t max_val=1;
      int sig_shift;
      
      y[0] = PUSH(stack, nsf, spx_word16_t);
      y[1] = PUSH(stack, nsf, spx_word16_t);
      y[2] = PUSH(stack, nsf, spx_word16_t);
      t = PUSH(stack, nsf, spx_word16_t);
      for (j=0;j<3;j++)
      {
         for (i=0;i<nsf;i++)
         {
            spx_sig_t tmp = x[j][i];
            if (tmp<0)
               tmp = -tmp;
            if (tmp > max_val)
               max_val = tmp;
         }
      }
      for (i=0;i<nsf;i++)
      {
         spx_sig_t tmp = target[i];
         if (tmp<0)
            tmp = -tmp;
         if (tmp > max_val)
            max_val = tmp;
      }

      sig_shift=0;
      while (max_val>16384)
      {
         sig_shift++;
         max_val >>= 1;
      }

      for (j=0;j<3;j++)
      {
         for (i=0;i<nsf;i++)
         {
            y[j][i] = SHR(x[j][i],sig_shift);
         }
      }     
      for (i=0;i<nsf;i++)
      {
         t[i] = SHR(target[i],sig_shift);
      }

      for (i=0;i<3;i++)
         corr[i]=inner_prod(y[i],t,nsf);
      
      for (i=0;i<3;i++)
         for (j=0;j<=i;j++)
            A[i][j]=A[j][i]=inner_prod(y[i],y[j],nsf);
   }
#else
   {
      for (i=0;i<3;i++)
         corr[i]=inner_prod(x[i],target,nsf);
      
      for (i=0;i<3;i++)
         for (j=0;j<=i;j++)
            A[i][j]=A[j][i]=inner_prod(x[i],x[j],nsf);
   }
#endif

   {
      spx_word32_t C[9];
      const signed char *ptr=gain_cdbk;
      int best_cdbk=0;
      spx_word32_t best_sum=0;
      C[0]=corr[2];
      C[1]=corr[1];
      C[2]=corr[0];
      C[3]=A[1][2];
      C[4]=A[0][1];
      C[5]=A[0][2];
      C[6]=A[2][2];
      C[7]=A[1][1];
      C[8]=A[0][0];
#ifndef FIXED_POINT
      C[6]*=.5;
      C[7]*=.5;
      C[8]*=.5;
#endif
      for (i=0;i<gain_cdbk_size;i++)
      {
         spx_word32_t sum=0;
         spx_word16_t g0,g1,g2;
         ptr = gain_cdbk+3*i;
         g0=ptr[0]+32;
         g1=ptr[1]+32;
         g2=ptr[2]+32;

         sum = ADD32(sum,MULT16_32_Q14(MULT16_16_16(g0,64),C[0]));
         sum = ADD32(sum,MULT16_32_Q14(MULT16_16_16(g1,64),C[1]));
         sum = ADD32(sum,MULT16_32_Q14(MULT16_16_16(g2,64),C[2]));
         sum = SUB32(sum,MULT16_32_Q14(MULT16_16_16(g0,g1),C[3]));
         sum = SUB32(sum,MULT16_32_Q14(MULT16_16_16(g2,g1),C[4]));
         sum = SUB32(sum,MULT16_32_Q14(MULT16_16_16(g2,g0),C[5]));
         sum = SUB32(sum,MULT16_32_Q15(MULT16_16_16(g0,g0),C[6]));
         sum = SUB32(sum,MULT16_32_Q15(MULT16_16_16(g1,g1),C[7]));
         sum = SUB32(sum,MULT16_32_Q15(MULT16_16_16(g2,g2),C[8]));

         /* We could force "safe" pitch values to handle packet loss better */

         if (sum>best_sum || i==0)
         {
            best_sum=sum;
            best_cdbk=i;
         }
      }
#ifdef FIXED_POINT
      gain[0] = 32+(spx_word16_t)gain_cdbk[best_cdbk*3];
      gain[1] = 32+(spx_word16_t)gain_cdbk[best_cdbk*3+1];
      gain[2] = 32+(spx_word16_t)gain_cdbk[best_cdbk*3+2];
#else
      gain[0] = 0.015625*gain_cdbk[best_cdbk*3]  + .5;
      gain[1] = 0.015625*gain_cdbk[best_cdbk*3+1]+ .5;
      gain[2] = 0.015625*gain_cdbk[best_cdbk*3+2]+ .5;
#endif
      *cdbk_index=best_cdbk;
   }

#ifdef FIXED_POINT
   for (i=0;i<nsf;i++)
     exc[i]=SHL(MULT16_32_Q15(SHL(gain[0],7),e[2][i])+MULT16_32_Q15(SHL(gain[1],7),e[1][i])+MULT16_32_Q15(SHL(gain[2],7),e[0][i]),2);
   
   err=0;
   for (i=0;i<nsf;i++)
   {
      spx_sig_t tmp = SHL((MULT16_32_Q15(SHL(gain[0],7),x[2][i])+MULT16_32_Q15(SHL(gain[1],7),x[1][i])+MULT16_32_Q15(SHL(gain[2],7),x[0][i])),2);
      spx_sig_t perr=target[i]-tmp;
      new_target[i] = target[i] - tmp;
      spx_word16_t perr2 = PSHR(perr,15);
      err = ADD64(err,MULT16_16(perr2,perr2));
      
   }
#else
   for (i=0;i<nsf;i++)
      exc[i]=gain[0]*e[2][i]+gain[1]*e[1][i]+gain[2]*e[0][i];
   
   err=0;
   for (i=0;i<nsf;i++)
      err+=(target[i]-gain[2]*x[0][i]-gain[1]*x[1][i]-gain[0]*x[2][i])
      * (target[i]-gain[2]*x[0][i]-gain[1]*x[1][i]-gain[0]*x[2][i]);
#endif

   return err;
}


/** Finds the best quantized 3-tap pitch predictor by analysis by synthesis */
int pitch_search_3tap(
spx_sig_t target[],                 /* Target vector */
spx_sig_t *sw,
spx_coef_t ak[],                     /* LPCs for this subframe */
spx_coef_t awk1[],                   /* Weighted LPCs #1 for this subframe */
spx_coef_t awk2[],                   /* Weighted LPCs #2 for this subframe */
spx_sig_t exc[],                    /* Excitation */
const void *par,
int   start,                    /* Smallest pitch value allowed */
int   end,                      /* Largest pitch value allowed */
spx_word16_t pitch_coef,               /* Voicing (pitch) coefficient */
int   p,                        /* Number of LPC coeffs */
int   nsf,                      /* Number of samples in subframe */
SpeexBits *bits,
char *stack,
spx_sig_t *exc2,
spx_sig_t *r,
int complexity,
int cdbk_offset
)
{
   int i,j;
   int cdbk_index, pitch=0, best_gain_index=0;
   spx_sig_t *best_exc, *new_target, *best_target;
   int best_pitch=0;
   spx_word64_t err, best_err=-1;
   int N;
   ltp_params *params;
   int *nbest;
   spx_word16_t *gains;

   N=complexity;
   if (N>10)
      N=10;
   if (N<1)
      N=1;

   nbest=PUSH(stack, N, int);
   gains = PUSH(stack, N, spx_word16_t);
   params = (ltp_params*) par;

   if (end<start)
   {
      speex_bits_pack(bits, 0, params->pitch_bits);
      speex_bits_pack(bits, 0, params->gain_bits);
      for (i=0;i<nsf;i++)
         exc[i]=0;
      return start;
   }
   
   best_exc=PUSH(stack,nsf, spx_sig_t);
   new_target=PUSH(stack,nsf, spx_sig_t);
   best_target=PUSH(stack,nsf, spx_sig_t);
   
   if (N>end-start+1)
      N=end-start+1;
   open_loop_nbest_pitch(sw, start, end, nsf, nbest, gains, N, stack);
   for (i=0;i<N;i++)
   {
      pitch=nbest[i];
      for (j=0;j<nsf;j++)
         exc[j]=0;
      err=pitch_gain_search_3tap(target, ak, awk1, awk2, exc, par, pitch, p, nsf,
                                 bits, stack, exc2, r, new_target, &cdbk_index, cdbk_offset);
      if (err<best_err || best_err<0)
      {
         for (j=0;j<nsf;j++)
            best_exc[j]=exc[j];
         for (j=0;j<nsf;j++)
            best_target[j]=new_target[j];
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
   for (i=0;i<nsf;i++)
      target[i]=best_target[i];

   return pitch;
}

void pitch_unquant_3tap(
spx_sig_t exc[],                    /* Excitation */
int   start,                    /* Smallest pitch value allowed */
int   end,                      /* Largest pitch value allowed */
spx_word16_t pitch_coef,               /* Voicing (pitch) coefficient */
const void *par,
int   nsf,                      /* Number of samples in subframe */
int *pitch_val,
spx_word16_t *gain_val,
SpeexBits *bits,
char *stack,
int count_lost,
int subframe_offset,
spx_word16_t last_pitch_gain,
int cdbk_offset
)
{
   int i;
   int pitch;
   int gain_index;
   spx_word16_t gain[3];
   const signed char *gain_cdbk;
   int gain_cdbk_size;
   ltp_params *params;

   params = (ltp_params*) par;
   gain_cdbk_size = 1<<params->gain_bits;
   gain_cdbk = params->gain_cdbk + 3*gain_cdbk_size*cdbk_offset;

   pitch = speex_bits_unpack_unsigned(bits, params->pitch_bits);
   pitch += start;
   gain_index = speex_bits_unpack_unsigned(bits, params->gain_bits);
   /*printf ("decode pitch: %d %d\n", pitch, gain_index);*/
#ifdef FIXED_POINT
   gain[0] = 32+(spx_word16_t)gain_cdbk[gain_index*3];
   gain[1] = 32+(spx_word16_t)gain_cdbk[gain_index*3+1];
   gain[2] = 32+(spx_word16_t)gain_cdbk[gain_index*3+2];
#else
   gain[0] = 0.015625*gain_cdbk[gain_index*3]+.5;
   gain[1] = 0.015625*gain_cdbk[gain_index*3+1]+.5;
   gain[2] = 0.015625*gain_cdbk[gain_index*3+2]+.5;
#endif

   if (count_lost && pitch > subframe_offset)
   {
      float gain_sum;
      if (1) {
	 float tmp = count_lost < 4 ? GAIN_SCALING_1*last_pitch_gain : 0.4 * GAIN_SCALING_1 * last_pitch_gain;
         if (tmp>.95)
            tmp=.95;
         gain_sum = GAIN_SCALING_1*gain_3tap_to_1tap(gain);

	 if (gain_sum > tmp) {
	    float fact = tmp/gain_sum;
	    for (i=0;i<3;i++)
	       gain[i]*=fact;

	 }

      }

   }

   *pitch_val = pitch;
   gain_val[0]=gain[0];
   gain_val[1]=gain[1];
   gain_val[2]=gain[2];

   {
      spx_sig_t *e[3];
      spx_sig_t *tmp2;
      tmp2=PUSH(stack, 3*nsf, spx_sig_t);
      e[0]=tmp2;
      e[1]=tmp2+nsf;
      e[2]=tmp2+2*nsf;
      
      for (i=0;i<3;i++)
      {
         int j;
         int pp=pitch+1-i;
#if 0
         for (j=0;j<nsf;j++)
         {
            if (j-pp<0)
               e[i][j]=exc[j-pp];
            else if (j-pp-pitch<0)
               e[i][j]=exc[j-pp-pitch];
            else
               e[i][j]=0;
         }
#else
         {
            int tmp1, tmp3;
            tmp1=nsf;
            if (tmp1>pp)
               tmp1=pp;
            for (j=0;j<tmp1;j++)
               e[i][j]=exc[j-pp];
            tmp3=nsf;
            if (tmp3>pp+pitch)
               tmp3=pp+pitch;
            for (j=tmp1;j<tmp3;j++)
               e[i][j]=exc[j-pp-pitch];
            for (j=tmp3;j<nsf;j++)
               e[i][j]=0;
         }
#endif
      }

#ifdef FIXED_POINT
      {
         for (i=0;i<nsf;i++)
            exc[i]=SHL(MULT16_32_Q15(SHL(gain[0],7),e[2][i])+MULT16_32_Q15(SHL(gain[1],7),e[1][i])+MULT16_32_Q15(SHL(gain[2],7),e[0][i]),2);
      }
#else
      for (i=0;i<nsf;i++)
         exc[i]=gain[0]*e[2][i]+gain[1]*e[1][i]+gain[2]*e[0][i];
#endif
   }
}


/** Forced pitch delay and gain */
int forced_pitch_quant(
spx_sig_t target[],                 /* Target vector */
spx_sig_t *sw,
spx_coef_t ak[],                     /* LPCs for this subframe */
spx_coef_t awk1[],                   /* Weighted LPCs #1 for this subframe */
spx_coef_t awk2[],                   /* Weighted LPCs #2 for this subframe */
spx_sig_t exc[],                    /* Excitation */
const void *par,
int   start,                    /* Smallest pitch value allowed */
int   end,                      /* Largest pitch value allowed */
spx_word16_t pitch_coef,               /* Voicing (pitch) coefficient */
int   p,                        /* Number of LPC coeffs */
int   nsf,                      /* Number of samples in subframe */
SpeexBits *bits,
char *stack,
spx_sig_t *exc2,
spx_sig_t *r,
int complexity,
int cdbk_offset
)
{
   int i;
   float coef = GAIN_SCALING_1*pitch_coef;
   if (coef>.99)
      coef=.99;
   for (i=0;i<nsf;i++)
   {
      exc[i]=exc[i-start]*coef;
   }
   return start;
}

/** Unquantize forced pitch delay and gain */
void forced_pitch_unquant(
spx_sig_t exc[],                    /* Excitation */
int   start,                    /* Smallest pitch value allowed */
int   end,                      /* Largest pitch value allowed */
spx_word16_t pitch_coef,               /* Voicing (pitch) coefficient */
const void *par,
int   nsf,                      /* Number of samples in subframe */
int *pitch_val,
spx_word16_t *gain_val,
SpeexBits *bits,
char *stack,
int count_lost,
int subframe_offset,
spx_word16_t last_pitch_gain,
int cdbk_offset
)
{
   int i;
   float coef = GAIN_SCALING_1*pitch_coef;
   if (coef>.99)
      coef=.99;
   for (i=0;i<nsf;i++)
   {
      exc[i]=exc[i-start]*coef;
   }
   *pitch_val = start;
   gain_val[0]=gain_val[2]=0;
   gain_val[1] = pitch_coef;
}
