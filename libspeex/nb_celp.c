/* Copyright (C) 2002 Jean-Marc Valin 
   File: nb_celp.c

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

#include <math.h>
#include "nb_celp.h"
#include "lpc.h"
#include "lsp.h"
#include "ltp.h"
#include "quant_lsp.h"
#include "cb_search.h"
#include "filters.h"
#include "stack_alloc.h"
#include "vq.h"
#include "speex_bits.h"
#include "vbr.h"
#include "misc.h"
#include "speex_callbacks.h"

#include <stdio.h>

#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif

#ifndef NULL
#define NULL 0
#endif

#define SUBMODE(x) st->submodes[st->submodeID]->x


#ifdef FIXED_POINT
spx_word32_t ol_gain_table[32]={18900, 25150, 33468, 44536, 59265, 78865, 104946, 139653, 185838, 247297, 329081, 437913, 582736, 775454, 1031906, 1373169, 1827293, 2431601, 3235761, 4305867, 5729870, 7624808, 10146425, 13501971, 17967238, 23909222, 31816294, 42338330, 56340132, 74972501, 99766822, 132760927};
spx_word16_t exc_gain_quant_scal3_bound[7]={1841, 3883, 6051, 8062, 10444, 13580, 18560};
spx_word16_t exc_gain_quant_scal3[8]={1002, 2680, 5086, 7016, 9108, 11781, 15380, 21740};
spx_word16_t exc_gain_quant_scal1_bound[1]={14385};
spx_word16_t exc_gain_quant_scal1[2]={11546, 17224};

#define LSP_MARGIN 16
#define LSP_DELTA 

#else

float exc_gain_quant_scal3_bound[7]={0.112338, 0.236980, 0.369316, 0.492054, 0.637471, 0.828874, 1.132784};
float exc_gain_quant_scal3[8]={0.061130, 0.163546, 0.310413, 0.428220, 0.555887, 0.719055, 0.938694, 1.326874};
float exc_gain_quant_scal1_bound[1]={0.87798};
float exc_gain_quant_scal1[2]={0.70469, 1.05127};

#define LSP_MARGIN .002
#define LSP_DELTA 

#endif




#define sqr(x) ((x)*(x))

void *nb_encoder_init(SpeexMode *m)
{
   EncState *st;
   SpeexNBMode *mode;
   int i;

   mode=(SpeexNBMode *)m->mode;
   st = (EncState*)speex_alloc(sizeof(EncState)+8000*sizeof(spx_sig_t));
   if (!st)
      return NULL;

   st->stack = ((char*)st) + sizeof(EncState);
   
   st->mode=m;

   st->frameSize = mode->frameSize;
   st->windowSize = st->frameSize*3/2;
   st->nbSubframes=mode->frameSize/mode->subframeSize;
   st->subframeSize=mode->subframeSize;
   st->lpcSize = mode->lpcSize;
   st->bufSize = mode->bufSize;
   st->gamma1=mode->gamma1;
   st->gamma2=mode->gamma2;
   st->min_pitch=mode->pitchStart;
   st->max_pitch=mode->pitchEnd;
   st->lag_factor=mode->lag_factor;
   st->lpc_floor = mode->lpc_floor;
  
   st->submodes=mode->submodes;
   st->submodeID=st->submodeSelect=mode->defaultSubmode;
   st->bounded_pitch = 1;

   st->encode_submode = 1;
#ifdef EPIC_48K
   st->lbr_48k=mode->lbr48k;
#endif

   /* Allocating input buffer */
   st->inBuf = PUSH(st->stack, st->bufSize, spx_sig_t);
   st->frame = st->inBuf + st->bufSize - st->windowSize;
   /* Allocating excitation buffer */
   st->excBuf = PUSH(st->stack, st->bufSize, spx_sig_t);
   st->exc = st->excBuf + st->bufSize - st->windowSize;
   st->swBuf = PUSH(st->stack, st->bufSize, spx_sig_t);
   st->sw = st->swBuf + st->bufSize - st->windowSize;

   st->exc2Buf = PUSH(st->stack, st->bufSize, spx_sig_t);
   st->exc2 = st->exc2Buf + st->bufSize - st->windowSize;

   st->innov = PUSH(st->stack, st->frameSize, spx_sig_t);

   /* Asymmetric "pseudo-Hamming" window */
   {
      int part1, part2;
      part1=st->frameSize - (st->subframeSize>>1);
      part2=(st->frameSize>>1) + (st->subframeSize>>1);
      st->window = PUSH(st->stack, st->windowSize, spx_word16_t);
      for (i=0;i<part1;i++)
         st->window[i]=(spx_word16_t)(SIG_SCALING*(.54-.46*cos(M_PI*i/part1)));
      for (i=0;i<part2;i++)
         st->window[part1+i]=(spx_word16_t)(SIG_SCALING*(.54+.46*cos(M_PI*i/part2)));
   }
   /* Create the window for autocorrelation (lag-windowing) */
   st->lagWindow = PUSH(st->stack, st->lpcSize+1, spx_word16_t);
   for (i=0;i<st->lpcSize+1;i++)
      st->lagWindow[i]=16384*exp(-.5*sqr(2*M_PI*st->lag_factor*i));

   st->autocorr = PUSH(st->stack, st->lpcSize+1, spx_word16_t);

   st->buf2 = PUSH(st->stack, st->windowSize, spx_sig_t);

   st->lpc = PUSH(st->stack, st->lpcSize+1, spx_coef_t);
   st->interp_lpc = PUSH(st->stack, st->lpcSize+1, spx_coef_t);
   st->interp_qlpc = PUSH(st->stack, st->lpcSize+1, spx_coef_t);
   st->bw_lpc1 = PUSH(st->stack, st->lpcSize+1, spx_coef_t);
   st->bw_lpc2 = PUSH(st->stack, st->lpcSize+1, spx_coef_t);

   st->lsp = PUSH(st->stack, st->lpcSize, spx_lsp_t);
   st->qlsp = PUSH(st->stack, st->lpcSize, spx_lsp_t);
   st->old_lsp = PUSH(st->stack, st->lpcSize, spx_lsp_t);
   st->old_qlsp = PUSH(st->stack, st->lpcSize, spx_lsp_t);
   st->interp_lsp = PUSH(st->stack, st->lpcSize, spx_lsp_t);
   st->interp_qlsp = PUSH(st->stack, st->lpcSize, spx_lsp_t);

   st->first = 1;
   for (i=0;i<st->lpcSize;i++)
   {
      st->lsp[i]=LSP_SCALING*(M_PI*((float)(i+1)))/(st->lpcSize+1);
   }

   st->mem_sp = PUSH(st->stack, st->lpcSize, spx_mem_t);
   st->mem_sw = PUSH(st->stack, st->lpcSize, spx_mem_t);
   st->mem_sw_whole = PUSH(st->stack, st->lpcSize, spx_mem_t);
   st->mem_exc = PUSH(st->stack, st->lpcSize, spx_mem_t);

   st->pi_gain = PUSH(st->stack, st->nbSubframes, spx_word32_t);

   st->pitch = PUSH(st->stack, st->nbSubframes, int);

   st->vbr = PUSHS(st->stack, VBRState);
   vbr_init(st->vbr);
   st->vbr_quality = 8;
   st->vbr_enabled = 0;
   st->vad_enabled = 0;
   st->dtx_enabled = 0;
   st->abr_enabled = 0;
   st->abr_drift = 0;

   st->complexity=2;
   st->sampling_rate=8000;
   st->dtx_count=0;

   return st;
}

void nb_encoder_destroy(void *state)
{
   EncState *st=(EncState *)state;
   /* Free all allocated memory */

   vbr_destroy(st->vbr);

   /*Free state memory... should be last*/
   speex_free(st);
}

int nb_encode(void *state, short *in, SpeexBits *bits)
{
   EncState *st;
   int i, sub, roots;
   int ol_pitch;
   float ol_pitch_coef;
   spx_word32_t ol_gain;
   spx_sig_t *res, *target;
   spx_mem_t *mem;
   char *stack;
   spx_sig_t *syn_resp;
   spx_sig_t *orig;
#ifdef EPIC_48K
   int pitch_half[2];
   int ol_pitch_id=0;
#endif

   st=(EncState *)state;
   stack=st->stack;

   /* Copy new data in input buffer */
   speex_move(st->inBuf, st->inBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(spx_sig_t));
   for (i=0;i<st->frameSize;i++)
      st->inBuf[st->bufSize-st->frameSize+i] = SHL((int)in[i], SIG_SHIFT);

   /* Move signals 1 frame towards the past */
   speex_move(st->exc2Buf, st->exc2Buf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(spx_sig_t));
   speex_move(st->excBuf, st->excBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(spx_sig_t));
   speex_move(st->swBuf, st->swBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(spx_sig_t));


   {
      spx_word16_t *w_sig;
      w_sig = PUSH(stack, st->windowSize, spx_word16_t);
      /* Window for analysis */
      for (i=0;i<st->windowSize;i++)
         w_sig[i] = SHR(MULT16_16(SHR((spx_word32_t)(st->frame[i]),SIG_SHIFT),st->window[i]),SIG_SHIFT);

      /* Compute auto-correlation */
      _spx_autocorr(w_sig, st->autocorr, st->lpcSize+1, st->windowSize);
   }
   st->autocorr[0] = (spx_word16_t) (st->autocorr[0]*st->lpc_floor); /* Noise floor in auto-correlation domain */

   /* Lag windowing: equivalent to filtering in the power-spectrum domain */
   for (i=0;i<st->lpcSize+1;i++)
      st->autocorr[i] = MULT16_16_Q14(st->autocorr[i],st->lagWindow[i]);

   /* Levinson-Durbin */
   _spx_lpc(st->lpc+1, st->autocorr, st->lpcSize);
   st->lpc[0]=(spx_coef_t)LPC_SCALING;

   /* LPC to LSPs (x-domain) transform */
   roots=lpc_to_lsp (st->lpc, st->lpcSize, st->lsp, 15, 0.2, stack);
   /* Check if we found all the roots */
   if (roots!=st->lpcSize)
   {
      /* Search again if we can afford it */
      if (st->complexity>1)
         roots = lpc_to_lsp (st->lpc, st->lpcSize, st->lsp, 11, 0.05, stack);
      if (roots!=st->lpcSize) 
      {
         /*If we can't find all LSP's, do some damage control and use previous filter*/
         for (i=0;i<st->lpcSize;i++)
         {
            st->lsp[i]=st->old_lsp[i];
         }
      }
   }



   /* Whole frame analysis (open-loop estimation of pitch and excitation gain) */
   {
      if (st->first)
         for (i=0;i<st->lpcSize;i++)
            st->interp_lsp[i] = st->lsp[i];
      else
         lsp_interpolate(st->old_lsp, st->lsp, st->interp_lsp, st->lpcSize, st->nbSubframes, st->nbSubframes<<1);

      lsp_enforce_margin(st->interp_lsp, st->lpcSize, LSP_MARGIN);

      /* Compute interpolated LPCs (unquantized) for whole frame*/
      lsp_to_lpc(st->interp_lsp, st->interp_lpc, st->lpcSize,stack);


      /*Open-loop pitch*/
      if (!st->submodes[st->submodeID] || st->vbr_enabled || st->vad_enabled || SUBMODE(forced_pitch_gain) ||
          SUBMODE(lbr_pitch) != -1)
      {
         int nol_pitch[6];
         spx_word16_t nol_pitch_coef[6];
         
         bw_lpc(st->gamma1, st->interp_lpc, st->bw_lpc1, st->lpcSize);
         bw_lpc(st->gamma2, st->interp_lpc, st->bw_lpc2, st->lpcSize);
         
         filter_mem2(st->frame, st->bw_lpc1, st->bw_lpc2, st->sw, st->frameSize, st->lpcSize, st->mem_sw_whole);

         open_loop_nbest_pitch(st->sw, st->min_pitch, st->max_pitch, st->frameSize, 
                               nol_pitch, nol_pitch_coef, 6, stack);
         ol_pitch=nol_pitch[0];
         ol_pitch_coef = nol_pitch_coef[0];
         /*Try to remove pitch multiples*/
         for (i=1;i<6;i++)
         {
#ifdef FIXED_POINT
            if ((nol_pitch_coef[i]>MULT16_16_Q15(nol_pitch_coef[0],27853)) && 
#else
            if ((nol_pitch_coef[i]>.85*nol_pitch_coef[0]) && 
#endif
                (ABS(2*nol_pitch[i]-ol_pitch)<=2 || ABS(3*nol_pitch[i]-ol_pitch)<=3 || 
                 ABS(4*nol_pitch[i]-ol_pitch)<=4 || ABS(5*nol_pitch[i]-ol_pitch)<=5))
            {
               /*ol_pitch_coef=nol_pitch_coef[i];*/
               ol_pitch = nol_pitch[i];
            }
         }
         /*if (ol_pitch>50)
           ol_pitch/=2;*/
         /*ol_pitch_coef = sqrt(ol_pitch_coef);*/

#ifdef EPIC_48K
         if (st->lbr_48k)
         {
            if (ol_pitch < st->min_pitch+2)
               ol_pitch = st->min_pitch+2;
            if (ol_pitch > st->max_pitch-2)
               ol_pitch = st->max_pitch-2;
            open_loop_nbest_pitch(st->sw, ol_pitch-2, ol_pitch+2, st->frameSize>>1, 
                                  &pitch_half[0], nol_pitch_coef, 1, stack);
            open_loop_nbest_pitch(st->sw+(st->frameSize>>1), pitch_half[0]-1, pitch_half[0]+2, st->frameSize>>1, 
                                  &pitch_half[1], nol_pitch_coef, 1, stack);
         }
#endif

#ifdef FIXED_POINT
         ol_pitch_coef *= 0.0039062;
#endif
      } else {
         ol_pitch=0;
         ol_pitch_coef=0;
      }
      /*Compute "real" excitation*/
      fir_mem2(st->frame, st->interp_lpc, st->exc, st->frameSize, st->lpcSize, st->mem_exc);

      /* Compute open-loop excitation gain */
#ifdef EPIC_48K
      if (st->lbr_48k)
      {
         float ol1=0,ol2=0;
         float ol_gain2;
         ol1 = compute_rms(st->exc, st->frameSize>>1);
         ol2 = compute_rms(st->exc+(st->frameSize>>1), st->frameSize>>1);
         ol1 *= ol1*(st->frameSize>>1);
         ol2 *= ol2*(st->frameSize>>1);

         ol_gain2=ol1;
         if (ol2>ol1)
            ol_gain2=ol2;
         ol_gain2 = sqrt(2*ol_gain2*(ol1+ol2))*1.3*(1-.5*ol_pitch_coef*ol_pitch_coef);
      
         ol_gain=SHR(sqrt(1+ol_gain2/st->frameSize),SIG_SHIFT);

      } else {
#endif
         ol_gain = SHL(compute_rms(st->exc, st->frameSize),SIG_SHIFT);
#ifdef EPIC_48K
      }
#endif
   }

   /*VBR stuff*/
   if (st->vbr && (st->vbr_enabled||st->vad_enabled))
   {
      float lsp_dist=0;
      for (i=0;i<st->lpcSize;i++)
         lsp_dist += (st->old_lsp[i] - st->lsp[i])*(st->old_lsp[i] - st->lsp[i]);
      lsp_dist /= LSP_SCALING*LSP_SCALING;
      
      if (st->abr_enabled)
      {
         float qual_change=0;
         if (st->abr_drift2 * st->abr_drift > 0)
         {
            /* Only adapt if long-term and short-term drift are the same sign */
            qual_change = -.00001*st->abr_drift/(1+st->abr_count);
            if (qual_change>.05)
               qual_change=.05;
            if (qual_change<-.05)
               qual_change=-.05;
         }
         st->vbr_quality += qual_change;
         if (st->vbr_quality>10)
            st->vbr_quality=10;
         if (st->vbr_quality<0)
            st->vbr_quality=0;
      }

      st->relative_quality = vbr_analysis(st->vbr, in, st->frameSize, ol_pitch, ol_pitch_coef);
      /*if (delta_qual<0)*/
      /*  delta_qual*=.1*(3+st->vbr_quality);*/
      if (st->vbr_enabled) 
      {
         int mode;
         int choice=0;
         float min_diff=100;
         mode = 8;
         while (mode)
         {
            int v1;
            float thresh;
            v1=(int)floor(st->vbr_quality);
            if (v1==10)
               thresh = vbr_nb_thresh[mode][v1];
            else
               thresh = (st->vbr_quality-v1)*vbr_nb_thresh[mode][v1+1] + (1+v1-st->vbr_quality)*vbr_nb_thresh[mode][v1];
            if (st->relative_quality > thresh && 
                st->relative_quality-thresh<min_diff)
            {
               choice = mode;
               min_diff = st->relative_quality-thresh;
            }
            mode--;
         }
         mode=choice;
         if (mode==0)
         {
            if (st->dtx_count==0 || lsp_dist>.05 || !st->dtx_enabled || st->dtx_count>20)
            {
               mode=1;
               st->dtx_count=1;
            } else {
               mode=0;
               st->dtx_count++;
            }
         } else {
            st->dtx_count=0;
         }

         speex_encoder_ctl(state, SPEEX_SET_MODE, &mode);

         if (st->abr_enabled)
         {
            int bitrate;
            speex_encoder_ctl(state, SPEEX_GET_BITRATE, &bitrate);
            st->abr_drift+=(bitrate-st->abr_enabled);
            st->abr_drift2 = .95*st->abr_drift2 + .05*(bitrate-st->abr_enabled);
            st->abr_count += 1.0;
         }

      } else {
         /*VAD only case*/
         int mode;
         if (st->relative_quality<2)
         {
            if (st->dtx_count==0 || lsp_dist>.05 || !st->dtx_enabled || st->dtx_count>20)
            {
               st->dtx_count=1;
               mode=1;
            } else {
               mode=0;
               st->dtx_count++;
            }
         } else {
            st->dtx_count = 0;
            mode=st->submodeSelect;
         }
         /*speex_encoder_ctl(state, SPEEX_SET_MODE, &mode);*/
         st->submodeID=mode;
      } 
   } else {
      st->relative_quality = -1;
   }

   if (st->encode_submode)
   {
#ifdef EPIC_48K
   if (!st->lbr_48k) {
#endif

   /* First, transmit a zero for narrowband */
   speex_bits_pack(bits, 0, 1);

   /* Transmit the sub-mode we use for this frame */
   speex_bits_pack(bits, st->submodeID, NB_SUBMODE_BITS);

#ifdef EPIC_48K
   }
#endif
   }

   /* If null mode (no transmission), just set a couple things to zero*/
   if (st->submodes[st->submodeID] == NULL)
   {
      for (i=0;i<st->frameSize;i++)
         st->exc[i]=st->exc2[i]=st->sw[i]=VERY_SMALL;

      for (i=0;i<st->lpcSize;i++)
         st->mem_sw[i]=0;
      st->first=1;
      st->bounded_pitch = 1;

      /* Final signal synthesis from excitation */
      iir_mem2(st->exc, st->interp_qlpc, st->frame, st->frameSize, st->lpcSize, st->mem_sp);

      for (i=0;i<st->frameSize;i++)
         in[i]=st->frame[i];

      return 0;

   }

   /* LSP Quantization */
   if (st->first)
   {
      for (i=0;i<st->lpcSize;i++)
         st->old_lsp[i] = st->lsp[i];
   }


   /*Quantize LSPs*/
#if 1 /*0 for unquantized*/
   SUBMODE(lsp_quant)(st->lsp, st->qlsp, st->lpcSize, bits);
#else
   for (i=0;i<st->lpcSize;i++)
     st->qlsp[i]=st->lsp[i];
#endif

#ifdef EPIC_48K
   if (st->lbr_48k) {
      speex_bits_pack(bits, pitch_half[0]-st->min_pitch, 7);
      speex_bits_pack(bits, pitch_half[1]-pitch_half[0]+1, 2);
      
      {
         int quant = (int)floor(.5+7.4*ol_pitch_coef);
         if (quant>7)
            quant=7;
         if (quant<0)
            quant=0;
         ol_pitch_id=quant;
         speex_bits_pack(bits, quant, 3);
         ol_pitch_coef=0.13514*quant;
         
      }
      {
         int qe = (int)(floor(.5+2.1*log(ol_gain*1.0/SIG_SCALING)))-2;
         if (qe<0)
            qe=0;
         if (qe>15)
            qe=15;
         ol_gain = exp((qe+2)/2.1)*SIG_SCALING;
         speex_bits_pack(bits, qe, 4);
      }

   } else {
#endif

   /*If we use low bit-rate pitch mode, transmit open-loop pitch*/
   if (SUBMODE(lbr_pitch)!=-1)
   {
      speex_bits_pack(bits, ol_pitch-st->min_pitch, 7);
   } 

   if (SUBMODE(forced_pitch_gain))
   {
      int quant;
      quant = (int)floor(.5+15*ol_pitch_coef);
      if (quant>15)
         quant=15;
      if (quant<0)
         quant=0;
      speex_bits_pack(bits, quant, 4);
      ol_pitch_coef=0.066667*quant;
   }
   
   
   /*Quantize and transmit open-loop excitation gain*/
#ifdef FIXED_POINT
   {
      int qe = scal_quant32(ol_gain, ol_gain_table, 32);
      /*ol_gain = exp(qe/3.5)*SIG_SCALING;*/
      ol_gain = MULT16_32_Q15(28406,ol_gain_table[qe]);
      speex_bits_pack(bits, qe, 5);
   }
#else
   {
      int qe = (int)(floor(.5+3.5*log(ol_gain*1.0/SIG_SCALING)));
      if (qe<0)
         qe=0;
      if (qe>31)
         qe=31;
      ol_gain = exp(qe/3.5)*SIG_SCALING;
      speex_bits_pack(bits, qe, 5);
   }
#endif


#ifdef EPIC_48K
   }
#endif


   /* Special case for first frame */
   if (st->first)
   {
      for (i=0;i<st->lpcSize;i++)
         st->old_qlsp[i] = st->qlsp[i];
   }

   /* Filter response */
   res = PUSH(stack, st->subframeSize, spx_sig_t);
   /* Target signal */
   target = PUSH(stack, st->subframeSize, spx_sig_t);
   syn_resp = PUSH(stack, st->subframeSize, spx_sig_t);
   mem = PUSH(stack, st->lpcSize, spx_mem_t);
   orig = PUSH(stack, st->frameSize, spx_sig_t);
   for (i=0;i<st->frameSize;i++)
      orig[i]=st->frame[i];

   /* Loop on sub-frames */
   for (sub=0;sub<st->nbSubframes;sub++)
   {
      int   offset;
      spx_sig_t *sp, *sw, *exc, *exc2;
      int pitch;

#ifdef EPIC_48K
      if (st->lbr_48k)
      {
         if (sub*2 < st->nbSubframes)
            ol_pitch = pitch_half[0];
         else
            ol_pitch = pitch_half[1];
      }
#endif

      /* Offset relative to start of frame */
      offset = st->subframeSize*sub;
      /* Original signal */
      sp=st->frame+offset;
      /* Excitation */
      exc=st->exc+offset;
      /* Weighted signal */
      sw=st->sw+offset;

      exc2=st->exc2+offset;


      /* LSP interpolation (quantized and unquantized) */
      lsp_interpolate(st->old_lsp, st->lsp, st->interp_lsp, st->lpcSize, sub, st->nbSubframes);
      lsp_interpolate(st->old_qlsp, st->qlsp, st->interp_qlsp, st->lpcSize, sub, st->nbSubframes);

      /* Make sure the filters are stable */
      lsp_enforce_margin(st->interp_lsp, st->lpcSize, LSP_MARGIN);
      lsp_enforce_margin(st->interp_qlsp, st->lpcSize, LSP_MARGIN);

      /* Compute interpolated LPCs (quantized and unquantized) */
      lsp_to_lpc(st->interp_lsp, st->interp_lpc, st->lpcSize,stack);

      lsp_to_lpc(st->interp_qlsp, st->interp_qlpc, st->lpcSize, stack);

      /* Compute analysis filter gain at w=pi (for use in SB-CELP) */
      {
         spx_word32_t pi_g=st->interp_qlpc[0];
         for (i=1;i<=st->lpcSize;i+=2)
         {
            pi_g += -st->interp_qlpc[i] +  st->interp_qlpc[i+1];
         }
         st->pi_gain[sub] = pi_g;
      }


      /* Compute bandwidth-expanded (unquantized) LPCs for perceptual weighting */
      bw_lpc(st->gamma1, st->interp_lpc, st->bw_lpc1, st->lpcSize);
      if (st->gamma2>=0)
         bw_lpc(st->gamma2, st->interp_lpc, st->bw_lpc2, st->lpcSize);
      else
      {
         st->bw_lpc2[0]=1;
         for (i=1;i<=st->lpcSize;i++)
            st->bw_lpc2[i]=0;
      }

      /* Compute impulse response of A(z/g1) / ( A(z)*A(z/g2) )*/
      for (i=0;i<st->subframeSize;i++)
         exc[i]=VERY_SMALL;
      exc[0]=SIG_SCALING;
      syn_percep_zero(exc, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2, syn_resp, st->subframeSize, st->lpcSize, stack);

      /* Reset excitation */
      for (i=0;i<st->subframeSize;i++)
         exc[i]=VERY_SMALL;
      for (i=0;i<st->subframeSize;i++)
         exc2[i]=VERY_SMALL;

      /* Compute zero response of A(z/g1) / ( A(z/g2) * A(z) ) */
      for (i=0;i<st->lpcSize;i++)
         mem[i]=st->mem_sp[i];
      iir_mem2(exc, st->interp_qlpc, exc, st->subframeSize, st->lpcSize, mem);
      
      for (i=0;i<st->lpcSize;i++)
         mem[i]=st->mem_sw[i];
      filter_mem2(exc, st->bw_lpc1, st->bw_lpc2, res, st->subframeSize, st->lpcSize, mem);
      
      /* Compute weighted signal */
      for (i=0;i<st->lpcSize;i++)
         mem[i]=st->mem_sw[i];
      filter_mem2(sp, st->bw_lpc1, st->bw_lpc2, sw, st->subframeSize, st->lpcSize, mem);
      
      /* Compute target signal */
      for (i=0;i<st->subframeSize;i++)
         target[i]=sw[i]-res[i];

      for (i=0;i<st->subframeSize;i++)
         exc[i]=exc2[i]=0;

      /* If we have a long-term predictor (otherwise, something's wrong) */
      if (SUBMODE(ltp_quant))
      {
         int pit_min, pit_max;
         /* Long-term prediction */
         if (SUBMODE(lbr_pitch) != -1)
         {
            /* Low bit-rate pitch handling */
            int margin;
            margin = SUBMODE(lbr_pitch);
            if (margin)
            {
               if (ol_pitch < st->min_pitch+margin-1)
                  ol_pitch=st->min_pitch+margin-1;
               if (ol_pitch > st->max_pitch-margin)
                  ol_pitch=st->max_pitch-margin;
               pit_min = ol_pitch-margin+1;
               pit_max = ol_pitch+margin;
            } else {
               pit_min=pit_max=ol_pitch;
            }
         } else {
            pit_min = st->min_pitch;
            pit_max = st->max_pitch;
         }
         
         /* Force pitch to use only the current frame if needed */
         if (st->bounded_pitch && pit_max>offset)
            pit_max=offset;

#ifdef EPIC_48K
         if (st->lbr_48k)
         {
            pitch = SUBMODE(ltp_quant)(target, sw, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2,
                                       exc, SUBMODE(ltp_params), pit_min, pit_max, ol_pitch_coef,
                                       st->lpcSize, st->subframeSize, bits, stack, 
                                       exc2, syn_resp, st->complexity, ol_pitch_id);
         } else {
#endif

         /* Perform pitch search */
         pitch = SUBMODE(ltp_quant)(target, sw, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2,
                                    exc, SUBMODE(ltp_params), pit_min, pit_max, ol_pitch_coef,
                                    st->lpcSize, st->subframeSize, bits, stack, 
                                    exc2, syn_resp, st->complexity, 0);
#ifdef EPIC_48K
         }
#endif

         st->pitch[sub]=pitch;
      } else {
         speex_error ("No pitch prediction, what's wrong");
      }

      /* Update target for adaptive codebook contribution */
      syn_percep_zero(exc, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2, res, st->subframeSize, st->lpcSize, stack);
      for (i=0;i<st->subframeSize;i++)
         target[i]-=res[i];


      /* Quantization of innovation */
      {
         spx_sig_t *innov;
         spx_word32_t ener=0;
         spx_word16_t fine_gain;

         innov = st->innov+sub*st->subframeSize;
         for (i=0;i<st->subframeSize;i++)
            innov[i]=0;
         
         residue_percep_zero(target, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2, st->buf2, st->subframeSize, st->lpcSize, stack);

         ener = SHL(compute_rms(st->buf2, st->subframeSize),SIG_SHIFT);

         /*for (i=0;i<st->subframeSize;i++)
            printf ("%f\n", st->buf2[i]/ener);
         */
         
         fine_gain = DIV32_16(ener,SHR(ol_gain,SIG_SHIFT));

         /* Calculate gain correction for the sub-frame (if any) */
         if (SUBMODE(have_subframe_gain)) 
         {
            int qe;
            if (SUBMODE(have_subframe_gain)==3)
            {
               qe = scal_quant(fine_gain, exc_gain_quant_scal3_bound, 8);
               speex_bits_pack(bits, qe, 3);
               ener=MULT16_32_Q14(exc_gain_quant_scal3[qe],ol_gain);
            } else {
               qe = scal_quant(fine_gain, exc_gain_quant_scal1_bound, 2);
               speex_bits_pack(bits, qe, 1);
               ener=MULT16_32_Q14(exc_gain_quant_scal1[qe],ol_gain);               
            }
         } else {
            ener=ol_gain;
         }

         /*printf ("%f %f\n", ener, ol_gain);*/

         /* Normalize innovation */
         signal_div(target, target, ener, st->subframeSize);

         /* Quantize innovation */
         if (SUBMODE(innovation_quant))
         {
            /* Codebook search */
            SUBMODE(innovation_quant)(target, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2, 
                                      SUBMODE(innovation_params), st->lpcSize, st->subframeSize, 
                                      innov, syn_resp, bits, stack, st->complexity);
            
            /* De-normalize innovation and update excitation */
            signal_mul(innov, innov, ener, st->subframeSize);

            for (i=0;i<st->subframeSize;i++)
               exc[i] += innov[i];
         } else {
            speex_error("No fixed codebook");
         }

         /* In some (rare) modes, we do a second search (more bits) to reduce noise even more */
         if (SUBMODE(double_codebook)) {
            char *tmp_stack=stack;
            spx_sig_t *innov2 = PUSH(tmp_stack, st->subframeSize, spx_sig_t);
            for (i=0;i<st->subframeSize;i++)
               innov2[i]=0;
            for (i=0;i<st->subframeSize;i++)
               target[i]*=2.2;
            SUBMODE(innovation_quant)(target, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2, 
                                      SUBMODE(innovation_params), st->lpcSize, st->subframeSize, 
                                      innov2, syn_resp, bits, tmp_stack, st->complexity);
            signal_mul(innov2, innov2, ener*(1/2.2), st->subframeSize);
            for (i=0;i<st->subframeSize;i++)
               exc[i] += innov2[i];
         }

         signal_mul(target, target, ener, st->subframeSize);
      }

      /*Keep the previous memory*/
      for (i=0;i<st->lpcSize;i++)
         mem[i]=st->mem_sp[i];
      /* Final signal synthesis from excitation */
      iir_mem2(exc, st->interp_qlpc, sp, st->subframeSize, st->lpcSize, st->mem_sp);

      /* Compute weighted signal again, from synthesized speech (not sure it's the right thing) */
      filter_mem2(sp, st->bw_lpc1, st->bw_lpc2, sw, st->subframeSize, st->lpcSize, st->mem_sw);
      for (i=0;i<st->subframeSize;i++)
         exc2[i]=exc[i];
   }

   /* Store the LSPs for interpolation in the next frame */
   if (st->submodeID>=1)
   {
      for (i=0;i<st->lpcSize;i++)
         st->old_lsp[i] = st->lsp[i];
      for (i=0;i<st->lpcSize;i++)
         st->old_qlsp[i] = st->qlsp[i];
   }

   if (st->submodeID==1)
   {
      if (st->dtx_count)
         speex_bits_pack(bits, 15, 4);
      else
         speex_bits_pack(bits, 0, 4);
   }

   /* The next frame will not be the first (Duh!) */
   st->first = 0;

   if (0) {
      float ener=0, err=0;
      float snr;
      for (i=0;i<st->frameSize;i++)
      {
         ener+=st->frame[i]*st->frame[i];
         err += (st->frame[i]-orig[i])*(st->frame[i]-orig[i]);
      }
      snr = 10*log10((ener+1)/(err+1));
      /*printf ("%f %f %f\n", snr, ener, err);*/
   }

   /* Replace input by synthesized speech */
   for (i=0;i<st->frameSize;i++)
   {
      spx_word32_t sig = PSHR(st->frame[i],SIG_SHIFT);
      if (sig>32767)
         sig = 32767;
      if (sig<-32767)
         sig = -32767;
     in[i]=sig;
   }

   if (SUBMODE(innovation_quant) == noise_codebook_quant || st->submodeID==0)
      st->bounded_pitch = 1;
   else
      st->bounded_pitch = 0;

   return 1;
}


void *nb_decoder_init(SpeexMode *m)
{
   DecState *st;
   SpeexNBMode *mode;
   int i;

   mode=(SpeexNBMode*)m->mode;
   st = (DecState *)speex_alloc(sizeof(DecState)+4000*sizeof(spx_sig_t));
   st->mode=m;

   st->stack = ((char*)st) + sizeof(DecState);

   st->encode_submode = 1;
#ifdef EPIC_48K
   st->lbr_48k=mode->lbr48k;
#endif

   st->first=1;
   /* Codec parameters, should eventually have several "modes"*/
   st->frameSize = mode->frameSize;
   st->windowSize = st->frameSize*3/2;
   st->nbSubframes=mode->frameSize/mode->subframeSize;
   st->subframeSize=mode->subframeSize;
   st->lpcSize = mode->lpcSize;
   st->bufSize = mode->bufSize;
   st->gamma1=mode->gamma1;
   st->gamma2=mode->gamma2;
   st->min_pitch=mode->pitchStart;
   st->max_pitch=mode->pitchEnd;

   st->submodes=mode->submodes;
   st->submodeID=mode->defaultSubmode;

   st->lpc_enh_enabled=0;


   st->inBuf = PUSH(st->stack, st->bufSize, spx_sig_t);
   st->frame = st->inBuf + st->bufSize - st->windowSize;
   st->excBuf = PUSH(st->stack, st->bufSize, spx_sig_t);
   st->exc = st->excBuf + st->bufSize - st->windowSize;
   for (i=0;i<st->bufSize;i++)
      st->inBuf[i]=0;
   for (i=0;i<st->bufSize;i++)
      st->excBuf[i]=0;
   st->innov = PUSH(st->stack, st->frameSize, spx_sig_t);

   st->interp_qlpc = PUSH(st->stack, st->lpcSize+1, spx_coef_t);
   st->qlsp = PUSH(st->stack, st->lpcSize, spx_lsp_t);
   st->old_qlsp = PUSH(st->stack, st->lpcSize, spx_lsp_t);
   st->interp_qlsp = PUSH(st->stack, st->lpcSize, spx_lsp_t);
   st->mem_sp = PUSH(st->stack, 5*st->lpcSize, spx_mem_t);
   st->comb_mem = PUSHS(st->stack, CombFilterMem);
   comb_filter_mem_init (st->comb_mem);

   st->pi_gain = PUSH(st->stack, st->nbSubframes, spx_word32_t);
   st->last_pitch = 40;
   st->count_lost=0;
   st->pitch_gain_buf[0] = st->pitch_gain_buf[1] = st->pitch_gain_buf[2] = 0;
   st->pitch_gain_buf_idx = 0;

   st->sampling_rate=8000;
   st->last_ol_gain = 0;

   st->user_callback.func = &speex_default_user_handler;
   st->user_callback.data = NULL;
   for (i=0;i<16;i++)
      st->speex_callbacks[i].func = NULL;

   st->voc_m1=st->voc_m2=st->voc_mean=0;
   st->voc_offset=0;
   st->dtx_enabled=0;
   return st;
}

void nb_decoder_destroy(void *state)
{
   DecState *st;
   st=(DecState*)state;
   
   speex_free(state);
}

#define median3(a, b, c)	((a) < (b) ? ((b) < (c) ? (b) : ((a) < (c) ? (c) : (a))) : ((c) < (b) ? (b) : ((c) < (a) ? (c) : (a))))

static void nb_decode_lost(DecState *st, short *out, char *stack)
{
   int i, sub;
   spx_coef_t *awk1, *awk2, *awk3;
   float pitch_gain, fact, gain_med;

   fact = exp(-.04*st->count_lost*st->count_lost);
   gain_med = median3(st->pitch_gain_buf[0], st->pitch_gain_buf[1], st->pitch_gain_buf[2]);
   if (gain_med < st->last_pitch_gain)
      st->last_pitch_gain = gain_med;
   
   pitch_gain = st->last_pitch_gain;
   if (pitch_gain>.95)
      pitch_gain=.95;

   pitch_gain *= fact;

   /* Shift all buffers by one frame */
   speex_move(st->inBuf, st->inBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(spx_sig_t));
   speex_move(st->excBuf, st->excBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(spx_sig_t));

   awk1=PUSH(stack, (st->lpcSize+1), spx_coef_t);
   awk2=PUSH(stack, (st->lpcSize+1), spx_coef_t);
   awk3=PUSH(stack, (st->lpcSize+1), spx_coef_t);

   for (sub=0;sub<st->nbSubframes;sub++)
   {
      int offset;
      spx_sig_t *sp, *exc;
      /* Offset relative to start of frame */
      offset = st->subframeSize*sub;
      /* Original signal */
      sp=st->frame+offset;
      /* Excitation */
      exc=st->exc+offset;
      /* Excitation after post-filter*/

      /* Calculate perceptually enhanced LPC filter */
      if (st->lpc_enh_enabled)
      {
         float r=.9;
         
         float k1,k2,k3;
         if (st->submodes[st->submodeID] != NULL)
         {
            k1=SUBMODE(lpc_enh_k1);
            k2=SUBMODE(lpc_enh_k2);
         } else {
            k1=k2=.7;
         }
         k3=(1-(1-r*k1)/(1-r*k2))/r;
         if (!st->lpc_enh_enabled)
         {
            k1=k2;
            k3=0;
         }
         bw_lpc(k1, st->interp_qlpc, awk1, st->lpcSize);
         bw_lpc(k2, st->interp_qlpc, awk2, st->lpcSize);
         bw_lpc(k3, st->interp_qlpc, awk3, st->lpcSize);
      }
        
      /* Make up a plausible excitation */
      /* THIS CAN BE IMPROVED */
      /*if (pitch_gain>.95)
        pitch_gain=.95;*/
      {
         float innov_gain=0;
         for (i=0;i<st->frameSize;i++)
            innov_gain += st->innov[i]*st->innov[i];
         innov_gain=sqrt(innov_gain/st->frameSize);
      for (i=0;i<st->subframeSize;i++)
      {
#if 0
         exc[i] = pitch_gain * exc[i - st->last_pitch] + fact*sqrt(1-pitch_gain)*st->innov[i+offset];
         /*Just so it give the same lost packets as with if 0*/
         /*rand();*/
#else
         /*exc[i]=pitch_gain*exc[i-st->last_pitch] +  fact*st->innov[i+offset];*/
         exc[i]=pitch_gain*exc[i-st->last_pitch] + 
         fact*sqrt(1-pitch_gain)*speex_rand(innov_gain);
#endif
      }
      }
      for (i=0;i<st->subframeSize;i++)
         sp[i]=exc[i];
      
      /* Signal synthesis */
      if (st->lpc_enh_enabled)
      {
         filter_mem2(sp, awk2, awk1, sp, st->subframeSize, st->lpcSize, 
                     st->mem_sp+st->lpcSize);
         filter_mem2(sp, awk3, st->interp_qlpc, sp, st->subframeSize, st->lpcSize, 
                     st->mem_sp);
      } else {
         for (i=0;i<st->lpcSize;i++)
            st->mem_sp[st->lpcSize+i] = 0;
         iir_mem2(sp, st->interp_qlpc, sp, st->subframeSize, st->lpcSize, 
                     st->mem_sp);
      }      
   }

   for (i=0;i<st->frameSize;i++)
   {
      spx_word32_t sig = PSHR(st->frame[i],SIG_SHIFT);
      if (sig>32767)
         sig = 32767;
      if (sig<-32767)
         sig = -32767;
     out[i]=sig;
   }
   
   st->first = 0;
   st->count_lost++;
   st->pitch_gain_buf[st->pitch_gain_buf_idx++] = pitch_gain;
   if (st->pitch_gain_buf_idx > 2) /* rollover */
      st->pitch_gain_buf_idx = 0;
}

int nb_decode(void *state, SpeexBits *bits, short *out)
{
   DecState *st;
   int i, sub;
   int pitch;
   float pitch_gain[3];
   spx_word32_t ol_gain=0;
   int ol_pitch=0;
   float ol_pitch_coef=0;
   int best_pitch=40;
   float best_pitch_gain=0;
   int wideband;
   int m;
   char *stack;
   spx_coef_t *awk1, *awk2, *awk3;
   float pitch_average=0;
#ifdef EPIC_48K
   int pitch_half[2];
   int ol_pitch_id=0;
#endif

   st=(DecState*)state;
   stack=st->stack;

   if (st->encode_submode)
   {
#ifdef EPIC_48K
   if (!st->lbr_48k) {
#endif

   /* Check if we're in DTX mode*/
   if (!bits && st->dtx_enabled)
   {
      st->submodeID=0;
   } else 
   {
      /* If bits is NULL, consider the packet to be lost (what could we do anyway) */
      if (!bits)
      {
         nb_decode_lost(st, out, stack);
         return 0;
      }

      /* Search for next narrowband block (handle requests, skip wideband blocks) */
      do {
         if (speex_bits_remaining(bits)<5)
            return -1;
         wideband = speex_bits_unpack_unsigned(bits, 1);
         if (wideband) /* Skip wideband block (for compatibility) */
         {
            int submode;
            int advance;
            advance = submode = speex_bits_unpack_unsigned(bits, SB_SUBMODE_BITS);
            speex_mode_query(&speex_wb_mode, SPEEX_SUBMODE_BITS_PER_FRAME, &advance);
            if (advance < 0)
            {
               speex_warning ("Invalid wideband mode encountered. Corrupted stream?");
               return -2;
            } 
            advance -= (SB_SUBMODE_BITS+1);
            speex_bits_advance(bits, advance);
            
            if (speex_bits_remaining(bits)<5)
               return -1;
            wideband = speex_bits_unpack_unsigned(bits, 1);
            if (wideband)
            {
               advance = submode = speex_bits_unpack_unsigned(bits, SB_SUBMODE_BITS);
               speex_mode_query(&speex_wb_mode, SPEEX_SUBMODE_BITS_PER_FRAME, &advance);
               if (advance < 0)
               {
                  speex_warning ("Invalid wideband mode encountered: corrupted stream?");
                  return -2;
               } 
               advance -= (SB_SUBMODE_BITS+1);
               speex_bits_advance(bits, advance);
               wideband = speex_bits_unpack_unsigned(bits, 1);
               if (wideband)
               {
                  speex_warning ("More than two wideband layers found: corrupted stream?");
                  return -2;
               }

            }
         }
         if (speex_bits_remaining(bits)<4)
            return -1;
         /* FIXME: Check for overflow */
         m = speex_bits_unpack_unsigned(bits, 4);
         if (m==15) /* We found a terminator */
         {
            return -1;
         } else if (m==14) /* Speex in-band request */
         {
            int ret = speex_inband_handler(bits, st->speex_callbacks, state);
            if (ret)
               return ret;
         } else if (m==13) /* User in-band request */
         {
            int ret = st->user_callback.func(bits, state, st->user_callback.data);
            if (ret)
               return ret;
         } else if (m>8) /* Invalid mode */
         {
            speex_warning("Invalid mode encountered: corrupted stream?");
            return -2;
         }
      
      } while (m>8);

      /* Get the sub-mode that was used */
      st->submodeID = m;

   }
#ifdef EPIC_48K
   }
#endif
   }

   /* Shift all buffers by one frame */
   speex_move(st->inBuf, st->inBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(spx_sig_t));
   speex_move(st->excBuf, st->excBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(spx_sig_t));

   /* If null mode (no transmission), just set a couple things to zero*/
   if (st->submodes[st->submodeID] == NULL)
   {
      spx_coef_t *lpc;
      lpc = PUSH(stack,11, spx_coef_t);
      bw_lpc(.93, st->interp_qlpc, lpc, 10);
      /*for (i=0;i<st->frameSize;i++)
        st->exc[i]=0;*/
      {
         float innov_gain=0;
         float pgain=st->last_pitch_gain;
         if (pgain>.6)
            pgain=.6;
         for (i=0;i<st->frameSize;i++)
            innov_gain += st->innov[i]*st->innov[i];
         innov_gain=sqrt(innov_gain/st->frameSize);
         for (i=0;i<st->frameSize;i++)
            st->exc[i]=0;
         speex_rand_vec(innov_gain, st->exc, st->frameSize);
      }


      st->first=1;

      /* Final signal synthesis from excitation */
      iir_mem2(st->exc, lpc, st->frame, st->frameSize, st->lpcSize, st->mem_sp);

      for (i=0;i<st->frameSize;i++)
      {
         spx_word32_t sig = PSHR(st->frame[i],SIG_SHIFT);
         if (sig>32767)
            sig = 32767;
         if (sig<-32767)
            sig = -32767;
         out[i]=sig;
      }

      st->count_lost=0;
      return 0;
   }

   /* Unquantize LSPs */
   SUBMODE(lsp_unquant)(st->qlsp, st->lpcSize, bits);

   /*Damp memory if a frame was lost and the LSP changed too much*/
   if (st->count_lost)
   {
      float lsp_dist=0, fact;
      for (i=0;i<st->lpcSize;i++)
         lsp_dist += fabs(st->old_qlsp[i] - st->qlsp[i]);
      lsp_dist /= LSP_SCALING*LSP_SCALING;
      fact = .6*exp(-.2*lsp_dist);
      for (i=0;i<2*st->lpcSize;i++)
         st->mem_sp[i] = (spx_mem_t)(st->mem_sp[i]*fact);
   }


   /* Handle first frame and lost-packet case */
   if (st->first || st->count_lost)
   {
      for (i=0;i<st->lpcSize;i++)
         st->old_qlsp[i] = st->qlsp[i];
   }

#ifdef EPIC_48K
   if (st->lbr_48k) {
      pitch_half[0] = st->min_pitch+speex_bits_unpack_unsigned(bits, 7);
      pitch_half[1] = pitch_half[0]+speex_bits_unpack_unsigned(bits, 2)-1;

      ol_pitch_id = speex_bits_unpack_unsigned(bits, 3);
      ol_pitch_coef=0.13514*ol_pitch_id;

      {
         int qe;
         qe = speex_bits_unpack_unsigned(bits, 4);
         ol_gain = SIG_SCALING*exp((qe+2)/2.1),SIG_SHIFT;
      }

   } else {
#endif

   /* Get open-loop pitch estimation for low bit-rate pitch coding */
   if (SUBMODE(lbr_pitch)!=-1)
   {
      ol_pitch = st->min_pitch+speex_bits_unpack_unsigned(bits, 7);
   } 
   
   if (SUBMODE(forced_pitch_gain))
   {
      int quant;
      quant = speex_bits_unpack_unsigned(bits, 4);
      ol_pitch_coef=0.066667*quant;
   }
   
   /* Get global excitation gain */
   {
      int qe;
      qe = speex_bits_unpack_unsigned(bits, 5);
#ifdef FIXED_POINT
      ol_gain = MULT16_32_Q15(28406,ol_gain_table[qe]);
#else
      ol_gain = SIG_SCALING*exp(qe/3.5);
#endif
   }
#ifdef EPIC_48K
   }
#endif

   awk1=PUSH(stack, st->lpcSize+1, spx_coef_t);
   awk2=PUSH(stack, st->lpcSize+1, spx_coef_t);
   awk3=PUSH(stack, st->lpcSize+1, spx_coef_t);

   if (st->submodeID==1)
   {
      int extra;
      extra = speex_bits_unpack_unsigned(bits, 4);

      if (extra==15)
         st->dtx_enabled=1;
      else
         st->dtx_enabled=0;
   }
   if (st->submodeID>1)
      st->dtx_enabled=0;

   /*Loop on subframes */
   for (sub=0;sub<st->nbSubframes;sub++)
   {
      int offset;
      spx_sig_t *sp, *exc;
      float tmp;

#ifdef EPIC_48K
      if (st->lbr_48k)
      {
         if (sub*2 < st->nbSubframes)
            ol_pitch = pitch_half[0];
         else
            ol_pitch = pitch_half[1];
      }
#endif

      /* Offset relative to start of frame */
      offset = st->subframeSize*sub;
      /* Original signal */
      sp=st->frame+offset;
      /* Excitation */
      exc=st->exc+offset;
      /* Excitation after post-filter*/

      /* LSP interpolation (quantized and unquantized) */
      lsp_interpolate(st->old_qlsp, st->qlsp, st->interp_qlsp, st->lpcSize, sub, st->nbSubframes);

      /* Make sure the LSP's are stable */
      lsp_enforce_margin(st->interp_qlsp, st->lpcSize, LSP_MARGIN);


      /* Compute interpolated LPCs (unquantized) */
      lsp_to_lpc(st->interp_qlsp, st->interp_qlpc, st->lpcSize, stack);

      /* Compute enhanced synthesis filter */
      if (st->lpc_enh_enabled)
      {
         float r=.9;
         
         float k1,k2,k3;
         k1=SUBMODE(lpc_enh_k1);
         k2=SUBMODE(lpc_enh_k2);
         k3=(1-(1-r*k1)/(1-r*k2))/r;
         if (!st->lpc_enh_enabled)
         {
            k1=k2;
            k3=0;
         }
         bw_lpc(k1, st->interp_qlpc, awk1, st->lpcSize);
         bw_lpc(k2, st->interp_qlpc, awk2, st->lpcSize);
         bw_lpc(k3, st->interp_qlpc, awk3, st->lpcSize);
         
      }

      /* Compute analysis filter at w=pi */
      {
         spx_word32_t pi_g=st->interp_qlpc[0];
         for (i=1;i<=st->lpcSize;i+=2)
         {
            pi_g += -st->interp_qlpc[i] +  st->interp_qlpc[i+1];
         }
         st->pi_gain[sub] = pi_g;
      }

      /* Reset excitation */
      for (i=0;i<st->subframeSize;i++)
         exc[i]=0;

      /*Adaptive codebook contribution*/
      if (SUBMODE(ltp_unquant))
      {
         int pit_min, pit_max;
         /* Handle pitch constraints if any */
         if (SUBMODE(lbr_pitch) != -1)
         {
            int margin;
            margin = SUBMODE(lbr_pitch);
            if (margin)
            {
/* GT - need optimization?
               if (ol_pitch < st->min_pitch+margin-1)
                  ol_pitch=st->min_pitch+margin-1;
               if (ol_pitch > st->max_pitch-margin)
                  ol_pitch=st->max_pitch-margin;
               pit_min = ol_pitch-margin+1;
               pit_max = ol_pitch+margin;
*/
               pit_min = ol_pitch-margin+1;
               if (pit_min < st->min_pitch)
		  pit_min = st->min_pitch;
               pit_max = ol_pitch+margin;
               if (pit_max > st->max_pitch)
		  pit_max = st->max_pitch;
            } else {
               pit_min = pit_max = ol_pitch;
            }
         } else {
            pit_min = st->min_pitch;
            pit_max = st->max_pitch;
         }


#ifdef EPIC_48K
         if (st->lbr_48k)
         {
             SUBMODE(ltp_unquant)(exc, pit_min, pit_max, ol_pitch_coef, SUBMODE(ltp_params), 
                                  st->subframeSize, &pitch, &pitch_gain[0], bits, stack, 
                                  st->count_lost, offset, st->last_pitch_gain, ol_pitch_id);
         } else {
#endif

             SUBMODE(ltp_unquant)(exc, pit_min, pit_max, ol_pitch_coef, SUBMODE(ltp_params), 
                                  st->subframeSize, &pitch, &pitch_gain[0], bits, stack, 
                                  st->count_lost, offset, st->last_pitch_gain, 0);

#ifdef EPIC_48K
         }
#endif

         
         /* If we had lost frames, check energy of last received frame */
         if (st->count_lost && ol_gain < st->last_ol_gain)
         {
            float fact = (float)ol_gain/(st->last_ol_gain+1);
            for (i=0;i<st->subframeSize;i++)
               exc[i]*=fact;
         }

         tmp = fabs(pitch_gain[0]+pitch_gain[1]+pitch_gain[2]);
         tmp = fabs(pitch_gain[1]);
         if (pitch_gain[0]>0)
            tmp += pitch_gain[0];
         else
            tmp -= .5*pitch_gain[0];
         if (pitch_gain[2]>0)
            tmp += pitch_gain[2];
         else
            tmp -= .5*pitch_gain[0];


         pitch_average += tmp;
         if (tmp>best_pitch_gain)
         {
            best_pitch = pitch;
	    best_pitch_gain = tmp;
	    /*            best_pitch_gain = tmp*.9;
	                if (best_pitch_gain>.85)
                        best_pitch_gain=.85;*/
         }
      } else {
         speex_error("No pitch prediction, what's wrong");
      }
      
      /* Unquantize the innovation */
      {
         int q_energy;
         spx_word32_t ener;
         spx_sig_t *innov;
         
         innov = st->innov+sub*st->subframeSize;
         for (i=0;i<st->subframeSize;i++)
            innov[i]=0;

         /* Decode sub-frame gain correction */
         if (SUBMODE(have_subframe_gain)==3)
         {
            q_energy = speex_bits_unpack_unsigned(bits, 3);
            ener = MULT16_32_Q14(exc_gain_quant_scal3[q_energy],ol_gain);
         } else if (SUBMODE(have_subframe_gain)==1)
         {
            q_energy = speex_bits_unpack_unsigned(bits, 1);
            ener = MULT16_32_Q14(exc_gain_quant_scal1[q_energy],ol_gain);
         } else {
            ener = ol_gain;
         }
                  
         if (SUBMODE(innovation_unquant))
         {
            /*Fixed codebook contribution*/
            SUBMODE(innovation_unquant)(innov, SUBMODE(innovation_params), st->subframeSize, bits, stack);
         } else {
            speex_error("No fixed codebook");
         }

         /* De-normalize innovation and update excitation */
#ifdef FIXED_POINT
         signal_mul(innov, innov, ener, st->subframeSize);
#else
         signal_mul(innov, innov, ener, st->subframeSize);
#endif
         /*Vocoder mode*/
         if (st->submodeID==1) 
         {
            float g=ol_pitch_coef;

            
            for (i=0;i<st->subframeSize;i++)
               exc[i]=0;
            while (st->voc_offset<st->subframeSize)
            {
               if (st->voc_offset>=0)
                  exc[st->voc_offset]=SIG_SCALING*sqrt(1.0*ol_pitch);
               st->voc_offset+=ol_pitch;
            }
            st->voc_offset -= st->subframeSize;

            g=.5+2*(g-.6);
            if (g<0)
               g=0;
            if (g>1)
               g=1;
            for (i=0;i<st->subframeSize;i++)
            {
               float exci=exc[i];
               exc[i]=.8*g*exc[i]*ol_gain/SIG_SCALING + .6*g*st->voc_m1*ol_gain/SIG_SCALING + .5*g*innov[i] - .5*g*st->voc_m2 + (1-g)*innov[i];
               st->voc_m1 = exci;
               st->voc_m2=innov[i];
               st->voc_mean = .95*st->voc_mean + .05*exc[i];
               exc[i]-=st->voc_mean;
            }
         } else {
            for (i=0;i<st->subframeSize;i++)
               exc[i]+=innov[i];
            /*print_vec(exc, 40, "innov");*/
         }
         /* Decode second codebook (only for some modes) */
         if (SUBMODE(double_codebook))
         {
            char *tmp_stack=stack;
            spx_sig_t *innov2 = PUSH(tmp_stack, st->subframeSize, spx_sig_t);
            for (i=0;i<st->subframeSize;i++)
               innov2[i]=0;
            SUBMODE(innovation_unquant)(innov2, SUBMODE(innovation_params), st->subframeSize, bits, tmp_stack);
            signal_mul(innov2, innov2, ener*(1/2.2), st->subframeSize);
            for (i=0;i<st->subframeSize;i++)
               exc[i] += innov2[i];
         }

      }

      for (i=0;i<st->subframeSize;i++)
         sp[i]=exc[i];

      /* Signal synthesis */
      if (st->lpc_enh_enabled && SUBMODE(comb_gain)>0)
         comb_filter(exc, sp, st->interp_qlpc, st->lpcSize, st->subframeSize,
                              pitch, pitch_gain, SUBMODE(comb_gain), st->comb_mem);

      if (st->lpc_enh_enabled)
      {
         /* Use enhanced LPC filter */
         filter_mem2(sp, awk2, awk1, sp, st->subframeSize, st->lpcSize, 
                     st->mem_sp+st->lpcSize);
         filter_mem2(sp, awk3, st->interp_qlpc, sp, st->subframeSize, st->lpcSize, 
                     st->mem_sp);
      } else {
         /* Use regular filter */
         for (i=0;i<st->lpcSize;i++)
            st->mem_sp[st->lpcSize+i] = 0;
         iir_mem2(sp, st->interp_qlpc, sp, st->subframeSize, st->lpcSize, 
                     st->mem_sp);
      }
   }
   
   /*Copy output signal*/   
   for (i=0;i<st->frameSize;i++)
   {
      spx_word32_t sig = PSHR(st->frame[i],SIG_SHIFT);
      if (sig>32767)
         sig = 32767;
      if (sig<-32767)
         sig = -32767;
     out[i]=sig;
   }

   /*for (i=0;i<st->frameSize;i++)
     printf ("%d\n", (int)st->frame[i]);*/

   /* Store the LSPs for interpolation in the next frame */
   for (i=0;i<st->lpcSize;i++)
      st->old_qlsp[i] = st->qlsp[i];

   /* The next frame will not be the first (Duh!) */
   st->first = 0;
   st->count_lost=0;
   st->last_pitch = best_pitch;
   st->last_pitch_gain = .25*pitch_average;
   st->pitch_gain_buf[st->pitch_gain_buf_idx++] = st->last_pitch_gain;
   if (st->pitch_gain_buf_idx > 2) /* rollover */
      st->pitch_gain_buf_idx = 0;

   st->last_ol_gain = ol_gain;

   return 0;
}

int nb_encoder_ctl(void *state, int request, void *ptr)
{
   EncState *st;
   st=(EncState*)state;     
   switch(request)
   {
   case SPEEX_GET_FRAME_SIZE:
      (*(int*)ptr) = st->frameSize;
      break;
   case SPEEX_SET_LOW_MODE:
   case SPEEX_SET_MODE:
      st->submodeSelect = st->submodeID = (*(int*)ptr);
      break;
   case SPEEX_GET_LOW_MODE:
   case SPEEX_GET_MODE:
      (*(int*)ptr) = st->submodeID;
      break;
   case SPEEX_SET_VBR:
      st->vbr_enabled = (*(int*)ptr);
      break;
   case SPEEX_GET_VBR:
      (*(int*)ptr) = st->vbr_enabled;
      break;
   case SPEEX_SET_VAD:
      st->vad_enabled = (*(int*)ptr);
      break;
   case SPEEX_GET_VAD:
      (*(int*)ptr) = st->vad_enabled;
      break;
   case SPEEX_SET_DTX:
      st->dtx_enabled = (*(int*)ptr);
      break;
   case SPEEX_GET_DTX:
      (*(int*)ptr) = st->dtx_enabled;
      break;
   case SPEEX_SET_ABR:
      st->abr_enabled = (*(int*)ptr);
      st->vbr_enabled = 1;
      {
         int i=10, rate, target;
         float vbr_qual;
         target = (*(int*)ptr);
         while (i>=0)
         {
            speex_encoder_ctl(st, SPEEX_SET_QUALITY, &i);
            speex_encoder_ctl(st, SPEEX_GET_BITRATE, &rate);
            if (rate <= target)
               break;
            i--;
         }
         vbr_qual=i;
         if (vbr_qual<0)
            vbr_qual=0;
         speex_encoder_ctl(st, SPEEX_SET_VBR_QUALITY, &vbr_qual);
         st->abr_count=0;
         st->abr_drift=0;
         st->abr_drift2=0;
      }
      
      break;
   case SPEEX_GET_ABR:
      (*(int*)ptr) = st->abr_enabled;
      break;
   case SPEEX_SET_VBR_QUALITY:
      st->vbr_quality = (*(float*)ptr);
      break;
   case SPEEX_GET_VBR_QUALITY:
      (*(float*)ptr) = st->vbr_quality;
      break;
   case SPEEX_SET_QUALITY:
      {
         int quality = (*(int*)ptr);
         if (quality < 0)
            quality = 0;
         if (quality > 10)
            quality = 10;
         st->submodeSelect = st->submodeID = ((SpeexNBMode*)(st->mode->mode))->quality_map[quality];
      }
      break;
   case SPEEX_SET_COMPLEXITY:
      st->complexity = (*(int*)ptr);
      if (st->complexity<1)
         st->complexity=1;
      break;
   case SPEEX_GET_COMPLEXITY:
      (*(int*)ptr) = st->complexity;
      break;
   case SPEEX_SET_BITRATE:
      {
         int i=10, rate, target;
         target = (*(int*)ptr);
         while (i>=0)
         {
            speex_encoder_ctl(st, SPEEX_SET_QUALITY, &i);
            speex_encoder_ctl(st, SPEEX_GET_BITRATE, &rate);
            if (rate <= target)
               break;
            i--;
         }
      }
      break;
   case SPEEX_GET_BITRATE:
      if (st->submodes[st->submodeID])
         (*(int*)ptr) = st->sampling_rate*SUBMODE(bits_per_frame)/st->frameSize;
      else
         (*(int*)ptr) = st->sampling_rate*(NB_SUBMODE_BITS+1)/st->frameSize;
      break;
   case SPEEX_SET_SAMPLING_RATE:
      st->sampling_rate = (*(int*)ptr);
      break;
   case SPEEX_GET_SAMPLING_RATE:
      (*(int*)ptr)=st->sampling_rate;
      break;
   case SPEEX_RESET_STATE:
      {
         int i;
         st->bounded_pitch = 1;
         st->first = 1;
         for (i=0;i<st->lpcSize;i++)
            st->lsp[i]=(M_PI*((float)(i+1)))/(st->lpcSize+1);
         for (i=0;i<st->lpcSize;i++)
            st->mem_sw[i]=st->mem_sw_whole[i]=st->mem_sp[i]=st->mem_exc[i]=0;
         for (i=0;i<st->bufSize;i++)
            st->excBuf[i]=st->swBuf[i]=st->inBuf[i]=st->exc2Buf[i]=0;
      }
      break;
   case SPEEX_SET_SUBMODE_ENCODING:
      st->encode_submode = (*(int*)ptr);
      break;
   case SPEEX_GET_SUBMODE_ENCODING:
      (*(int*)ptr) = st->encode_submode;
      break;
   case SPEEX_GET_PI_GAIN:
      {
         int i;
         spx_word32_t *g = (spx_word32_t*)ptr;
         for (i=0;i<st->nbSubframes;i++)
            g[i]=st->pi_gain[i];
      }
      break;
   case SPEEX_GET_EXC:
      {
         int i;
         spx_sig_t *e = (spx_sig_t*)ptr;
         for (i=0;i<st->frameSize;i++)
            e[i]=st->exc[i];
      }
      break;
   case SPEEX_GET_INNOV:
      {
         int i;
         spx_sig_t *e = (spx_sig_t*)ptr;
         for (i=0;i<st->frameSize;i++)
            e[i]=st->innov[i];
      }
      break;
   case SPEEX_GET_RELATIVE_QUALITY:
      (*(float*)ptr)=st->relative_quality;
      break;
   default:
      speex_warning_int("Unknown nb_ctl request: ", request);
      return -1;
   }
   return 0;
}

int nb_decoder_ctl(void *state, int request, void *ptr)
{
   DecState *st;
   st=(DecState*)state;
   switch(request)
   {
   case SPEEX_SET_LOW_MODE:
   case SPEEX_SET_MODE:
      st->submodeID = (*(int*)ptr);
      break;
   case SPEEX_GET_LOW_MODE:
   case SPEEX_GET_MODE:
      (*(int*)ptr) = st->submodeID;
      break;
   case SPEEX_SET_ENH:
      st->lpc_enh_enabled = *((int*)ptr);
      break;
   case SPEEX_GET_ENH:
      *((int*)ptr) = st->lpc_enh_enabled;
      break;
   case SPEEX_GET_FRAME_SIZE:
      (*(int*)ptr) = st->frameSize;
      break;
   case SPEEX_GET_BITRATE:
      if (st->submodes[st->submodeID])
         (*(int*)ptr) = st->sampling_rate*SUBMODE(bits_per_frame)/st->frameSize;
      else
         (*(int*)ptr) = st->sampling_rate*(NB_SUBMODE_BITS+1)/st->frameSize;
      break;
   case SPEEX_SET_SAMPLING_RATE:
      st->sampling_rate = (*(int*)ptr);
      break;
   case SPEEX_GET_SAMPLING_RATE:
      (*(int*)ptr)=st->sampling_rate;
      break;
   case SPEEX_SET_HANDLER:
      {
         SpeexCallback *c = (SpeexCallback*)ptr;
         st->speex_callbacks[c->callback_id].func=c->func;
         st->speex_callbacks[c->callback_id].data=c->data;
         st->speex_callbacks[c->callback_id].callback_id=c->callback_id;
      }
      break;
   case SPEEX_SET_USER_HANDLER:
      {
         SpeexCallback *c = (SpeexCallback*)ptr;
         st->user_callback.func=c->func;
         st->user_callback.data=c->data;
         st->user_callback.callback_id=c->callback_id;
      }
      break;
   case SPEEX_RESET_STATE:
      {
         int i;
         for (i=0;i<2*st->lpcSize;i++)
            st->mem_sp[i]=0;
         for (i=0;i<st->bufSize;i++)
            st->excBuf[i]=st->inBuf[i]=0;
      }
      break;
   case SPEEX_SET_SUBMODE_ENCODING:
      st->encode_submode = (*(int*)ptr);
      break;
   case SPEEX_GET_SUBMODE_ENCODING:
      (*(int*)ptr) = st->encode_submode;
      break;
   case SPEEX_GET_PI_GAIN:
      {
         int i;
         spx_word32_t *g = (spx_word32_t*)ptr;
         for (i=0;i<st->nbSubframes;i++)
            g[i]=st->pi_gain[i];
      }
      break;
   case SPEEX_GET_EXC:
      {
         int i;
         spx_sig_t *e = (spx_sig_t*)ptr;
         for (i=0;i<st->frameSize;i++)
            e[i]=st->exc[i];
      }
      break;
   case SPEEX_GET_INNOV:
      {
         int i;
         spx_sig_t *e = (spx_sig_t*)ptr;
         for (i=0;i<st->frameSize;i++)
            e[i]=st->innov[i];
      }
      break;
   case SPEEX_GET_DTX_STATUS:
      *((int*)ptr) = st->dtx_enabled;
      break;
   default:
      speex_warning_int("Unknown nb_ctl request: ", request);
      return -1;
   }
   return 0;
}
