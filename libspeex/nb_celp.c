/* Copyright (C) 2002 Jean-Marc Valin 
   File: nb_celp.c

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

#include <stdlib.h>
#include <stdio.h>
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

extern int training_weight;
#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif

#define SUBMODE(x) st->submodes[st->submodeID]->x

float exc_gain_quant_scal3[8]={-2.794750, -1.810660, -1.169850, -0.848119, -0.587190, -0.329818, -0.063266, 0.282826};

float exc_gain_quant_scal1[2]={-0.35, 0.05};
/*float exc_gain_quant_scal1[2]={-0.35, 0.05};*/

#define sqr(x) ((x)*(x))
#define min(a,b) ((a) < (b) ? (a) : (b))

void *nb_encoder_init(SpeexMode *m)
{
   EncState *st;
   SpeexNBMode *mode;
   int i;

   mode=m->mode;
   st = speex_alloc(sizeof(EncState));
   st->mode=m;
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
   st->lag_factor=mode->lag_factor;
   st->lpc_floor = mode->lpc_floor;
   st->preemph = mode->preemph;
  
   st->submodes=mode->submodes;
   st->submodeID=mode->defaultSubmode;
   st->pre_mem=0;
   st->pre_mem2=0;

   /* Allocating input buffer */
   st->inBuf = speex_alloc(st->bufSize*sizeof(float));
   st->frame = st->inBuf + st->bufSize - st->windowSize;
   /* Allocating excitation buffer */
   st->excBuf = speex_alloc(st->bufSize*sizeof(float));
   st->exc = st->excBuf + st->bufSize - st->windowSize;
   st->swBuf = speex_alloc(st->bufSize*sizeof(float));
   st->sw = st->swBuf + st->bufSize - st->windowSize;

   st->exc2Buf = speex_alloc(st->bufSize*sizeof(float));
   st->exc2 = st->exc2Buf + st->bufSize - st->windowSize;

   st->innov = speex_alloc(st->frameSize*sizeof(float));

   /* Asymetric "pseudo-Hamming" window */
   {
      int part1, part2;
      part1 = st->subframeSize*7/2;
      part2 = st->subframeSize*5/2;
      st->window = speex_alloc(st->windowSize*sizeof(float));
      for (i=0;i<part1;i++)
         st->window[i]=.54-.46*cos(M_PI*i/part1);
      for (i=0;i<part2;i++)
         st->window[part1+i]=.54+.46*cos(M_PI*i/part2);
   }
   /* Create the window for autocorrelation (lag-windowing) */
   st->lagWindow = speex_alloc((st->lpcSize+1)*sizeof(float));
   for (i=0;i<st->lpcSize+1;i++)
      st->lagWindow[i]=exp(-.5*sqr(2*M_PI*st->lag_factor*i));

   st->autocorr = speex_alloc((st->lpcSize+1)*sizeof(float));

   st->stack = speex_alloc(20000*sizeof(float));

   st->buf2 = speex_alloc(st->windowSize*sizeof(float));

   st->lpc = speex_alloc((st->lpcSize+1)*sizeof(float));
   st->interp_lpc = speex_alloc((st->lpcSize+1)*sizeof(float));
   st->interp_qlpc = speex_alloc((st->lpcSize+1)*sizeof(float));
   st->bw_lpc1 = speex_alloc((st->lpcSize+1)*sizeof(float));
   st->bw_lpc2 = speex_alloc((st->lpcSize+1)*sizeof(float));

   st->lsp = speex_alloc(st->lpcSize*sizeof(float));
   st->qlsp = speex_alloc(st->lpcSize*sizeof(float));
   st->old_lsp = speex_alloc(st->lpcSize*sizeof(float));
   st->old_qlsp = speex_alloc(st->lpcSize*sizeof(float));
   st->interp_lsp = speex_alloc(st->lpcSize*sizeof(float));
   st->interp_qlsp = speex_alloc(st->lpcSize*sizeof(float));
   st->rc = speex_alloc(st->lpcSize*sizeof(float));
   st->first = 1;

   st->mem_sp = speex_alloc(st->lpcSize*sizeof(float));
   st->mem_sw = speex_alloc(st->lpcSize*sizeof(float));
   st->mem_sw_whole = speex_alloc(st->lpcSize*sizeof(float));
   st->mem_exc = speex_alloc(st->lpcSize*sizeof(float));

   st->pi_gain = speex_alloc(st->nbSubframes*sizeof(float));

   st->pitch = speex_alloc(st->nbSubframes*sizeof(int));

   if (1) {
      st->vbr = speex_alloc(sizeof(VBRState));
      vbr_init(st->vbr);
      st->vbr_quality = 8;
      st->vbr_enabled = 0;
   } else {
      st->vbr = 0;
   }
   st->complexity=2;

   return st;
}

void nb_encoder_destroy(void *state)
{
   EncState *st=state;
   /* Free all allocated memory */
   speex_free(st->inBuf);
   speex_free(st->excBuf);
   speex_free(st->swBuf);
   speex_free(st->exc2Buf);
   speex_free(st->innov);
   speex_free(st->stack);

   speex_free(st->window);
   speex_free(st->buf2);
   speex_free(st->lpc);
   speex_free(st->interp_lpc);
   speex_free(st->interp_qlpc);
   
   speex_free(st->bw_lpc1);
   speex_free(st->bw_lpc2);
   speex_free(st->autocorr);
   speex_free(st->lagWindow);
   speex_free(st->lsp);
   speex_free(st->qlsp);
   speex_free(st->old_lsp);
   speex_free(st->interp_lsp);
   speex_free(st->old_qlsp);
   speex_free(st->interp_qlsp);
   speex_free(st->rc);

   speex_free(st->mem_sp);
   speex_free(st->mem_sw);
   speex_free(st->mem_sw_whole);
   speex_free(st->mem_exc);
   speex_free(st->pi_gain);
   speex_free(st->pitch);

   vbr_destroy(st->vbr);
   speex_free(st->vbr);

   /*Free state memory... should be last*/
   speex_free(st);
}

void nb_encode(void *state, float *in, SpeexBits *bits)
{
   EncState *st;
   int i, sub, roots;
   float error;
   int ol_pitch;
   float ol_pitch_coef;
   float ol_gain;
   float delta_qual=0;

   st=state;
   
   /* Copy new data in input buffer */
   speex_move(st->inBuf, st->inBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(float));
   st->inBuf[st->bufSize-st->frameSize] = in[0] - st->preemph*st->pre_mem;
   for (i=1;i<st->frameSize;i++)
      st->inBuf[st->bufSize-st->frameSize+i] = in[i] - st->preemph*in[i-1];
   st->pre_mem = in[st->frameSize-1];

   /* Move signals 1 frame towards the past */
   speex_move(st->exc2Buf, st->exc2Buf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(float));
   speex_move(st->excBuf, st->excBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(float));
   speex_move(st->swBuf, st->swBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(float));


   /* Window for analysis */
   for (i=0;i<st->windowSize;i++)
      st->buf2[i] = st->frame[i] * st->window[i];

   /* Compute auto-correlation */
   autocorr(st->buf2, st->autocorr, st->lpcSize+1, st->windowSize);

   st->autocorr[0] += 10;        /* prevents NANs */
   st->autocorr[0] *= st->lpc_floor; /* Noise floor in auto-correlation domain */

   /* Lag windowing: equivalent to filtering in the power-spectrum domain */
   for (i=0;i<st->lpcSize+1;i++)
      st->autocorr[i] *= st->lagWindow[i];

   /* Levinson-Durbin */
   error = wld(st->lpc+1, st->autocorr, st->rc, st->lpcSize);
   st->lpc[0]=1;

   /* LPC to LSPs (x-domain) transform */
   roots=lpc_to_lsp (st->lpc, st->lpcSize, st->lsp, 15, 0.2, st->stack);
   if (roots!=st->lpcSize)
   {
      roots = lpc_to_lsp (st->lpc, st->lpcSize, st->lsp, 11, 0.02, st->stack);
      if (roots!=st->lpcSize) {
         /*fprintf (stderr, "roots!=st->lpcSize (found only %d roots)\n", roots);*/

         /*If we can't find all LSP's, do some damage control and use a flat filter*/
         for (i=0;i<st->lpcSize;i++)
         {
            st->lsp[i]=cos(M_PI*((float)(i+1))/(st->lpcSize+1));
         }
      }
   }


   /* LSP x-domain to angle domain*/
   for (i=0;i<st->lpcSize;i++)
      st->lsp[i] = acos(st->lsp[i]);
   /*print_vec(st->lsp, 10, "LSP:");*/
   /* LSP Quantization */
   if (st->first)
   {
      for (i=0;i<st->lpcSize;i++)
         st->old_lsp[i] = st->lsp[i];
   }

   if (0) {
      float dd=0;
      for (i=0;i<st->lpcSize;i++)
         dd += fabs(st->old_lsp[i] - st->lsp[i]);
      printf ("lspdist = %f\n", dd);
   }

   /* Whole frame analysis (open-loop estimation of pitch and excitation gain) */
   {
      for (i=0;i<st->lpcSize;i++)
         st->interp_lsp[i] = .5*st->old_lsp[i] + .5*st->lsp[i];

      lsp_enforce_margin(st->interp_lsp, st->lpcSize, .002);

      /* Compute interpolated LPCs (unquantized) for whole frame*/
      for (i=0;i<st->lpcSize;i++)
         st->interp_lsp[i] = cos(st->interp_lsp[i]);
      lsp_to_lpc(st->interp_lsp, st->interp_lpc, st->lpcSize,st->stack);

      bw_lpc(st->gamma1, st->interp_lpc, st->bw_lpc1, st->lpcSize);
      bw_lpc(st->gamma2, st->interp_lpc, st->bw_lpc2, st->lpcSize);

      filter_mem2(st->frame, st->bw_lpc1, st->bw_lpc2, st->sw, st->frameSize, st->lpcSize, st->mem_sw_whole);

      /*Open-loop pitch*/
      {
         int nol_pitch[4];
         float nol_pitch_coef[4];
         open_loop_nbest_pitch(st->sw, st->min_pitch, st->max_pitch, st->frameSize, 
                               nol_pitch, nol_pitch_coef, 4, st->stack);
         ol_pitch=nol_pitch[0];
         ol_pitch_coef = nol_pitch_coef[0];
         /*Try to remove pitch multiples*/
         for (i=1;i<4;i++)
         {
            if ((nol_pitch_coef[i] > .85*ol_pitch_coef) && 
                (fabs(2*nol_pitch[i]-ol_pitch)<=2 || fabs(3*nol_pitch[i]-ol_pitch)<=4 || 
                 fabs(4*nol_pitch[i]-ol_pitch)<=6 || fabs(5*nol_pitch[i]-ol_pitch)<=8))
            {
               /*ol_pitch_coef=nol_pitch_coef[i];*/
               ol_pitch = nol_pitch[i];
            }
         }
         /*ol_pitch_coef = sqrt(ol_pitch_coef);*/
         /*printf ("ol_pitch: %d %f\n", ol_pitch, ol_pitch_coef);*/
      }
      /*Compute "real" excitation*/
      fir_mem2(st->frame, st->interp_lpc, st->exc, st->frameSize, st->lpcSize, st->mem_exc);

      /* Compute open-loop excitation gain */
      ol_gain=0;
      for (i=0;i<st->frameSize;i++)
         ol_gain += st->exc[i]*st->exc[i];
      
      ol_gain=sqrt(1+ol_gain/st->frameSize);
   }

   /*Experimental VBR stuff*/
   if (st->vbr)
   {
      delta_qual = vbr_analysis(st->vbr, in, st->frameSize, ol_pitch, ol_pitch_coef);
      /*if (delta_qual<0)*/
         delta_qual*=.1*(3+st->vbr_quality);
      if (st->vbr_enabled) 
      {
         int qual = (int)floor(st->vbr_quality+delta_qual+.5);
         if (qual<1 && delta_qual>-3.5)
            qual=1;
         if (qual<0)
            qual=0;
         if (qual>10)
            qual=10;
         if (qual==10 && st->vbr_quality<10)
            qual=9;
         speex_encoder_ctl(state, SPEEX_SET_QUALITY, &qual);
      }
   }
   /*printf ("VBR quality = %f\n", vbr_qual);*/

   /* First, transmit a zero for narrowband */
   speex_bits_pack(bits, 0, 1);

   /* Transmit the sub-mode we use for this frame */
   speex_bits_pack(bits, st->submodeID, NB_SUBMODE_BITS);


   /* If null mode (no transmission), just set a couple things to zero*/
   if (st->submodes[st->submodeID] == NULL)
   {
      for (i=0;i<st->frameSize;i++)
         st->exc[i]=st->exc2[i]=st->sw[i]=0;

      for (i=0;i<st->lpcSize;i++)
         st->mem_sw[i]=0;
      st->first=1;

      /* Final signal synthesis from excitation */
      iir_mem2(st->exc, st->interp_qlpc, st->frame, st->subframeSize, st->lpcSize, st->mem_sp);

      in[0] = st->frame[0] + st->preemph*st->pre_mem2;
      for (i=1;i<st->frameSize;i++)
         in[i]=st->frame[i] + st->preemph*in[i-1];
      st->pre_mem2=in[st->frameSize-1];

      return;

   }

   /*Quantize LSPs*/
#if 1 /*0 for unquantized*/
   SUBMODE(lsp_quant)(st->lsp, st->qlsp, st->lpcSize, bits);
#else
   for (i=0;i<st->lpcSize;i++)
     st->qlsp[i]=st->lsp[i];
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
         quant=0;
      if (quant<0)
         quant=0;
      speex_bits_pack(bits, quant, 4);
      ol_pitch_coef=0.066667*quant;
   }
   
   
   /*Quantize and transmit open-loop excitation gain*/
   {
      int qe = (int)(floor(3.5*log(ol_gain)));
      if (qe<0)
         qe=0;
      if (qe>31)
         qe=31;
      ol_gain = exp(qe/3.5);
      speex_bits_pack(bits, qe, 5);
   }

   /* Special case for first frame */
   if (st->first)
   {
      for (i=0;i<st->lpcSize;i++)
         st->old_qlsp[i] = st->qlsp[i];
   }

   /* Loop on sub-frames */
   for (sub=0;sub<st->nbSubframes;sub++)
   {
      float esig, enoise, snr, tmp;
      int   offset;
      float *sp, *sw, *res, *exc, *target, *mem, *exc2;
      int pitch;

      /* Offset relative to start of frame */
      offset = st->subframeSize*sub;
      /* Original signal */
      sp=st->frame+offset;
      /* Excitation */
      exc=st->exc+offset;
      /* Weighted signal */
      sw=st->sw+offset;

      exc2=st->exc2+offset;

      /* Filter response */
      res = PUSH(st->stack, st->subframeSize);
      /* Target signal */
      target = PUSH(st->stack, st->subframeSize);
      mem = PUSH(st->stack, st->lpcSize);

      /* LSP interpolation (quantized and unquantized) */
      tmp = (1.0 + sub)/st->nbSubframes;
      for (i=0;i<st->lpcSize;i++)
         st->interp_lsp[i] = (1-tmp)*st->old_lsp[i] + tmp*st->lsp[i];
      for (i=0;i<st->lpcSize;i++)
         st->interp_qlsp[i] = (1-tmp)*st->old_qlsp[i] + tmp*st->qlsp[i];

      /* Make sure the filters are stable */
      lsp_enforce_margin(st->interp_lsp, st->lpcSize, .002);
      lsp_enforce_margin(st->interp_qlsp, st->lpcSize, .002);

      /* Compute interpolated LPCs (quantized and unquantized) */
      for (i=0;i<st->lpcSize;i++)
         st->interp_lsp[i] = cos(st->interp_lsp[i]);
      lsp_to_lpc(st->interp_lsp, st->interp_lpc, st->lpcSize,st->stack);

      for (i=0;i<st->lpcSize;i++)
         st->interp_qlsp[i] = cos(st->interp_qlsp[i]);
      lsp_to_lpc(st->interp_qlsp, st->interp_qlpc, st->lpcSize, st->stack);

      /* Compute analysis filter gain at w=pi (for use in SB-CELP) */
      tmp=1;
      st->pi_gain[sub]=0;
      for (i=0;i<=st->lpcSize;i++)
      {
         st->pi_gain[sub] += tmp*st->interp_qlpc[i];
         tmp = -tmp;
      }
     

      /* Compute bandwidth-expanded (unquantized) LPCs for perceptual weighting */
      bw_lpc(st->gamma1, st->interp_lpc, st->bw_lpc1, st->lpcSize);
      if (st->gamma2>=0)
         bw_lpc(st->gamma2, st->interp_lpc, st->bw_lpc2, st->lpcSize);
      else
      {
         st->bw_lpc2[0]=1;
         st->bw_lpc2[1]=-st->preemph;
         for (i=2;i<=st->lpcSize;i++)
            st->bw_lpc2[i]=0;
      }

      /* Reset excitation */
      for (i=0;i<st->subframeSize;i++)
         exc[i]=0;
      for (i=0;i<st->subframeSize;i++)
         exc2[i]=0;

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

      esig=0;
      for (i=0;i<st->subframeSize;i++)
         esig+=sw[i]*sw[i];
      
      /* Compute target signal */
      for (i=0;i<st->subframeSize;i++)
         target[i]=sw[i]-res[i];

      for (i=0;i<st->subframeSize;i++)
         exc[i]=exc2[i]=0;

      /* If we have a long-term predictor (not all sub-modes have one) */
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

         pitch = SUBMODE(ltp_quant)(target, sw, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2,
                                    exc, SUBMODE(ltp_params), pit_min, pit_max, ol_pitch_coef,
                                    st->lpcSize, st->subframeSize, bits, st->stack, 
                                    exc2, st->complexity);

         /*printf ("cl_pitch: %d\n", pitch);*/
         st->pitch[sub]=pitch;
      } else {
         fprintf (stderr, "No pitch prediction, what's wrong\n");
      }

      /* Update target for adaptive codebook contribution */
      residue_zero(exc, st->bw_lpc1, res, st->subframeSize, st->lpcSize);
      syn_filt_zero(res, st->interp_qlpc, res, st->subframeSize, st->lpcSize);
      syn_filt_zero(res, st->bw_lpc2, res, st->subframeSize, st->lpcSize);
      for (i=0;i<st->subframeSize;i++)
        target[i]-=res[i];

      /* Compute noise energy and SNR */
      enoise=0;
      for (i=0;i<st->subframeSize;i++)
         enoise += target[i]*target[i];
      snr = 10*log10((esig+1)/(enoise+1));
      /*st->pitch[sub]=(int)snr;*/
#ifdef DEBUG
      printf ("pitch SNR = %f\n", snr);
#endif


      /* Quantization of innovation */
      {
         float *innov;
         float ener=0, ener_1;
         /*innov=PUSH(st->stack, st->subframeSize);*/
         innov = st->innov+sub*st->subframeSize;
         for (i=0;i<st->subframeSize;i++)
            innov[i]=0;
         syn_filt_zero(target, st->bw_lpc1, res, st->subframeSize, st->lpcSize);
         residue_zero(res, st->interp_qlpc, st->buf2, st->subframeSize, st->lpcSize);
         residue_zero(st->buf2, st->bw_lpc2, st->buf2, st->subframeSize, st->lpcSize);
         for (i=0;i<st->subframeSize;i++)
            ener+=st->buf2[i]*st->buf2[i];
         ener=sqrt(.1+ener/st->subframeSize);

         
         ener /= ol_gain;

         if (0)
            printf ("ener: %f %f %f\n", ener, ol_gain, ol_pitch_coef);

         if (SUBMODE(have_subframe_gain)) 
         {
            int qe;
            ener=log(ener);
            if (SUBMODE(have_subframe_gain)==3)
            {
               qe = vq_index(&ener, exc_gain_quant_scal3, 1, 8);
               speex_bits_pack(bits, qe, 3);
               ener=exc_gain_quant_scal3[qe];
            } else {
               qe = vq_index(&ener, exc_gain_quant_scal1, 1, 2);
               speex_bits_pack(bits, qe, 1);
               ener=exc_gain_quant_scal1[qe];               
            }
            ener=exp(ener);
            /*printf ("encode gain: %d %f\n", qe, ener);*/
         } else {
            ener=1;
         }

         ener*=ol_gain;
         /*printf ("transmit gain: %f\n", ener);*/
         ener_1 = 1/ener;

         if (0) {
            int start=rand()%35;
            printf ("norm_exc: ");
            for (i=start;i<start+5;i++)
               printf ("%f ", ener_1*st->buf2[i]);
            printf ("\n");
         }
         
         for (i=0;i<st->subframeSize;i++)
            target[i]*=ener_1;
         
         if (SUBMODE(innovation_quant))
         {
            /* Normal quantization */
            SUBMODE(innovation_quant)(target, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2, 
                                      SUBMODE(innovation_params), st->lpcSize, st->subframeSize, 
                                      innov, bits, st->stack, st->complexity);
            for (i=0;i<st->subframeSize;i++)
               innov[i]*=ener;
            for (i=0;i<st->subframeSize;i++)
               exc[i] += innov[i];
         } else {
            fprintf(stderr, "No fixed codebook\n");
         }

         if (SUBMODE(double_codebook)) {
            float *innov2 = PUSH(st->stack, st->subframeSize);
            for (i=0;i<st->subframeSize;i++)
               innov2[i]=0;
            for (i=0;i<st->subframeSize;i++)
               target[i]*=2.2;
            SUBMODE(innovation_quant)(target, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2, 
                                      SUBMODE(innovation_params), st->lpcSize, st->subframeSize, 
                                      innov2, bits, st->stack, st->complexity);
            for (i=0;i<st->subframeSize;i++)
               innov2[i]*=ener*(1/2.2);
            for (i=0;i<st->subframeSize;i++)
               exc[i] += innov2[i];
            POP(st->stack);
         }

         /*POP(st->stack);*/
         for (i=0;i<st->subframeSize;i++)
            target[i]*=ener;

      }

      /* Compute weighted noise energy and SNR */
      enoise=0;
      for (i=0;i<st->subframeSize;i++)
         enoise += target[i]*target[i];
      snr = 10*log10((esig+1)/(enoise+1));
#ifdef DEBUG
      printf ("seg SNR = %f\n", snr);
#endif

      /*Keep the previous memory*/
      for (i=0;i<st->lpcSize;i++)
         mem[i]=st->mem_sp[i];
      /* Final signal synthesis from excitation */
      iir_mem2(exc, st->interp_qlpc, sp, st->subframeSize, st->lpcSize, st->mem_sp);

      /* Compute weighted signal again, from synthesized speech (not sure it's the right thing) */
      filter_mem2(sp, st->bw_lpc1, st->bw_lpc2, sw, st->subframeSize, st->lpcSize, st->mem_sw);
      for (i=0;i<st->subframeSize;i++)
         exc2[i]=exc[i];

      POP(st->stack);
      POP(st->stack);
      POP(st->stack);
   }

   /* Store the LSPs for interpolation in the next frame */
   for (i=0;i<st->lpcSize;i++)
      st->old_lsp[i] = st->lsp[i];
   for (i=0;i<st->lpcSize;i++)
      st->old_qlsp[i] = st->qlsp[i];

   /* The next frame will not be the first (Duh!) */
   st->first = 0;

   /* Replace input by synthesized speech */
   in[0] = st->frame[0] + st->preemph*st->pre_mem2;
   for (i=1;i<st->frameSize;i++)
     in[i]=st->frame[i] + st->preemph*in[i-1];
   st->pre_mem2=in[st->frameSize-1];

}


void *nb_decoder_init(SpeexMode *m)
{
   DecState *st;
   SpeexNBMode *mode;
   int i;

   mode=m->mode;
   st = speex_alloc(sizeof(DecState));
   st->mode=m;

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
   st->preemph = mode->preemph;

   st->submodes=mode->submodes;
   st->submodeID=mode->defaultSubmode;

   st->pre_mem=0;
   st->lpc_enh_enabled=0;

   st->stack = speex_alloc(20000*sizeof(float));

   st->inBuf = speex_alloc(st->bufSize*sizeof(float));
   st->frame = st->inBuf + st->bufSize - st->windowSize;
   st->excBuf = speex_alloc(st->bufSize*sizeof(float));
   st->exc = st->excBuf + st->bufSize - st->windowSize;
   for (i=0;i<st->bufSize;i++)
      st->inBuf[i]=0;
   for (i=0;i<st->bufSize;i++)
      st->excBuf[i]=0;
   st->innov = speex_alloc(st->frameSize*sizeof(float));

   st->interp_qlpc = speex_alloc((st->lpcSize+1)*sizeof(float));
   st->qlsp = speex_alloc(st->lpcSize*sizeof(float));
   st->old_qlsp = speex_alloc(st->lpcSize*sizeof(float));
   st->interp_qlsp = speex_alloc(st->lpcSize*sizeof(float));
   st->mem_sp = speex_alloc(5*st->lpcSize*sizeof(float));

   st->pi_gain = speex_alloc(st->nbSubframes*sizeof(float));
   st->last_pitch = 40;
   st->count_lost=0;


   st->user_callback.func = &speex_default_user_handler;
   st->user_callback.data = NULL;
   for (i=0;i<16;i++)
      st->speex_callbacks[i].func = NULL;

   return st;
}

void nb_decoder_destroy(void *state)
{
   DecState *st;
   st=state;
   speex_free(st->inBuf);
   speex_free(st->excBuf);
   speex_free(st->innov);
   speex_free(st->interp_qlpc);
   speex_free(st->qlsp);
   speex_free(st->old_qlsp);
   speex_free(st->interp_qlsp);
   speex_free(st->stack);
   speex_free(st->mem_sp);
   speex_free(st->pi_gain);
   
   speex_free(state);
}

static void nb_decode_lost(DecState *st, float *out)
{
   int i, sub;
   /*float exc_ener=0,g;*/
   /* Shift all buffers by one frame */
   speex_move(st->inBuf, st->inBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(float));
   speex_move(st->excBuf, st->excBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(float));

   for (sub=0;sub<st->nbSubframes;sub++)
   {
      int offset;
      float *sp, *exc;
      float *num, *den;
      /* Offset relative to start of frame */
      offset = st->subframeSize*sub;
      /* Original signal */
      sp=st->frame+offset;
      /* Excitation */
      exc=st->exc+offset;
      /* Excitation after post-filter*/

      num=PUSH(st->stack, ((st->lpcSize<<1)+1));
      den=PUSH(st->stack, ((st->lpcSize<<1)+1));
      if (st->lpc_enh_enabled)
      {
         enh_lpc(st->interp_qlpc, st->lpcSize, num, den, 
                 SUBMODE(lpc_enh_k1), SUBMODE(lpc_enh_k2), st->stack);
      } else {
         enh_lpc(st->interp_qlpc, st->lpcSize, num, den, 
                 SUBMODE(lpc_enh_k2), SUBMODE(lpc_enh_k2), st->stack);
      }
        
      for (i=0;i<st->subframeSize;i++)
      {
         exc[i]=st->last_pitch_gain*exc[i-st->last_pitch] + 
         .8*st->innov[i+offset];
      }

      for (i=0;i<st->subframeSize;i++)
         sp[i]=exc[i];
      
      /*pole_zero_mem(sp, num, den, sp, st->subframeSize, (st->lpcSize<<1), 
                    st->mem_sp+st->lpcSize, st->stack);*/
      filter_mem2(sp, num, den, sp, st->subframeSize, (st->lpcSize<<1), 
        st->mem_sp+st->lpcSize);
      iir_mem2(sp, st->interp_qlpc, sp, st->subframeSize, st->lpcSize, 
        st->mem_sp);
      
      POP(st->stack);
      POP(st->stack);
  
   }

   out[0] = st->frame[0] + st->preemph*st->pre_mem;
   for (i=1;i<st->frameSize;i++)
      out[i]=st->frame[i] + st->preemph*out[i-1];
   st->pre_mem=out[st->frameSize-1];
   
   st->first = 0;
   st->count_lost++;
}


int nb_decode(void *state, SpeexBits *bits, float *out)
{
   DecState *st;
   int i, sub;
   int pitch;
   float pitch_gain[3];
   float ol_gain;
   int ol_pitch=0;
   float ol_pitch_coef=0;
   int best_pitch=40;
   float best_pitch_gain=-1;
   int wideband;
   int m;

   st=state;

   if (!bits)
   {
      nb_decode_lost(st, out);
      return 0;
   }

   m = speex_bits_peek_unsigned(bits, 5);
   if (m==15)
   {
      speex_bits_unpack_unsigned(bits, 5);
      return -1;
   } else if (m==14)
   {
      int ret = speex_inband_handler(bits, st->speex_callbacks, state);
      if (ret)
         return ret;
   } else if (m==13)
   {
      int ret = st->user_callback.func(bits, state, st->user_callback.data);
      if (ret)
         return ret;
   } else if (m>7)
   {
      return -2;
   }

   wideband = speex_bits_unpack_unsigned(bits, 1);
   if (wideband)
   {
      int submode;
      int advance;
      submode = speex_bits_unpack_unsigned(bits, SB_SUBMODE_BITS);
      advance = submode;
      speex_mode_query(&speex_wb_mode, SPEEX_SUBMODE_BITS_PER_FRAME, &advance);
      advance -= (SB_SUBMODE_BITS+1);
      speex_bits_advance(bits, advance);
      wideband = speex_bits_unpack_unsigned(bits, 1);
      if (wideband)
      {
         fprintf (stderr, "Corrupted stream\n");
      }
   }

   /* Get the sub-mode that was used */
   st->submodeID = speex_bits_unpack_unsigned(bits, NB_SUBMODE_BITS);

   /* Shift all buffers by one frame */
   speex_move(st->inBuf, st->inBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(float));
   speex_move(st->excBuf, st->excBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(float));

   /* If null mode (no transmission), just set a couple things to zero*/
   if (st->submodes[st->submodeID] == NULL)
   {
      for (i=0;i<st->frameSize;i++)
         st->exc[i]=0;
      st->first=1;
      
      /* Final signal synthesis from excitation */
      iir_mem2(st->exc, st->interp_qlpc, st->frame, st->subframeSize, st->lpcSize, st->mem_sp);

      out[0] = st->frame[0] + st->preemph*st->pre_mem;
      for (i=1;i<st->frameSize;i++)
         out[i]=st->frame[i] + st->preemph*out[i-1];
      st->pre_mem=out[st->frameSize-1];
      st->count_lost=0;
      return 0;
   }

   /* Unquantize LSPs */
   SUBMODE(lsp_unquant)(st->qlsp, st->lpcSize, bits);

   /* Handle first frame and lost-packet case */
   if (st->first || st->count_lost)
   {
      for (i=0;i<st->lpcSize;i++)
         st->old_qlsp[i] = st->qlsp[i];
   }

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
      /*fprintf (stderr, "unquant pitch coef: %f\n", ol_pitch_coef);*/
   }
   
   /* Get global excitation gain */
   {
      int qe;
      qe = speex_bits_unpack_unsigned(bits, 5);
      ol_gain = exp(qe/3.5);
      /*printf ("decode_ol_gain: %f\n", ol_gain);*/
   }

   /*Loop on subframes */
   for (sub=0;sub<st->nbSubframes;sub++)
   {
      int offset;
      float *sp, *exc, tmp;
      float *num, *den;
      /* Offset relative to start of frame */
      offset = st->subframeSize*sub;
      /* Original signal */
      sp=st->frame+offset;
      /* Excitation */
      exc=st->exc+offset;
      /* Excitation after post-filter*/

      /* LSP interpolation (quantized and unquantized) */
      tmp = (1.0 + sub)/st->nbSubframes;
      for (i=0;i<st->lpcSize;i++)
         st->interp_qlsp[i] = (1-tmp)*st->old_qlsp[i] + tmp*st->qlsp[i];

      lsp_enforce_margin(st->interp_qlsp, st->lpcSize, .002);


      /* Compute interpolated LPCs (unquantized) */
      for (i=0;i<st->lpcSize;i++)
         st->interp_qlsp[i] = cos(st->interp_qlsp[i]);
      lsp_to_lpc(st->interp_qlsp, st->interp_qlpc, st->lpcSize, st->stack);

      num=PUSH(st->stack, ((st->lpcSize<<1)+1));
      den=PUSH(st->stack, ((st->lpcSize<<1)+1));
      if (st->lpc_enh_enabled)
      {
         enh_lpc(st->interp_qlpc, st->lpcSize, num, den, 
                 SUBMODE(lpc_enh_k1), SUBMODE(lpc_enh_k2), st->stack);
      } else {
         enh_lpc(st->interp_qlpc, st->lpcSize, num, den, 
                 SUBMODE(lpc_enh_k2), SUBMODE(lpc_enh_k2), st->stack);
      }
      /* Compute analysis filter at w=pi */
      tmp=1;
      st->pi_gain[sub]=0;
      for (i=0;i<=st->lpcSize;i++)
      {
         st->pi_gain[sub] += tmp*st->interp_qlpc[i];
         tmp = -tmp;
      }

      /* Reset excitation */
      for (i=0;i<st->subframeSize;i++)
         exc[i]=0;

      /*Adaptive codebook contribution*/
      if (SUBMODE(ltp_unquant))
      {
         int pit_min, pit_max;
         if (SUBMODE(lbr_pitch) != -1)
         {
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

         SUBMODE(ltp_unquant)(exc, pit_min, pit_max, ol_pitch_coef, SUBMODE(ltp_params), 
                              st->subframeSize, &pitch, &pitch_gain[0], bits, st->stack, st->count_lost);
         
         tmp = (pitch_gain[0]+pitch_gain[1]+pitch_gain[2]);
         if (tmp>best_pitch_gain)
         {
            best_pitch = pitch;
            /*while (best_pitch+pitch<st->max_pitch)
            {
               best_pitch+=pitch;
               }*/
            best_pitch_gain = tmp*.9;
            if (best_pitch_gain>.85)
               best_pitch_gain=.85;
         }
      } else {
         fprintf (stderr, "No pitch prediction, what's wrong\n");
      }
      
      /* Unquantize the innovation */
      {
         int q_energy;
         float ener;
         float *innov;
         
         /*innov = PUSH(st->stack, st->subframeSize);*/
         innov = st->innov+sub*st->subframeSize;
         for (i=0;i<st->subframeSize;i++)
            innov[i]=0;

         if (SUBMODE(have_subframe_gain)==3)
         {
            q_energy = speex_bits_unpack_unsigned(bits, 3);
            ener = ol_gain*exp(exc_gain_quant_scal3[q_energy]);
         } else if (SUBMODE(have_subframe_gain)==1)
         {
            q_energy = speex_bits_unpack_unsigned(bits, 1);
            ener = ol_gain*exp(exc_gain_quant_scal1[q_energy]);
         } else {
            ener = ol_gain;
         }
         
         /*printf ("unquant_energy: %d %f\n", q_energy, ener);*/
         
         if (SUBMODE(innovation_unquant))
         {
            /*Fixed codebook contribution*/
            SUBMODE(innovation_unquant)(innov, SUBMODE(innovation_params), st->subframeSize, bits, st->stack);
         } else {
            fprintf(stderr, "No fixed codebook\n");
         }

         for (i=0;i<st->subframeSize;i++)
            innov[i]*=ener;
         for (i=0;i<st->subframeSize;i++)
            exc[i]+=innov[i];

         if (SUBMODE(double_codebook))
         {
            float *innov2 = PUSH(st->stack, st->subframeSize);
            for (i=0;i<st->subframeSize;i++)
               innov2[i]=0;
            SUBMODE(innovation_unquant)(innov2, SUBMODE(innovation_params), st->subframeSize, bits, st->stack);
            for (i=0;i<st->subframeSize;i++)
               innov2[i]*=ener*(1/2.2);
            for (i=0;i<st->subframeSize;i++)
               exc[i] += innov2[i];
            POP(st->stack);
         }

         /*POP(st->stack);*/
      }

      for (i=0;i<st->subframeSize;i++)
         sp[i]=exc[i];

      if (st->lpc_enh_enabled && SUBMODE(comb_gain>0))
         comb_filter(exc, sp, st->interp_qlpc, st->lpcSize, st->subframeSize,
                              pitch, pitch_gain, .5);
      /*pole_zero_mem(sp, num, den, sp, st->subframeSize, (st->lpcSize<<1), 
        st->mem_sp+st->lpcSize, st->stack);*/
      filter_mem2(sp, num, den, sp, st->subframeSize, (st->lpcSize<<1), 
        st->mem_sp+st->lpcSize);
      iir_mem2(sp, st->interp_qlpc, sp, st->subframeSize, st->lpcSize, 
        st->mem_sp);
      
      POP(st->stack);
      POP(st->stack);
   }
   
   /*Copy output signal*/
   out[0] = st->frame[0] + st->preemph*st->pre_mem;
   for (i=1;i<st->frameSize;i++)
     out[i]=st->frame[i] + st->preemph*out[i-1];
   st->pre_mem=out[st->frameSize-1];


   /* Store the LSPs for interpolation in the next frame */
   for (i=0;i<st->lpcSize;i++)
      st->old_qlsp[i] = st->qlsp[i];

   /* The next frame will not be the first (Duh!) */
   st->first = 0;
   st->count_lost=0;
   st->last_pitch = best_pitch;
   st->last_pitch_gain = best_pitch_gain;

   return 0;
}

void nb_encoder_ctl(void *state, int request, void *ptr)
{
   EncState *st;
   st=state;     
   switch(request)
   {
   case SPEEX_GET_FRAME_SIZE:
      (*(int*)ptr) = st->frameSize;
      break;
   case SPEEX_SET_MODE:
      st->submodeID = (*(int*)ptr);
      break;
   case SPEEX_GET_MODE:
      (*(int*)ptr) = st->submodeID;
      break;
   case SPEEX_SET_VBR:
      st->vbr_enabled = (*(int*)ptr);
      break;
   case SPEEX_GET_VBR:
      (*(int*)ptr) = st->vbr_enabled;
      break;
   case SPEEX_SET_VBR_QUALITY:
      st->vbr_quality = (*(int*)ptr);
      break;
   case SPEEX_GET_VBR_QUALITY:
      (*(int*)ptr) = st->vbr_quality;
      break;
   case SPEEX_SET_QUALITY:
      {
         int quality = (*(int*)ptr);
         if (quality<=0)
            st->submodeID = 0;
         else if (quality<=1)
            st->submodeID = 1;
         else if (quality<=2)
            st->submodeID = 2;
         else if (quality<=4)
            st->submodeID = 3;
         else if (quality<=6)
            st->submodeID = 4;
         else if (quality<=8)
            st->submodeID = 5;
         else if (quality<=9)
            st->submodeID = 6;
         else if (quality<=10)
            st->submodeID = 7;
         else
            fprintf(stderr, "Unknown nb_ctl quality: %d\n", quality);
      }
      break;
   case SPEEX_SET_COMPLEXITY:
      st->complexity = (*(int*)ptr);
      break;
   case SPEEX_GET_COMPLEXITY:
      (*(int*)ptr) = st->complexity;
      break;
   case SPEEX_GET_BITRATE:
      if (st->submodes[st->submodeID])
         (*(int*)ptr) = 50*SUBMODE(bits_per_frame);
      else
         (*(int*)ptr) = 50*(NB_SUBMODE_BITS+1);
      break;
   default:
      fprintf(stderr, "Unknown nb_ctl request: %d\n", request);
   }
}

void nb_decoder_ctl(void *state, int request, void *ptr)
{
   DecState *st;
   st=state;
   switch(request)
   {
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
         (*(int*)ptr) = 50*SUBMODE(bits_per_frame);
      else
         (*(int*)ptr) = 50*(NB_SUBMODE_BITS+1);
      break;
   default:
      fprintf(stderr, "Unknown nb_ctl request: %d\n", request);
   }
}
