/* Copyright (C) 2002 Jean-Marc Valin 
   File: speex.c

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
#include <string.h>
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
#include "post_filter.h"
#include "vbr.h"

#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif

#define SUBMODE(x) st->submodes[st->submodeID]->x

float exc_gain_quant_scal[8]={-2.794750, -1.810660, -1.169850, -0.848119, -0.587190, -0.329818, -0.063266, 0.282826};

#define sqr(x) ((x)*(x))
#define min(a,b) ((a) < (b) ? (a) : (b))

void *nb_encoder_init(SpeexMode *m)
{
   EncState *st;
   SpeexNBMode *mode;
   int i;

   mode=m->mode;
   st = malloc(sizeof(EncState));
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
   st->inBuf = calloc(st->bufSize,sizeof(float));
   st->frame = st->inBuf + st->bufSize - st->windowSize;
   /* Allocating excitation buffer */
   st->excBuf = calloc(st->bufSize,sizeof(float));
   st->exc = st->excBuf + st->bufSize - st->windowSize;
   st->swBuf = calloc(st->bufSize,sizeof(float));
   st->sw = st->swBuf + st->bufSize - st->windowSize;

   st->exc2Buf = calloc(st->bufSize,sizeof(float));
   st->exc2 = st->exc2Buf + st->bufSize - st->windowSize;

   /* Asymetric "pseudo-Hamming" window */
   {
      int part1, part2;
      part1 = st->subframeSize*7/2;
      part2 = st->subframeSize*5/2;
      st->window = malloc(st->windowSize*sizeof(float));
      for (i=0;i<part1;i++)
         st->window[i]=.54-.46*cos(M_PI*i/part1);
      for (i=0;i<part2;i++)
         st->window[part1+i]=.54+.46*cos(M_PI*i/part2);
   }
   /* Create the window for autocorrelation (lag-windowing) */
   st->lagWindow = malloc((st->lpcSize+1)*sizeof(float));
   for (i=0;i<st->lpcSize+1;i++)
      st->lagWindow[i]=exp(-.5*sqr(2*M_PI*st->lag_factor*i));

   st->autocorr = malloc((st->lpcSize+1)*sizeof(float));

   st->stack = calloc(20000, sizeof(float));

   st->buf2 = malloc(st->windowSize*sizeof(float));

   st->lpc = malloc((st->lpcSize+1)*sizeof(float));
   st->interp_lpc = malloc((st->lpcSize+1)*sizeof(float));
   st->interp_qlpc = malloc((st->lpcSize+1)*sizeof(float));
   st->bw_lpc1 = malloc((st->lpcSize+1)*sizeof(float));
   st->bw_lpc2 = malloc((st->lpcSize+1)*sizeof(float));

   st->lsp = malloc(st->lpcSize*sizeof(float));
   st->qlsp = malloc(st->lpcSize*sizeof(float));
   st->old_lsp = malloc(st->lpcSize*sizeof(float));
   st->old_qlsp = malloc(st->lpcSize*sizeof(float));
   st->interp_lsp = malloc(st->lpcSize*sizeof(float));
   st->interp_qlsp = malloc(st->lpcSize*sizeof(float));
   st->rc = malloc(st->lpcSize*sizeof(float));
   st->first = 1;
   
   st->mem_sp = calloc(st->lpcSize, sizeof(float));
   st->mem_sw = calloc(st->lpcSize, sizeof(float));

   st->pi_gain = calloc(st->nbSubframes, sizeof(float));

   st->pitch = calloc(st->nbSubframes, sizeof(int));

   if (1) {
      st->vbr = malloc(sizeof(VBRState));
      vbr_init(st->vbr);
   } else {
      st->vbr = 0;
   }

   return st;
}

void nb_encoder_destroy(void *state)
{
   EncState *st=state;
   /* Free all allocated memory */
   free(st->inBuf);
   free(st->excBuf);
   free(st->swBuf);
   free(st->exc2Buf);
   free(st->stack);

   free(st->window);
   free(st->buf2);
   free(st->lpc);
   free(st->interp_lpc);
   free(st->interp_qlpc);
   
   free(st->bw_lpc1);
   free(st->bw_lpc2);
   free(st->autocorr);
   free(st->lagWindow);
   free(st->lsp);
   free(st->qlsp);
   free(st->old_lsp);
   free(st->interp_lsp);
   free(st->old_qlsp);
   free(st->interp_qlsp);
   free(st->rc);

   free(st->mem_sp);
   free(st->mem_sw);
   free(st->pi_gain);
   free(st->pitch);

   vbr_destroy(st->vbr);
   free(st->vbr);

   /*Free state memory... should be last*/
   free(st);
}

void nb_encode(void *state, float *in, SpeexBits *bits)
{
   EncState *st;
   int i, sub, roots;
   float error;
   int ol_pitch;
   float ol_gain;
   float vbr_qual=0;

   st=state;
   
   if (st->vbr)
      vbr_qual = vbr_analysis(st->vbr, in, st->frameSize);
   if (0) {
      int qual = (int)floor(8+1*vbr_qual+.5);
      if (qual<0)
         qual=0;
      if (qual>10)
         qual=10;
      speex_encoder_ctl(state, SPEEX_SET_QUALITY, &qual);
   }
   printf ("VBR quality = %f\n", vbr_qual);

   /* First, transmit the sub-mode we use for this frame */
   speex_bits_pack(bits, st->submodeID, NB_SUBMODE_BITS);

   /* Copy new data in input buffer */
   memmove(st->inBuf, st->inBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(float));
   st->inBuf[st->bufSize-st->frameSize] = in[0] - st->preemph*st->pre_mem;
   for (i=1;i<st->frameSize;i++)
      st->inBuf[st->bufSize-st->frameSize+i] = in[i] - st->preemph*in[i-1];
   st->pre_mem = in[st->frameSize-1];

   memmove(st->exc2Buf, st->exc2Buf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(float));
   memmove(st->excBuf, st->excBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(float));
   memmove(st->swBuf, st->swBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(float));

   /* Window for analysis */
   for (i=0;i<st->windowSize;i++)
      st->buf2[i] = st->frame[i] * st->window[i];

   /* Compute auto-correlation */
   autocorr(st->buf2, st->autocorr, st->lpcSize+1, st->windowSize);

   st->autocorr[0] += 1;        /* prevents NANs */
   st->autocorr[0] *= st->lpc_floor; /* Noise floor in auto-correlation domain */
   /* Lag windowing: equivalent to filtering in the power-spectrum domain */
   for (i=0;i<st->lpcSize+1;i++)
      st->autocorr[i] *= st->lagWindow[i];

   /* Levinson-Durbin */
   error = wld(st->lpc+1, st->autocorr, st->rc, st->lpcSize);
   st->lpc[0]=1;

   /* LPC to LSPs (x-domain) transform */
   roots=lpc_to_lsp (st->lpc, st->lpcSize, st->lsp, 6, 0.002, st->stack);
   if (roots!=st->lpcSize)
   {
      fprintf (stderr, "roots!=st->lpcSize (found only %d roots)\n", roots);
      exit(1);
   }

   /* x-domain to angle domain*/
   for (i=0;i<st->lpcSize;i++)
      st->lsp[i] = acos(st->lsp[i]);
   /*print_vec(st->lsp, 10, "LSP:");*/
   /* LSP Quantization */
#if 1
   SUBMODE(lsp_quant)(st->lsp, st->qlsp, st->lpcSize, bits);
#else
   for (i=0;i<st->lpcSize;i++)
     st->qlsp[i]=st->lsp[i];
#endif

   /* Special case for first frame */
   if (st->first)
   {
      for (i=0;i<st->lpcSize;i++)
         st->old_lsp[i] = st->lsp[i];
      for (i=0;i<st->lpcSize;i++)
         st->old_qlsp[i] = st->qlsp[i];
   }


   /* Whole frame analysis (some open-loop estimations) */
   {
      for (i=0;i<st->lpcSize;i++)
         st->interp_lsp[i] = .5*st->old_lsp[i] + .5*st->lsp[i];

      lsp_enforce_margin(st->interp_lsp, st->lpcSize, .002);

      /* Compute interpolated LPCs (unquantized) */
      for (i=0;i<st->lpcSize;i++)
         st->interp_lsp[i] = cos(st->interp_lsp[i]);
      lsp_to_lpc(st->interp_lsp, st->interp_lpc, st->lpcSize,st->stack);

      bw_lpc(st->gamma1, st->interp_lpc, st->bw_lpc1, st->lpcSize);
      bw_lpc(st->gamma2, st->interp_lpc, st->bw_lpc2, st->lpcSize);

      residue(st->frame, st->bw_lpc1, st->exc, st->frameSize, st->lpcSize);
      syn_filt(st->exc, st->bw_lpc2, st->sw, st->frameSize, st->lpcSize);
      
      if (SUBMODE(lbr_pitch) && SUBMODE(ltp_params))
      {
         open_loop_nbest_pitch(st->sw, st->min_pitch, st->max_pitch, st->frameSize, &ol_pitch, 1, st->stack);
         speex_bits_pack(bits, ol_pitch-st->min_pitch, 7);
      } else 
         ol_pitch = 0;

      residue(st->frame, st->interp_lpc, st->exc, st->frameSize, st->lpcSize);

      /* Compute open-loop excitation gain */
      ol_gain=0;
      for (i=0;i<st->frameSize;i++)
         ol_gain += st->exc[i]*st->exc[i];
      
      ol_gain=sqrt(1+ol_gain/st->frameSize);

      /* Quantize open-loop gain */
      /*printf ("ol_gain: %f\n", ol_gain);*/
      if (1) {
         int qe = (int)(floor(3.5*log(ol_gain)));
         if (qe<0)
            qe=0;
         if (qe>31)
            qe=31;
         ol_gain = exp(qe/3.5);
         speex_bits_pack(bits, qe, 5);
      }

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

      /* Compute zero response of A(z/g1) / ( A(z/g2) * Aq(z) ) */
      for (i=0;i<st->lpcSize;i++)
         mem[i]=st->mem_sp[i];
      syn_filt_mem(exc, st->interp_qlpc, exc, st->subframeSize, st->lpcSize, mem);
      for (i=0;i<st->lpcSize;i++)
         mem[i]=st->mem_sp[i];
      residue_mem(exc, st->bw_lpc1, res, st->subframeSize, st->lpcSize, mem);
      for (i=0;i<st->lpcSize;i++)
         mem[i]=st->mem_sw[i];
      syn_filt_mem(res, st->bw_lpc2, res, st->subframeSize, st->lpcSize, mem);

      /* Compute weighted signal */
      for (i=0;i<st->lpcSize;i++)
         mem[i]=st->mem_sp[i];
      residue_mem(sp, st->bw_lpc1, sw, st->subframeSize, st->lpcSize, mem);
      for (i=0;i<st->lpcSize;i++)
         mem[i]=st->mem_sw[i];
      syn_filt_mem(sw, st->bw_lpc2, sw, st->subframeSize, st->lpcSize, mem);
      
      esig=0;
      for (i=0;i<st->subframeSize;i++)
         esig+=sw[i]*sw[i];
      
      /* Compute target signal */
      for (i=0;i<st->subframeSize;i++)
         target[i]=sw[i]-res[i];

      for (i=0;i<st->subframeSize;i++)
         exc[i]=exc2[i]=0;

      /* If we have a long-term predictor (not all sub-modes have one) */
      if (SUBMODE(ltp_params))
      {
         /* Long-term prediction */
         if (SUBMODE(lbr_pitch) != -1)
         {
            /* Low bit-rate pitch handling */
            int pit_min, pit_max;
            int margin;
            margin = SUBMODE(lbr_pitch);
            if (ol_pitch < st->min_pitch+margin-1)
               ol_pitch=st->min_pitch+margin-1;
            if (ol_pitch > st->max_pitch-margin)
               ol_pitch=st->max_pitch-margin;
            pit_min = ol_pitch-margin+1;
            pit_max = ol_pitch+margin;
            pitch = SUBMODE(ltp_quant)(target, sw, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2,
                                       exc, SUBMODE(ltp_params), pit_min, pit_max, 
                                       st->lpcSize, st->subframeSize, bits, st->stack, exc2);
         } else {
            /* Normal pitch handling */
            pitch = SUBMODE(ltp_quant)(target, sw, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2,
                                       exc, SUBMODE(ltp_params), st->min_pitch, st->max_pitch, 
                                       st->lpcSize, st->subframeSize, bits, st->stack, exc2);
         }
         /*printf ("cl_pitch: %d\n", pitch);*/
         st->pitch[sub]=pitch;
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


#if 0 /*If set to 1, compute "real innovation" i.e. cheat to get perfect reconstruction*/
      syn_filt_zero(target, st->bw_lpc1, res, st->subframeSize, st->lpcSize);
      residue_zero(res, st->interp_qlpc, st->buf2, st->subframeSize, st->lpcSize);
      residue_zero(st->buf2, st->bw_lpc2, st->buf2, st->subframeSize, st->lpcSize);
      /*if (1||(snr>9 && (rand()%6==0)))
      {
         float ener=0;
         printf ("exc ");
         for (i=0;i<st->subframeSize;i++)
         {
            ener+=st->buf2[i]*st->buf2[i];
            if (i && i%5==0)
               printf ("\nexc ");
            printf ("%f ", st->buf2[i]);
         }
         printf ("\n");
      printf ("innovation_energy = %f\n", ener);
      }*/
      if (rand()%5==0 && snr>5)
      {
         float ener=0, sign=1;
         if (rand()%2)
            sign=-1;
         for (i=0;i<st->subframeSize;i++)
         {
            ener+=st->buf2[i]*st->buf2[i];
         }
         ener=sign/sqrt(.01+ener/st->subframeSize);
         for (i=0;i<st->subframeSize;i++)
         {
            if (i%10==0)
               printf ("\nexc ");
            printf ("%f ", ener*st->buf2[i]);
         }
         printf ("\n");
      }

      for (i=0;i<st->subframeSize;i++)
         exc[i]+=st->buf2[i];
#else
      /* Quantization of innovation */
      {
         float *innov;
         float ener=0, ener_1;
         innov=PUSH(st->stack, st->subframeSize);
         for (i=0;i<st->subframeSize;i++)
            innov[i]=0;
         syn_filt_zero(target, st->bw_lpc1, res, st->subframeSize, st->lpcSize);
         residue_zero(res, st->interp_qlpc, st->buf2, st->subframeSize, st->lpcSize);
         residue_zero(st->buf2, st->bw_lpc2, st->buf2, st->subframeSize, st->lpcSize);
         for (i=0;i<st->subframeSize;i++)
            ener+=st->buf2[i]*st->buf2[i];
         ener=sqrt(.1+ener/st->subframeSize);

         ener /= ol_gain;
         if (SUBMODE(have_subframe_gain)) 
         {
            int qe;
            ener=log(ener);
            qe = vq_index(&ener, exc_gain_quant_scal, 1, 8);
            speex_bits_pack(bits, qe, 3);
            ener=exc_gain_quant_scal[qe];
            ener=exp(ener);
            /*printf ("encode gain: %d %f\n", qe, ener);*/
         } else {
            ener=1;
         }
         ener*=ol_gain;
         /*printf ("transmit gain: %f\n", ener);*/
         ener_1 = 1/ener;
         
         for (i=0;i<st->subframeSize;i++)
            target[i]*=ener_1;
         
         if (SUBMODE(innovation_quant))
         {
            /* Normal quantization */
            SUBMODE(innovation_quant)(target, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2, 
                                      SUBMODE(innovation_params), st->lpcSize, st->subframeSize, 
                                      innov, bits, st->stack);
            
            for (i=0;i<st->subframeSize;i++)
               exc[i] += innov[i]*ener;
         } else {
            /* This is the "real" (cheating) excitation in the encoder but the decoder will
               use white noise */
            for (i=0;i<st->subframeSize;i++)
               exc[i] += st->buf2[i];
         }
         POP(st->stack);
         for (i=0;i<st->subframeSize;i++)
            target[i]*=ener;

      }
#endif
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
      syn_filt_mem(exc, st->interp_qlpc, sp, st->subframeSize, st->lpcSize, st->mem_sp);

      /* Compute weighted signal again, from synthesized speech (not sure it's the right thing) */
      residue_mem(sp, st->bw_lpc1, sw, st->subframeSize, st->lpcSize, mem);
      syn_filt_mem(sw, st->bw_lpc2, sw, st->subframeSize, st->lpcSize, st->mem_sw);

#if 0
      /*for (i=0;i<st->subframeSize;i++)
        exc2[i]=.75*exc[i]+.2*exc[i-pitch]+.05*exc[i-2*pitch];*/
      {
         float max_exc=0;
         for (i=0;i<st->subframeSize;i++)
            if (fabs(exc[i])>max_exc)
               max_exc=fabs(exc[i]);
         max_exc=1/(max_exc+.01);
         for (i=0;i<st->subframeSize;i++)
         {
            float xx=max_exc*exc[i];
            exc2[i]=exc[i]*(1-exp(-100*xx*xx));
         }
      }
#else
      for (i=0;i<st->subframeSize;i++)
         exc2[i]=exc[i];
#endif
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
   st = malloc(sizeof(DecState));
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
   st->pf_enabled=0;

   st->stack = calloc(10000, sizeof(float));

   st->inBuf = malloc(st->bufSize*sizeof(float));
   st->frame = st->inBuf + st->bufSize - st->windowSize;
   st->excBuf = malloc(st->bufSize*sizeof(float));
   st->exc = st->excBuf + st->bufSize - st->windowSize;
   st->exc2Buf = malloc(st->bufSize*sizeof(float));
   st->exc2 = st->exc2Buf + st->bufSize - st->windowSize;
   for (i=0;i<st->bufSize;i++)
      st->inBuf[i]=0;
   for (i=0;i<st->bufSize;i++)
      st->excBuf[i]=0;
   for (i=0;i<st->bufSize;i++)
      st->exc2Buf[i]=0;

   st->interp_qlpc = malloc((st->lpcSize+1)*sizeof(float));
   st->qlsp = malloc(st->lpcSize*sizeof(float));
   st->old_qlsp = malloc(st->lpcSize*sizeof(float));
   st->interp_qlsp = malloc(st->lpcSize*sizeof(float));
   st->mem_sp = calloc(st->lpcSize, sizeof(float));
   st->mem_pf = calloc(st->lpcSize, sizeof(float));
   st->mem_pf2 = calloc(st->lpcSize, sizeof(float));

   st->pi_gain = calloc(st->nbSubframes, sizeof(float));
   st->last_pitch = 40;
   st->count_lost=0;
   return st;
}

void nb_decoder_destroy(void *state)
{
   DecState *st;
   st=state;
   free(st->inBuf);
   free(st->excBuf);
   free(st->exc2Buf);
   free(st->interp_qlpc);
   free(st->qlsp);
   free(st->old_qlsp);
   free(st->interp_qlsp);
   free(st->stack);
   free(st->mem_sp);
   free(st->mem_pf);
   free(st->mem_pf2);
   free(st->pi_gain);
   
   free(state);
}

void nb_decode(void *state, SpeexBits *bits, float *out, int lost)
{
   DecState *st;
   int i, sub;
   int pitch;
   float pitch_gain[3];
   float ol_gain;
   int ol_pitch=0;
   int best_pitch=40;
   float best_pitch_gain=-1;
   st=state;

   /* Get the sub-mode that was used */
   st->submodeID = speex_bits_unpack_unsigned(bits, NB_SUBMODE_BITS);

   /* Shift all buffers by one frame */
   memmove(st->inBuf, st->inBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(float));
   memmove(st->excBuf, st->excBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(float));
   memmove(st->exc2Buf, st->exc2Buf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(float));

   /* Unquantize LSPs */
   SUBMODE(lsp_unquant)(st->qlsp, st->lpcSize, bits);

   /* Handle first frame and lost-packet case */
   if (st->first || st->count_lost)
   {
      for (i=0;i<st->lpcSize;i++)
         st->old_qlsp[i] = st->qlsp[i];
   }

   /* Get open-loop pitch estimation for low bit-rate pitch coding */
   if (SUBMODE(lbr_pitch) && SUBMODE(ltp_params))
      ol_pitch = st->min_pitch+speex_bits_unpack_unsigned(bits, 7);
   
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
      float *sp, *exc, *exc2, tmp;
      
      /* Offset relative to start of frame */
      offset = st->subframeSize*sub;
      /* Original signal */
      sp=st->frame+offset;
      /* Excitation */
      exc=st->exc+offset;
      /* Excitation after post-filter*/
      exc2=st->exc2+offset;

      /* LSP interpolation (quantized and unquantized) */
      tmp = (1.0 + sub)/st->nbSubframes;
      for (i=0;i<st->lpcSize;i++)
         st->interp_qlsp[i] = (1-tmp)*st->old_qlsp[i] + tmp*st->qlsp[i];

      lsp_enforce_margin(st->interp_qlsp, st->lpcSize, .002);


      /* Compute interpolated LPCs (unquantized) */
      for (i=0;i<st->lpcSize;i++)
         st->interp_qlsp[i] = cos(st->interp_qlsp[i]);
      lsp_to_lpc(st->interp_qlsp, st->interp_qlpc, st->lpcSize, st->stack);


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
         if (SUBMODE(lbr_pitch) != -1)
         {
            int pit_min, pit_max;
            int margin;
            margin = SUBMODE(lbr_pitch);
            if (ol_pitch < st->min_pitch+margin-1)
               ol_pitch=st->min_pitch+margin-1;
            if (ol_pitch > st->max_pitch-margin)
               ol_pitch=st->max_pitch-margin;
            pit_min = ol_pitch-margin+1;
            pit_max = ol_pitch+margin;
            SUBMODE(ltp_unquant)(exc, pit_min, pit_max, SUBMODE(ltp_params), st->subframeSize, &pitch, &pitch_gain[0], bits, st->stack, 0);
         } else {
            SUBMODE(ltp_unquant)(exc, st->min_pitch, st->max_pitch, SUBMODE(ltp_params), st->subframeSize, &pitch, &pitch_gain[0], bits, st->stack, 0);
         }
         
         if (!lost)
         {
            /* If the frame was not lost... */
            tmp = fabs(pitch_gain[0])+fabs(pitch_gain[1])+fabs(pitch_gain[2]);
            tmp = fabs(pitch_gain[0]+pitch_gain[1]+pitch_gain[2]);
            if (tmp>best_pitch_gain)
            {
               best_pitch = pitch;
               while (best_pitch+pitch<st->max_pitch)
               {
                  best_pitch+=pitch;
               }
               best_pitch_gain = tmp*.9;
               if (best_pitch_gain>.85)
                  best_pitch_gain=.85;
            }
         } else {
            /* What to do with pitch if we lost the frame */
            for (i=0;i<st->subframeSize;i++)
               exc[i]=0;
            /*printf ("best_pitch: %d %f\n", st->last_pitch, st->last_pitch_gain);*/
            for (i=0;i<st->subframeSize;i++)
               exc[i]=st->last_pitch_gain*exc[i-st->last_pitch];
         }
      }
      
      /* Unquantize the innovation */
      {
         int q_energy;
         float ener;
         float *innov;
         
         innov = PUSH(st->stack, st->subframeSize);
         for (i=0;i<st->subframeSize;i++)
            innov[i]=0;

         if (SUBMODE(have_subframe_gain))
         {
            q_energy = speex_bits_unpack_unsigned(bits, 3);
            ener = ol_gain*exp(exc_gain_quant_scal[q_energy]);
         } else {
            ener = ol_gain;
         }
         
         /*printf ("unquant_energy: %d %f\n", q_energy, ener);*/
         
         if (SUBMODE(innovation_unquant))
         {
            /*Fixed codebook contribution*/
            SUBMODE(innovation_unquant)(innov, SUBMODE(innovation_params), st->subframeSize, bits, st->stack);
         } else {
            for (i=0;i<st->subframeSize;i++)
               innov[i] = 3*((((float)rand())/RAND_MAX)-.5);
            
         }

         if (st->count_lost)
            ener*=pow(.8,st->count_lost);

         for (i=0;i<st->subframeSize;i++)
            exc[i]+=ener*innov[i];

         POP(st->stack);
      }

      for (i=0;i<st->subframeSize;i++)
         exc2[i]=exc[i];

      /* Apply post-filter */
      if (st->pf_enabled && SUBMODE(post_filter_func))
         SUBMODE(post_filter_func)(exc, exc2, st->interp_qlpc, st->lpcSize, st->subframeSize,
                              pitch, pitch_gain, SUBMODE(post_filter_params), st->mem_pf, 
                              st->mem_pf2, st->stack);
      
      /* Apply synthesis filter */
      syn_filt_mem(exc2, st->interp_qlpc, sp, st->subframeSize, st->lpcSize, st->mem_sp);

   }
   
   /*Copy output signal*/
   for (i=0;i<st->frameSize;i++)
      out[i]=st->frame[i];

   out[0] = st->frame[0] + st->preemph*st->pre_mem;
   for (i=1;i<st->frameSize;i++)
     out[i]=st->frame[i] + st->preemph*out[i-1];
   st->pre_mem=out[st->frameSize-1];


   /* Store the LSPs for interpolation in the next frame */
   for (i=0;i<st->lpcSize;i++)
      st->old_qlsp[i] = st->qlsp[i];

   /* The next frame will not be the first (Duh!) */
   st->first = 0;
   if (!lost)
      st->count_lost=0;
   else
      st->count_lost++;
   if (!lost)
   {
      st->last_pitch = best_pitch;
      st->last_pitch_gain = best_pitch_gain;
   }
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
   case SPEEX_SET_QUALITY:
      {
         int quality = (*(int*)ptr);
         if (quality<=0)
            st->submodeID = 1;
         else if (quality<=2)
            st->submodeID = 1;
         else if (quality<=4)
            st->submodeID = 2;
         else if (quality<=6)
            st->submodeID = 3;
         else if (quality<=8)
            st->submodeID = 4;
         else if (quality<=10)
            st->submodeID = 5;
         else
            fprintf(stderr, "Unknown nb_ctl quality: %d\n", quality);
      }
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
   case SPEEX_SET_PF:
      st->pf_enabled = *((int*)ptr);
      break;
   case SPEEX_GET_FRAME_SIZE:
      (*(int*)ptr) = st->frameSize;
      break;
   default:
      fprintf(stderr, "Unknown nb_ctl request: %d\n", request);
   }
}
