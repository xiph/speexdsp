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
#include "speex.h"
#include "lpc.h"
#include "lsp.h"
#include "ltp.h"
#include "quant_lsp.h"
#include "cb_search.h"
#include "filters.h"
#include "stack_alloc.h"

extern float stoc[];
extern float exc_table[][8];
#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif

#define sqr(x) ((x)*(x))
#define min(a,b) ((a) < (b) ? (a) : (b))

void encoder_init(EncState *st, SpeexMode *mode)
{
   int i;
   float tmp;
   /* Codec parameters, should eventually have several "modes"*/
   st->frameSize = mode->frameSize;
   st->windowSize = mode->windowSize;
   st->nbSubframes=mode->frameSize/mode->subframeSize;
   st->subframeSize=mode->subframeSize;
   st->lpcSize = mode->lpcSize;
   st->bufSize = mode->bufSize;
   st->gamma1=mode->gamma1;
   st->gamma2=mode->gamma2;
   st->min_pitch=mode->pitchStart;
   st->max_pitch=mode->pitchEnd;

   /* Over-sampling filter (fractional pitch)*/
   st->os_fact=4;
   st->os_filt_ord2=4*st->os_fact;
   st->os_filt = malloc((1+2*st->os_filt_ord2)*sizeof(float));
   st->os_filt[st->os_filt_ord2] = 1;
   for (i=1;i<=st->os_filt_ord2;i++)
   {
      float x=M_PI*i/st->os_fact;
      st->os_filt[st->os_filt_ord2-i] = st->os_filt[st->os_filt_ord2+i]=sin(x)/x*(.5+.5*cos(M_PI*i/st->os_filt_ord2));
   }
   /* Normalizing the over-sampling filter */
   tmp=0;
   for (i=0;i<2*st->os_filt_ord2+1;i++)
      tmp += st->os_filt[i];
   tmp=1/tmp;
   for (i=0;i<2*st->os_filt_ord2+1;i++)
      st->os_filt[i] *= tmp;

   /*for (i=0;i<2*st->os_filt_ord2+1;i++)
      printf ("%f ", st->os_filt[i]);
      printf ("\n");*/

   /* Allocating input buffer */
   st->inBuf = calloc(st->bufSize,sizeof(float));
   st->frame = st->inBuf + st->bufSize - st->windowSize;
   /* Allocating excitation buffer */
   st->excBuf = calloc(st->bufSize,sizeof(float));
   st->exc = st->excBuf + st->bufSize - st->windowSize;
   st->swBuf = calloc(st->bufSize,sizeof(float));
   st->sw = st->swBuf + st->bufSize - st->windowSize;

   /* Hanning window */
   st->window = malloc(st->windowSize*sizeof(float));
   for (i=0;i<st->windowSize;i++)
      st->window[i]=.5*(1-cos(2*M_PI*i/st->windowSize));

   /* Create the window for autocorrelation (lag-windowing) */
   st->lagWindow = malloc((st->lpcSize+1)*sizeof(float));
   for (i=0;i<st->lpcSize+1;i++)
      st->lagWindow[i]=exp(-.5*sqr(2*M_PI*.01*i));

   st->autocorr = malloc((st->lpcSize+1)*sizeof(float));

   st->stack = calloc(10000, sizeof(float));

   st->buf2 = malloc(st->windowSize*sizeof(float));

   st->lpc = malloc((st->lpcSize+1)*sizeof(float));
   st->interp_lpc = malloc((st->lpcSize+1)*sizeof(float));
   st->interp_qlpc = malloc((st->lpcSize+1)*sizeof(float));
   st->bw_lpc1 = malloc((st->lpcSize+1)*sizeof(float));
   st->bw_lpc2 = malloc((st->lpcSize+1)*sizeof(float));
   st->bw_az = malloc((st->lpcSize*2+1)*sizeof(float));

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
}

void encoder_destroy(EncState *st)
{
   /* Free all allocated memory */
   free(st->inBuf);
   free(st->excBuf);
   free(st->swBuf);
   
   free(st->stack);

   free(st->window);
   free(st->buf2);
   free(st->lpc);
   free(st->interp_lpc);
   free(st->interp_qlpc);

   free(st->bw_lpc1);
   free(st->bw_lpc2);
   free(st->bw_az);
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
}

void encode(EncState *st, float *in, FrameBits *bits)
{
   int i, j, sub, roots;
   float error;

   /* Copy new data in input buffer */
   memmove(st->inBuf, st->inBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(float));
   for (i=0;i<st->frameSize;i++)
      st->inBuf[st->bufSize-st->frameSize+i] = in[i];
   memmove(st->excBuf, st->excBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(float));
   memmove(st->swBuf, st->swBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(float));

   /* Window for analysis */
   for (i=0;i<st->windowSize;i++)
      st->buf2[i] = st->frame[i] * st->window[i];

   /* Compute auto-correlation */
   autocorr(st->buf2, st->autocorr, st->lpcSize+1, st->windowSize);

   st->autocorr[0] += 1;        /* prevents NANs */
   st->autocorr[0] *= 1.0001;   /* 40 dB noise floor */
   /* Lag windowing: equivalent to filtering in the power-spectrum domain */
   for (i=0;i<st->lpcSize+1;i++)
      st->autocorr[i] *= st->lagWindow[i];

   /* Levinson-Durbin */
   error = wld(st->lpc+1, st->autocorr, st->rc, st->lpcSize);
   st->lpc[0]=1;

   /* LPC to LSPs (x-domain) transform */
   roots=lpc_to_lsp (st->lpc, st->lpcSize, st->lsp, 6, 0.02, st->stack);
   if (roots!=st->lpcSize)
   {
      fprintf (stderr, "roots!=st->lpcSize\n");
      exit(1);
   }

   /* x-domain to angle domain*/
   for (i=0;i<st->lpcSize;i++)
      st->lsp[i] = acos(st->lsp[i]);
   
   /* LSP Quantization */
   {
      unsigned int id;
      for (i=0;i<st->lpcSize;i++)
         st->qlsp[i]=st->lsp[i];
      id=lsp_quant_nb(st->qlsp,10 );
      lsp_unquant_nb(st->qlsp,10,id);
   }

   /*Find open-loop pitch for the whole frame*/
   {
      float *mem = PUSH(st->stack, st->lpcSize);
      
      for (i=0;i<st->lpcSize;i++)
         st->interp_lsp[i] = .5*st->old_lsp[i] + .5*st->lsp[i];
      for (i=0;i<st->lpcSize;i++)
         st->interp_lsp[i] = cos(st->interp_lsp[i]);
      lsp_to_lpc(st->interp_lsp, st->interp_lpc, st->lpcSize,st->stack);
      bw_lpc(st->gamma1, st->interp_lpc, st->bw_lpc1, st->lpcSize);
      bw_lpc(st->gamma2, st->interp_lpc, st->bw_lpc2, st->lpcSize);
      for (i=0;i<st->lpcSize;i++)
         mem[i]=st->mem_sp[i];
      residue_mem(st->frame, st->bw_lpc1, st->sw, st->frameSize, st->lpcSize, mem);
      for (i=0;i<st->lpcSize;i++)
         mem[i]=st->mem_sw[i];
         syn_filt_mem(st->sw, st->bw_lpc2, st->sw, st->frameSize, st->lpcSize, mem);
      open_loop_pitch(st->sw, st->min_pitch, st->max_pitch, st->frameSize, &st->ol_pitch, &st->ol_voiced);
      printf ("Open-loop pitch = %d\n", st->ol_pitch);
      POP(st->stack);
   }

   /* Loop on sub-frames */
   for (sub=0;sub<st->nbSubframes;sub++)
   {
      float tmp, tmp1,tmp2,gain[3];
      float esig=0, enoise=0, snr;
      int pitch, offset, pitch_gain_index;
      float *sp, *sw, *res, *exc, *target, *mem;
      
      /* Offset relative to start of frame */
      offset = st->subframeSize*sub;
      /* Original signal */
      sp=st->frame+offset;
      /* Excitation */
      exc=st->exc+offset;
      /* Weighted signal */
      sw=st->sw+offset;
      /* Filter response */
      res = PUSH(st->stack, st->subframeSize);
      /* Target signal */
      target = PUSH(st->stack, st->subframeSize);
      mem = PUSH(st->stack, st->lpcSize);

      /* LSP interpolation (quantized and unquantized) */
      tmp = (.5 + sub)/st->nbSubframes;
      for (i=0;i<st->lpcSize;i++)
         st->interp_lsp[i] = (1-tmp)*st->old_lsp[i] + tmp*st->lsp[i];
      for (i=0;i<st->lpcSize;i++)
         st->interp_qlsp[i] = (1-tmp)*st->old_qlsp[i] + tmp*st->qlsp[i];

      /* Compute interpolated LPCs (quantized and unquantized) */
      for (i=0;i<st->lpcSize;i++)
         st->interp_lsp[i] = cos(st->interp_lsp[i]);
      lsp_to_lpc(st->interp_lsp, st->interp_lpc, st->lpcSize,st->stack);

      for (i=0;i<st->lpcSize;i++)
         st->interp_qlsp[i] = cos(st->interp_qlsp[i]);
      lsp_to_lpc(st->interp_qlsp, st->interp_qlpc, st->lpcSize, st->stack);

      /* Compute bandwidth-expanded (unquantized) LPCs for perceptual weighting */
      bw_lpc(st->gamma1, st->interp_lpc, st->bw_lpc1, st->lpcSize);
      bw_lpc(st->gamma2, st->interp_lpc, st->bw_lpc2, st->lpcSize);
      
      /* Reset excitation */
      for (i=0;i<st->subframeSize;i++)
         exc[i]=0;

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

      for (i=0;i<st->subframeSize;i++)
         esig+=sw[i]*sw[i];
      
      /* Compute target signal */
      for (i=0;i<st->subframeSize;i++)
         target[i]=sw[i]-res[i];

      for (i=0;i<st->subframeSize;i++)
         exc[i]=0;
#if 1 /*If set to 0, we compute the excitation directly from the target, i.e. we're cheating */

      /* Perform adaptive codebook search (3-tap pitch predictor) */
      pitch = st->ol_pitch;
#if 0 /* 1 for fractional pitch, 0 for integer pitch */
      closed_loop_fractional_pitch(target, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2,
                                   exc, st->os_filt, st->os_filt_ord2, st->os_fact, 20, 147, 
                                   &gain[0], &pitch, st->lpcSize,
                                   st->subframeSize, st->stack);
#else
      pitch_search_3tap(target, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2,
                        exc, 20, 147, &gain[0], &pitch, &pitch_gain_index, st->lpcSize,
                        st->subframeSize);
      for (i=0;i<st->subframeSize;i++)
        exc[i]=gain[0]*exc[i-pitch]+gain[1]*exc[i-pitch-1]+gain[2]*exc[i-pitch-2];
      printf ("3-tap pitch = %d, gains = [%f %f %f]\n",pitch, gain[0], gain[1], gain[2]);
#endif
      /* Update target for adaptive codebook contribution */
      residue_zero(exc, st->bw_lpc1, res, st->subframeSize, st->lpcSize);
      syn_filt_zero(res, st->interp_qlpc, res, st->subframeSize, st->lpcSize);
      syn_filt_zero(res, st->bw_lpc2, res, st->subframeSize, st->lpcSize);
      for (i=0;i<st->subframeSize;i++)
        target[i]-=res[i];

      enoise=0;
      for (i=0;i<st->subframeSize;i++)
         enoise += target[i]*target[i];
      snr = 10*log10((esig+1)/(enoise+1));
      printf ("pitch SNR = %f\n", snr);
#if 0 /* 1 for stochastic excitation, 0 for split-VQ*/
      for(j=0;j<1;j++){
         /*float stoc2[1080];*/
         float *stoc2 = PUSH(st->stack,1080);
         for (i=0;i<1080;i++)
         {
            stoc2[i]=stoc[i];
            if (i-(pitch-1)>=0)
               stoc2[i] += .0*stoc[i-(pitch-1)];
         }
         POP(st->stack);
      /* Perform stochastic codebook search */
      overlap_cb_search(target, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2,
                        stoc2, 1024, &gain[0], &pitch, st->lpcSize,
                        st->subframeSize);
      printf ("gain = %f index = %d energy = %f\n",gain[0], pitch, esig);
      for (i=0;i<st->subframeSize;i++)
         exc[i]+=gain[0]*stoc2[i+pitch];
      
      /* Update target for adaptive codebook contribution (Useless for now)*/
      residue_zero(stoc2+pitch, st->bw_lpc1, res, st->subframeSize, st->lpcSize);
      syn_filt_zero(res, st->interp_qlpc, res, st->subframeSize, st->lpcSize);
      syn_filt_zero(res, st->bw_lpc2, res, st->subframeSize, st->lpcSize);
      for (i=0;i<st->subframeSize;i++)
         target[i]-=gain[0]*res[i];
      }
#else
      split_cb_search(target, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2,
                        exc_table, 64, &gain[0], &pitch, st->lpcSize,
                        st->subframeSize, exc);
#endif
      /* Compute weighted noise energy, SNR */
      enoise=0;
      for (i=0;i<st->subframeSize;i++)
         enoise += target[i]*target[i];
      snr = 10*log10((esig+1)/(enoise+1));
      printf ("seg SNR = %f\n", snr);

#else

#if 1 /* Code to calculate the exact excitation after pitch prediction  */
      for (i=0;i<st->subframeSize;i++)
         st->buf2[i]=target[i];
#if 0 /* 0 for fractional pitch, 1 for integer */
      pitch_search_3tap(target, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2,
                                exc, 20, 147, &gain[0], &pitch, &pitch_gain_index, st->lpcSize,
                        st->subframeSize);
      for (i=0;i<st->subframeSize;i++)
        exc[i]=gain[0]*exc[i-pitch]+gain[1]*exc[i-pitch-1]+gain[2]*exc[i-pitch-2];
      printf ("3-tap pitch = %d, gains = [%f %f %f]\n",pitch, gain[0], gain[1], gain[2]);
#else
      pitch = st->ol_pitch;
      closed_loop_fractional_pitch(target, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2,
                                   exc, st->os_filt, st->os_filt_ord2, st->os_fact, 20, 147, 
                                   &gain[0], &pitch, st->lpcSize,
                                   st->subframeSize, st->stack);
#endif
      /* Update target for adaptive codebook contribution */
      residue_zero(exc, st->bw_lpc1, res, st->subframeSize, st->lpcSize);
      syn_filt_zero(res, st->interp_qlpc, res, st->subframeSize, st->lpcSize);
      syn_filt_zero(res, st->bw_lpc2, res, st->subframeSize, st->lpcSize);
      for (i=0;i<st->subframeSize;i++)
        target[i]-=res[i];

      enoise=0;
      for (i=0;i<st->subframeSize;i++)
         enoise += target[i]*target[i];
      snr = 10*log10((esig+1)/(enoise+1));
      printf ("pitch SNR = %f\n", snr);

      syn_filt_zero(target, st->bw_lpc1, res, st->subframeSize, st->lpcSize);
      residue_zero(res, st->interp_qlpc, exc, st->subframeSize, st->lpcSize);
      residue_zero(exc, st->bw_lpc2, exc, st->subframeSize, st->lpcSize);
      if (snr>5)
      {
         for (i=0;i<st->subframeSize;i++)
         {
            if (i%8==0&&i)
               printf("\n");
            printf ("%f ", exc[i]);
         }
         printf ("\n");
      }
      for (i=0;i<st->subframeSize;i++)
         target[i]=st->buf2[i];
#endif

      /* We're cheating to get perfect reconstruction */
      syn_filt_zero(target, st->bw_lpc1, res, st->subframeSize, st->lpcSize);
      residue_zero(res, st->interp_qlpc, exc, st->subframeSize, st->lpcSize);
      residue_zero(exc, st->bw_lpc2, exc, st->subframeSize, st->lpcSize);
#endif

      /*Keep the previous memory*/
      for (i=0;i<st->lpcSize;i++)
         mem[i]=st->mem_sp[i];
      /* Final signal synthesis from excitation */
      syn_filt_mem(exc, st->interp_qlpc, sp, st->subframeSize, st->lpcSize, st->mem_sp);

      /* Compute weighted signal again, from synthesized speech (not sure it's the right thing) */
      residue_mem(sp, st->bw_lpc1, sw, st->subframeSize, st->lpcSize, mem);
      syn_filt_mem(sw, st->bw_lpc2, sw, st->subframeSize, st->lpcSize, st->mem_sw);

      POP(st->stack);
      POP(st->stack);
      POP(st->stack);
   }

   /* Store the LSPs for interpolation in the next frame */
   for (i=0;i<st->lpcSize;i++)
      st->old_lsp[i] = st->lsp[i];
   for (i=0;i<st->lpcSize;i++)
      st->old_qlsp[i] = st->qlsp[i];

   /* The next frame will not by the first (Duh!) */
   st->first = 0;

   /* Replace input by synthesized speech */
   for (i=0;i<st->frameSize;i++)
     in[i]=st->frame[i];
}


void decoder_init(DecState *st, SpeexMode *mode)
{
   int i;
   /* Codec parameters, should eventually have several "modes"*/
   st->frameSize = mode->frameSize;
   st->windowSize = mode->windowSize;
   st->nbSubframes=mode->frameSize/mode->subframeSize;
   st->subframeSize=mode->subframeSize;
   st->lpcSize = mode->lpcSize;
   st->bufSize = mode->bufSize;
   st->gamma1=mode->gamma1;
   st->gamma2=mode->gamma2;


   st->inBuf = malloc(st->bufSize*sizeof(float));
   st->frame = st->inBuf + st->bufSize - st->windowSize;
   st->excBuf = malloc(st->bufSize*sizeof(float));
   st->exc = st->excBuf + st->bufSize - st->windowSize;
   for (i=0;i<st->bufSize;i++)
      st->inBuf[i]=0;
   for (i=0;i<st->bufSize;i++)
      st->excBuf[i]=0;
}

void decoder_destroy(DecState *st)
{
   free(st->inBuf);
   free(st->excBuf);
}

void decode(DecState *st, FrameBits *bits, float *out)
{
}
