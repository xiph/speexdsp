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
   st->lag_factor=mode->lag_factor;
   st->lpc_floor = mode->lpc_floor;
   st->preemph = mode->preemph;
  

   st->lsp_quant = mode->lsp_quant;
   st->ltp_quant = mode->ltp_quant;
   st->ltp_params = mode->ltp_params;
   st->innovation_quant = mode->innovation_quant;
   st->innovation_params = mode->innovation_params;

   st->pre_mem=0;
   st->pre_mem2=0;

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
   int i, sub, roots;
   float error;

   /* Copy new data in input buffer */
   memmove(st->inBuf, st->inBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(float));
   st->inBuf[st->bufSize-st->frameSize] = in[0] - st->preemph*st->pre_mem;
   for (i=1;i<st->frameSize;i++)
      st->inBuf[st->bufSize-st->frameSize+i] = in[i] - st->preemph*in[i-1];
   st->pre_mem = in[st->frameSize-1];

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
   
   /* LSP Quantization */
   st->lsp_quant(st->lsp, st->qlsp, st->lpcSize, bits);
   /*for (i=0;i<st->lpcSize;i++)
     st->qlsp[i]=st->lsp[i];*/
   /*printf ("LSP ");
   for (i=0;i<st->lpcSize;i++)
      printf ("%f ", st->lsp[i]);
   printf ("\n");
   printf ("QLSP ");
   for (i=0;i<st->lpcSize;i++)
      printf ("%f ", st->qlsp[i]);
   printf ("\n");*/
   /* Special case for first frame */
   if (st->first)
   {
      for (i=0;i<st->lpcSize;i++)
         st->old_lsp[i] = st->lsp[i];
      for (i=0;i<st->lpcSize;i++)
         st->old_qlsp[i] = st->qlsp[i];
   }

   /* Loop on sub-frames */
   for (sub=0;sub<st->nbSubframes;sub++)
   {
      float esig, enoise, snr, tmp;
      int   offset;
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
      if (st->gamma2>=0)
         bw_lpc(st->gamma2, st->interp_lpc, st->bw_lpc2, st->lpcSize);
      else
      {
         st->bw_lpc2[0]=1;
         st->bw_lpc2[1]=-st->preemph;
         for (i=2;i<=st->lpcSize;i++)
            st->bw_lpc2[i]=0;
      }
      printf ("\nlpc0 ");
      for (i=0;i<=st->lpcSize;i++)
         printf ("%f ", st->interp_lpc[i]);
      printf ("\nlpc1 ");
      for (i=0;i<=st->lpcSize;i++)
         printf ("%f ", st->bw_lpc1[i]);
      printf ("\nlpc2 ");
      for (i=0;i<=st->lpcSize;i++)
         printf ("%f ", st->bw_lpc2[i]);
      printf ("\n\n");

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
      
      esig=0;
      for (i=0;i<st->subframeSize;i++)
         esig+=sw[i]*sw[i];
      
      /* Compute target signal */
      for (i=0;i<st->subframeSize;i++)
         target[i]=sw[i]-res[i];

      for (i=0;i<st->subframeSize;i++)
         exc[i]=0;

      /* Long-term prediction */
#if 1
      st->ltp_quant(target, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2,
                    exc, st->ltp_params, st->min_pitch, st->max_pitch, 
                    st->lpcSize, st->subframeSize, bits, st->stack);
#else
      {
         float gain[3];
         int pitch;
         closed_loop_fractional_pitch(target, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2,                                      exc, st->os_filt, st->os_filt_ord2, st->os_fact, 35, 290,                                      &gain[0], &pitch, st->lpcSize,                                      st->subframeSize, st->stack);
      }
#endif
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
      printf ("pitch SNR = %f\n", snr);

#if 0 /*If set to 1, compute "real innovation" i.e. cheat to get perfect reconstruction*/
      syn_filt_zero(target, st->bw_lpc1, res, st->subframeSize, st->lpcSize);
      residue_zero(res, st->interp_qlpc, st->buf2, st->subframeSize, st->lpcSize);
      residue_zero(st->buf2, st->bw_lpc2, st->buf2, st->subframeSize, st->lpcSize);
      if (1||(snr>9 && (rand()%10==0)))
      {
         printf ("exc ");
         for (i=0;i<st->subframeSize;i++)
         {
            if (0&&i && i%8==0)
               printf ("\nexc ");
            printf ("%f ", st->buf2[i]);
         }
         printf ("\n");
      }
      for (i=0;i<st->subframeSize;i++)
         exc[i]+=st->buf2[i];
#else
      /* Perform a split-codebook search */
      st->innovation_quant(target, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2,
                           st->innovation_params, st->lpcSize,
                           st->subframeSize, exc, bits, st->stack);

#endif
      /* Compute weighted noise energy and SNR */
      enoise=0;
      for (i=0;i<st->subframeSize;i++)
         enoise += target[i]*target[i];
      snr = 10*log10((esig+1)/(enoise+1));
      printf ("seg SNR = %f\n", snr);

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

   /* The next frame will not be the first (Duh!) */
   st->first = 0;

   /* Replace input by synthesized speech */
   in[0] = st->frame[0] + st->preemph*st->pre_mem2;
   for (i=1;i<st->frameSize;i++)
     in[i]=st->frame[i] + st->preemph*in[i-1];
   st->pre_mem2=in[st->frameSize-1];
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
   st->min_pitch=mode->pitchStart;
   st->max_pitch=mode->pitchEnd;
   st->preemph = mode->preemph;

   st->pre_mem=0;
   st->lsp_unquant = mode->lsp_unquant;
   st->ltp_unquant = mode->ltp_unquant;
   st->ltp_params = mode->ltp_params;

   st->innovation_unquant = mode->innovation_unquant;
   st->innovation_params = mode->innovation_params;

   st->stack = calloc(10000, sizeof(float));

   st->inBuf = malloc(st->bufSize*sizeof(float));
   st->frame = st->inBuf + st->bufSize - st->windowSize;
   st->excBuf = malloc(st->bufSize*sizeof(float));
   st->exc = st->excBuf + st->bufSize - st->windowSize;
   for (i=0;i<st->bufSize;i++)
      st->inBuf[i]=0;
   for (i=0;i<st->bufSize;i++)
      st->excBuf[i]=0;

   st->interp_qlpc = malloc((st->lpcSize+1)*sizeof(float));
   st->qlsp = malloc(st->lpcSize*sizeof(float));
   st->old_qlsp = malloc(st->lpcSize*sizeof(float));
   st->interp_qlsp = malloc(st->lpcSize*sizeof(float));
   st->mem_sp = calloc(st->lpcSize, sizeof(float));

}

void decoder_destroy(DecState *st)
{
   free(st->inBuf);
   free(st->excBuf);
   free(st->interp_qlpc);
   free(st->qlsp);
   free(st->old_qlsp);
   free(st->interp_qlsp);
   free(st->stack);
   free(st->mem_sp);
}

void decode(DecState *st, FrameBits *bits, float *out)
{
   int i, sub;

   memmove(st->inBuf, st->inBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(float));
   memmove(st->excBuf, st->excBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(float));


   st->lsp_unquant(st->qlsp, st->lpcSize, bits);
   if (st->first)
   {
      for (i=0;i<st->lpcSize;i++)
         st->old_qlsp[i] = st->qlsp[i];
   }

   /*Loop on subframes */
   for (sub=0;sub<st->nbSubframes;sub++)
   {
      int offset;
      float *sp, *exc, tmp;
      
      /* Offset relative to start of frame */
      offset = st->subframeSize*sub;
      /* Original signal */
      sp=st->frame+offset;
      /* Excitation */
      exc=st->exc+offset;

      /* LSP interpolation (quantized and unquantized) */
      tmp = (.5 + sub)/st->nbSubframes;
      for (i=0;i<st->lpcSize;i++)
         st->interp_qlsp[i] = (1-tmp)*st->old_qlsp[i] + tmp*st->qlsp[i];

      /* Compute interpolated LPCs (unquantized) */
      for (i=0;i<st->lpcSize;i++)
         st->interp_qlsp[i] = cos(st->interp_qlsp[i]);
      lsp_to_lpc(st->interp_qlsp, st->interp_qlpc, st->lpcSize, st->stack);

      /* Reset excitation */
      for (i=0;i<st->subframeSize;i++)
         exc[i]=0;

      /*Adaptive codebook contribution*/
      st->ltp_unquant(exc, st->min_pitch, st->max_pitch, st->ltp_params, st->subframeSize, bits, st->stack);
       
      /*Fixed codebook contribution*/
      st->innovation_unquant(exc, st->innovation_params, st->subframeSize, bits, st->stack);

      /*Compute decoded signal*/
      syn_filt_mem(exc, st->interp_qlpc, sp, st->subframeSize, st->lpcSize, st->mem_sp);

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

}
