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

#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif


#define sqr(x) ((x)*(x))
#define min(a,b) ((a) < (b) ? (a) : (b))

void *nb_encoder_init(SpeexMode *m)
{
   EncState *st;
   SpeexNBMode *mode;
   int i;
   float tmp;

   mode=m->mode;
   st = malloc(sizeof(EncState));
   st->mode=m;
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

   st->exc2Buf = calloc(st->bufSize,sizeof(float));
   st->exc2 = st->exc2Buf + st->bufSize - st->windowSize;

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

   st->pi_gain = calloc(st->nbSubframes, sizeof(float));

   st->pitch = calloc(st->nbSubframes, sizeof(int));
   return st;
}

void nb_encoder_destroy(void *state)
{
   EncState *st=state;
   /* Free all allocated memory */
   free(st->inBuf);
   free(st->excBuf);
   free(st->swBuf);
   free(st->os_filt);
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
   
   free(st);
}

void nb_encode(void *state, float *in, SpeexBits *bits)
{
   EncState *st;
   int i, sub, roots;
   float error;
   
   st=state;
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
      tmp = (.5 + sub)/st->nbSubframes;
      for (i=0;i<st->lpcSize;i++)
         st->interp_lsp[i] = (1-tmp)*st->old_lsp[i] + tmp*st->lsp[i];
      for (i=0;i<st->lpcSize;i++)
         st->interp_qlsp[i] = (1-tmp)*st->old_qlsp[i] + tmp*st->qlsp[i];

      if (0) {
         float *h=PUSH(st->stack, 8);
         for (i=0;i<8;i++)
            h[i]=0;
         h[0]=1;
         
         residue_zero(h, st->bw_lpc1, h, 8, st->lpcSize);
         syn_filt_zero(h, st->interp_qlpc, h, 8, st->lpcSize);
         syn_filt_zero(h, st->bw_lpc2, h, 8, st->lpcSize);
         print_vec(h, 8, "lpc_resp");
         POP(st->stack);
      }
      
      /* Compute interpolated LPCs (quantized and unquantized) */
      for (i=0;i<st->lpcSize;i++)
         st->interp_lsp[i] = cos(st->interp_lsp[i]);
      lsp_to_lpc(st->interp_lsp, st->interp_lpc, st->lpcSize,st->stack);

      for (i=0;i<st->lpcSize;i++)
         st->interp_qlsp[i] = cos(st->interp_qlsp[i]);
      lsp_to_lpc(st->interp_qlsp, st->interp_qlpc, st->lpcSize, st->stack);

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
#ifdef DEBUG
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
#endif
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
         exc[i]=0;

      /* Long-term prediction */
      pitch = st->ltp_quant(target, sw, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2,
                    exc, st->ltp_params, st->min_pitch, st->max_pitch, 
                    st->lpcSize, st->subframeSize, bits, st->stack, exc2);
      st->pitch[sub]=pitch;

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
      if (0)
      {
      /* Perform innovation search */
      st->innovation_quant(target, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2,
                           st->innovation_params, st->lpcSize,
                           st->subframeSize, exc, bits, st->stack);
      }
      else
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

         {
            int qe = (int)(floor(7*log(ener)));
            if (qe<0)
               qe=0;
            if (qe>63)
               qe=63;
            ener = exp(qe/7.0);
            speex_bits_pack(bits, qe, 6);
         }
         ener_1 = 1/ener;
         
         for (i=0;i<st->subframeSize;i++)
            target[i]*=ener_1;
#if 1
         st->innovation_quant(target, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2, 
                                st->innovation_params, st->lpcSize, st->subframeSize, 
                                innov, bits, st->stack);
         
         for (i=0;i<st->subframeSize;i++)
            exc[i] += innov[i]*ener;
#else
         for (i=0;i<st->subframeSize;i++)
            exc[i] += st->buf2[i];
#endif
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

      printf ("seg SNR = %f\n", snr);

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

   st->post_filter_func = mode->post_filter_func;
   st->post_filter_params = mode->post_filter_params;
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

   st->pi_gain = calloc(st->nbSubframes, sizeof(float));
   
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
   free(st->pi_gain);
   
   free(state);
}

void nb_decode(void *state, SpeexBits *bits, float *out, int lost)
{
   DecState *st;
   int i, sub;
   int pitch;
   float pitch_gain[3];

   st=state;

   memmove(st->inBuf, st->inBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(float));
   memmove(st->excBuf, st->excBuf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(float));
   memmove(st->exc2Buf, st->exc2Buf+st->frameSize, (st->bufSize-st->frameSize)*sizeof(float));


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
      tmp = (.5 + sub)/st->nbSubframes;
      for (i=0;i<st->lpcSize;i++)
         st->interp_qlsp[i] = (1-tmp)*st->old_qlsp[i] + tmp*st->qlsp[i];

      /* Compute interpolated LPCs (unquantized) */
      for (i=0;i<st->lpcSize;i++)
         st->interp_qlsp[i] = cos(st->interp_qlsp[i]);
      lsp_to_lpc(st->interp_qlsp, st->interp_qlpc, st->lpcSize, st->stack);

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
      st->ltp_unquant(exc, st->min_pitch, st->max_pitch, st->ltp_params, st->subframeSize, &pitch, &pitch_gain[0], bits, st->stack, lost);
      

      {
         int q_energy;
         float ener;
         float *innov;
         
         innov = PUSH(st->stack, st->subframeSize);
         for (i=0;i<st->subframeSize;i++)
            innov[i]=0;

         q_energy = speex_bits_unpack_unsigned(bits, 6);
         ener = exp(q_energy/7.0);
         /*printf ("unquant_energy: %d %f\n", q_energy, ener);*/
         
         /*Fixed codebook contribution*/
         st->innovation_unquant(innov, st->innovation_params, st->subframeSize, bits, st->stack);

         for (i=0;i<st->subframeSize;i++)
            exc[i]+=ener*innov[i];

         POP(st->stack);
      }

      for (i=0;i<st->subframeSize;i++)
         exc2[i]=exc[i];

      if (st->pf_enabled)
         st->post_filter_func(exc, exc2, st->interp_qlpc, st->lpcSize, st->subframeSize,
                              pitch, pitch_gain, st->post_filter_params, st->stack);

      /*Compute decoded signal*/
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
