/* Copyright (C) 2002 Jean-Marc Valin 
   File: sb_celp.c

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


#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "nb_celp.h"
#include "sb_celp.h"
#include "stdlib.h"
#include "filters.h"
#include "lpc.h"
#include "lsp.h"
#include "stack_alloc.h"
#include "cb_search.h"
#include "quant_lsp.h"
#include "vq.h"
#include "ltp.h"
#include "misc.h"

#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif

#define sqr(x) ((x)*(x))

#define SUBMODE(x) st->submodes[st->submodeID]->x

#define QMF_ORDER 64
static float h0[64] = {
   3.596189e-05, -0.0001123515,
   -0.0001104587, 0.0002790277,
   0.0002298438, -0.0005953563,
   -0.0003823631, 0.00113826,
   0.0005308539, -0.001986177,
   -0.0006243724, 0.003235877,
   0.0005743159, -0.004989147,
   -0.0002584767, 0.007367171,
   -0.0004857935, -0.01050689,
   0.001894714, 0.01459396,
   -0.004313674, -0.01994365,
   0.00828756, 0.02716055,
   -0.01485397, -0.03764973,
   0.026447, 0.05543245,
   -0.05095487, -0.09779096,
   0.1382363, 0.4600981,
   0.4600981, 0.1382363,
   -0.09779096, -0.05095487,
   0.05543245, 0.026447,
   -0.03764973, -0.01485397,
   0.02716055, 0.00828756,
   -0.01994365, -0.004313674,
   0.01459396, 0.001894714,
   -0.01050689, -0.0004857935,
   0.007367171, -0.0002584767,
   -0.004989147, 0.0005743159,
   0.003235877, -0.0006243724,
   -0.001986177, 0.0005308539,
   0.00113826, -0.0003823631,
   -0.0005953563, 0.0002298438,
   0.0002790277, -0.0001104587,
   -0.0001123515, 3.596189e-05
};

static float h1[64] = {
   3.596189e-05, 0.0001123515,
   -0.0001104587, -0.0002790277,
   0.0002298438, 0.0005953563,
   -0.0003823631, -0.00113826,
   0.0005308539, 0.001986177,
   -0.0006243724, -0.003235877,
   0.0005743159, 0.004989147,
   -0.0002584767, -0.007367171,
   -0.0004857935, 0.01050689,
   0.001894714, -0.01459396,
   -0.004313674, 0.01994365,
   0.00828756, -0.02716055,
   -0.01485397, 0.03764973,
   0.026447, -0.05543245,
   -0.05095487, 0.09779096,
   0.1382363, -0.4600981,
   0.4600981, -0.1382363,
   -0.09779096, 0.05095487,
   0.05543245, -0.026447,
   -0.03764973, 0.01485397,
   0.02716055, -0.00828756,
   -0.01994365, 0.004313674,
   0.01459396, -0.001894714,
   -0.01050689, 0.0004857935,
   0.007367171, 0.0002584767,
   -0.004989147, -0.0005743159,
   0.003235877, 0.0006243724,
   -0.001986177, -0.0005308539,
   0.00113826, 0.0003823631,
   -0.0005953563, -0.0002298438,
   0.0002790277, 0.0001104587,
   -0.0001123515, -3.596189e-05
};

void *sb_encoder_init(SpeexMode *m)
{
   int i;
   SBEncState *st;
   SpeexSBMode *mode;

   st = (SBEncState*)speex_alloc(sizeof(SBEncState));
   st->mode = m;
   mode = (SpeexSBMode*)m->mode;

   st->st_low = speex_encoder_init(mode->nb_mode);
   st->full_frame_size = 2*mode->frameSize;
   st->frame_size = mode->frameSize;
   st->subframeSize = mode->subframeSize;
   st->nbSubframes = mode->frameSize/mode->subframeSize;
   st->windowSize = st->frame_size*3/2;
   st->lpcSize=mode->lpcSize;
   st->bufSize=mode->bufSize;

   st->submodes=mode->submodes;
   st->submodeID=mode->defaultSubmode;
   
   {
      /* FIXME: Should do this without explicit reference to a mode */
      int mod=6;
      speex_encoder_ctl(st->st_low, SPEEX_SET_MODE, &mod);
   }

   st->lag_factor = mode->lag_factor;
   st->lpc_floor = mode->lpc_floor;
   st->gamma1=mode->gamma1;
   st->gamma2=mode->gamma2;
   st->first=1;
   st->stack = speex_alloc(20000*sizeof(float));

   st->x0d=(float*)speex_alloc(st->frame_size*sizeof(float));
   st->x1d=(float*)speex_alloc(st->frame_size*sizeof(float));
   st->high=(float*)speex_alloc(st->full_frame_size*sizeof(float));
   st->y0=(float*)speex_alloc(st->full_frame_size*sizeof(float));
   st->y1=(float*)speex_alloc(st->full_frame_size*sizeof(float));

   st->h0_mem=(float*)speex_alloc(QMF_ORDER*sizeof(float));
   st->h1_mem=(float*)speex_alloc(QMF_ORDER*sizeof(float));
   st->g0_mem=(float*)speex_alloc(QMF_ORDER*sizeof(float));
   st->g1_mem=(float*)speex_alloc(QMF_ORDER*sizeof(float));

   st->buf=(float*)speex_alloc(st->windowSize*sizeof(float));
   st->excBuf=(float*)speex_alloc(st->bufSize*sizeof(float));
   st->exc = st->excBuf + st->bufSize - st->windowSize;

   st->res=(float*)speex_alloc(st->frame_size*sizeof(float));
   st->sw=(float*)speex_alloc(st->frame_size*sizeof(float));
   st->target=(float*)speex_alloc(st->frame_size*sizeof(float));
   /*Asymetric "pseudo-Hamming" window*/
   {
      int part1, part2;
      part1 = st->subframeSize*7/2;
      part2 = st->subframeSize*5/2;
      st->window = (float*)speex_alloc(st->windowSize*sizeof(float));
      for (i=0;i<part1;i++)
         st->window[i]=.54-.46*cos(M_PI*i/part1);
      for (i=0;i<part2;i++)
         st->window[part1+i]=.54+.46*cos(M_PI*i/part2);
   }

   st->lagWindow = (float*)speex_alloc((st->lpcSize+1)*sizeof(float));
   for (i=0;i<st->lpcSize+1;i++)
      st->lagWindow[i]=exp(-.5*sqr(2*M_PI*st->lag_factor*i));

   st->rc = (float*)speex_alloc(st->lpcSize*sizeof(float));
   st->autocorr = (float*)speex_alloc((st->lpcSize+1)*sizeof(float));
   st->lpc = (float*)speex_alloc((st->lpcSize+1)*sizeof(float));
   st->bw_lpc1 = (float*)speex_alloc((st->lpcSize+1)*sizeof(float));
   st->bw_lpc2 = (float*)speex_alloc((st->lpcSize+1)*sizeof(float));
   st->lsp = (float*)speex_alloc(st->lpcSize*sizeof(float));
   st->qlsp = (float*)speex_alloc(st->lpcSize*sizeof(float));
   st->old_lsp = (float*)speex_alloc(st->lpcSize*sizeof(float));
   st->old_qlsp = (float*)speex_alloc(st->lpcSize*sizeof(float));
   st->interp_lsp = (float*)speex_alloc(st->lpcSize*sizeof(float));
   st->interp_qlsp = (float*)speex_alloc(st->lpcSize*sizeof(float));
   st->interp_lpc = (float*)speex_alloc((st->lpcSize+1)*sizeof(float));
   st->interp_qlpc = (float*)speex_alloc((st->lpcSize+1)*sizeof(float));
   st->pi_gain = (float*)speex_alloc(st->nbSubframes*sizeof(float));

   st->mem_sp = (float*)speex_alloc(st->lpcSize*sizeof(float));
   st->mem_sp2 = (float*)speex_alloc(st->lpcSize*sizeof(float));
   st->mem_sw = (float*)speex_alloc(st->lpcSize*sizeof(float));
   st->complexity=2;

   return st;
}

void sb_encoder_destroy(void *state)
{
   SBEncState *st=(SBEncState*)state;

   speex_encoder_destroy(st->st_low);
   speex_free(st->x0d);
   speex_free(st->x1d);
   speex_free(st->high);
   speex_free(st->y0);
   speex_free(st->y1);
   speex_free(st->h0_mem);
   speex_free(st->h1_mem);
   speex_free(st->g0_mem);
   speex_free(st->g1_mem);
   
   speex_free(st->buf);
   speex_free(st->window);
   speex_free(st->excBuf);
   speex_free(st->sw);
   speex_free(st->res);
   speex_free(st->target);
   speex_free(st->lagWindow);
   speex_free(st->rc);
   speex_free(st->autocorr);
   speex_free(st->lpc);
   speex_free(st->bw_lpc1);
   speex_free(st->bw_lpc2);
   speex_free(st->lsp);
   speex_free(st->qlsp);
   speex_free(st->old_lsp);
   speex_free(st->old_qlsp);
   speex_free(st->interp_lsp);
   speex_free(st->interp_qlsp);
   speex_free(st->interp_lpc);
   speex_free(st->interp_qlpc);

   speex_free(st->mem_sp);
   speex_free(st->mem_sp2);
   speex_free(st->mem_sw);
   speex_free(st->pi_gain);

   speex_free(st->stack);

   speex_free(st);
}


void sb_encode(void *state, float *in, SpeexBits *bits)
{
   SBEncState *st;
   int i, roots, sub;
   void *stack;
   float *mem, *innov, *syn_resp;
   float *low_pi_gain, *low_exc, *low_innov;

   st = (SBEncState*)state;
   stack=st->stack;

   /* Compute the two sub-bands by filtering with h0 and h1*/
   qmf_decomp(in, h0, st->x0d, st->x1d, st->full_frame_size, QMF_ORDER, st->h0_mem, stack);
    
   /* Encode the narrowband part*/
   speex_encode(st->st_low, st->x0d, bits);

   speex_bits_pack(bits, 1, 1);
   speex_bits_pack(bits, st->submodeID, SB_SUBMODE_BITS);

   /* High-band buffering / sync with low band */
   for (i=0;i<st->windowSize-st->frame_size;i++)
      st->high[i] = st->high[st->frame_size+i];
   for (i=0;i<st->frame_size;i++)
      st->high[st->windowSize-st->frame_size+i]=st->x1d[i];

   speex_move(st->excBuf, st->excBuf+st->frame_size, (st->bufSize-st->frame_size)*sizeof(float));


   low_pi_gain = PUSH(stack, st->nbSubframes, float);
   low_exc = PUSH(stack, st->frame_size, float);
   low_innov = PUSH(stack, st->frame_size, float);
   speex_encoder_ctl(st->st_low, SPEEX_GET_PI_GAIN, low_pi_gain);
   speex_encoder_ctl(st->st_low, SPEEX_GET_EXC, low_exc);
   speex_encoder_ctl(st->st_low, SPEEX_GET_INNOV, low_innov);
   
   /* Start encoding the high-band */
   for (i=0;i<st->windowSize;i++)
      st->buf[i] = st->high[i] * st->window[i];

   /* Compute auto-correlation */
   _spx_autocorr(st->buf, st->autocorr, st->lpcSize+1, st->windowSize);

   st->autocorr[0] += 1;        /* prevents NANs */
   st->autocorr[0] *= st->lpc_floor; /* Noise floor in auto-correlation domain */
   /* Lag windowing: equivalent to filtering in the power-spectrum domain */
   for (i=0;i<st->lpcSize+1;i++)
      st->autocorr[i] *= st->lagWindow[i];

   /* Levinson-Durbin */
   wld(st->lpc+1, st->autocorr, st->rc, st->lpcSize);
   st->lpc[0]=1;

   /* LPC to LSPs (x-domain) transform */
   roots=lpc_to_lsp (st->lpc, st->lpcSize, st->lsp, 15, 0.2, stack);
   if (roots!=st->lpcSize)
   {
      roots = lpc_to_lsp (st->lpc, st->lpcSize, st->lsp, 11, 0.02, stack);
      if (roots!=st->lpcSize) {
         /*fprintf (stderr, "roots!=st->lpcSize (found only %d roots)\n", roots);*/

         /*If we can't find all LSP's, do some damage control and use a flat filter*/
         for (i=0;i<st->lpcSize;i++)
         {
            st->lsp[i]=cos(M_PI*((float)(i+1))/(st->lpcSize+1));
         }
      }
   }

   /* x-domain to angle domain*/
   for (i=0;i<st->lpcSize;i++)
      st->lsp[i] = acos(st->lsp[i]);

   /* If null mode (no transmission), just set a couple things to zero*/
   if (st->submodes[st->submodeID] == NULL)
   {
      for (i=0;i<st->frame_size;i++)
         st->exc[i]=st->sw[i]=0;

      for (i=0;i<st->lpcSize;i++)
         st->mem_sw[i]=0;
      st->first=1;

      /* Final signal synthesis from excitation */
      iir_mem2(st->exc, st->interp_qlpc, st->high, st->subframeSize, st->lpcSize, st->mem_sp);

#ifndef RELEASE

      /* Reconstruct the original */
      fir_mem_up(st->x0d, h0, st->y0, st->full_frame_size, QMF_ORDER, st->g0_mem, stack);
      fir_mem_up(st->high, h1, st->y1, st->full_frame_size, QMF_ORDER, st->g1_mem, stack);

      for (i=0;i<st->full_frame_size;i++)
         in[i]=2*(st->y0[i]-st->y1[i]);
#endif

      return;

   }


   /* LSP quantization */
   SUBMODE(lsp_quant)(st->lsp, st->qlsp, st->lpcSize, bits);
   
   /*printf ("high_lsp:");
   for (i=0;i<st->lpcSize;i++)
      printf (" %f", st->lsp[i]);
      printf ("\n");*/
   /*for (i=0;i<st->lpcSize;i++)
     st->qlsp[i]=st->lsp[i];*/
   

   if (st->first)
   {
      for (i=0;i<st->lpcSize;i++)
         st->old_lsp[i] = st->lsp[i];
      for (i=0;i<st->lpcSize;i++)
         st->old_qlsp[i] = st->qlsp[i];
   }
   
   mem=PUSH(stack, st->lpcSize, float);
   syn_resp=PUSH(stack, st->subframeSize, float);
   innov = PUSH(stack, st->subframeSize, float);

   for (sub=0;sub<st->nbSubframes;sub++)
   {
      float *exc, *sp, *res, *target, *sw, tmp, filter_ratio;
      int offset;
      float rl, rh, eh=0, el=0;
      int fold;

      offset = st->subframeSize*sub;
      sp=st->high+offset;
      exc=st->exc+offset;
      res=st->res+offset;
      target=st->target+offset;
      sw=st->sw+offset;
      
      /* LSP interpolation (quantized and unquantized) */
      tmp = (1.0 + sub)/st->nbSubframes;
      for (i=0;i<st->lpcSize;i++)
         st->interp_lsp[i] = (1-tmp)*st->old_lsp[i] + tmp*st->lsp[i];
      for (i=0;i<st->lpcSize;i++)
         st->interp_qlsp[i] = (1-tmp)*st->old_qlsp[i] + tmp*st->qlsp[i];
      
      /* Compute interpolated LPCs (quantized and unquantized) */
      for (i=0;i<st->lpcSize;i++)
         st->interp_lsp[i] = cos(st->interp_lsp[i]);
      for (i=0;i<st->lpcSize;i++)
         st->interp_qlsp[i] = cos(st->interp_qlsp[i]);

      lsp_enforce_margin(st->interp_lsp, st->lpcSize, .002);
      lsp_enforce_margin(st->interp_qlsp, st->lpcSize, .002);

      lsp_to_lpc(st->interp_lsp, st->interp_lpc, st->lpcSize,stack);
      lsp_to_lpc(st->interp_qlsp, st->interp_qlpc, st->lpcSize, stack);

      bw_lpc(st->gamma1, st->interp_lpc, st->bw_lpc1, st->lpcSize);
      bw_lpc(st->gamma2, st->interp_lpc, st->bw_lpc2, st->lpcSize);

      /* Compute mid-band (4000 Hz for wideband) response of low-band and high-band
         filters */
      rl=rh=0;
      tmp=1;
      st->pi_gain[sub]=0;
      for (i=0;i<=st->lpcSize;i++)
      {
         rh += tmp*st->interp_qlpc[i];
         tmp = -tmp;
         st->pi_gain[sub]+=st->interp_qlpc[i];
      }
      rl = low_pi_gain[sub];
      rl=1/(fabs(rl)+.01);
      rh=1/(fabs(rh)+.01);
      /* Compute ratio, will help predict the gain */
      filter_ratio=fabs(.01+rh)/(.01+fabs(rl));

      fold = filter_ratio<5;
      /*printf ("filter_ratio %f\n", filter_ratio);*/
      fold=0;

      /* Compute "real excitation" */
      fir_mem2(sp, st->interp_qlpc, exc, st->subframeSize, st->lpcSize, st->mem_sp2);
      /* Compute energy of low-band and high-band excitation */
      for (i=0;i<st->subframeSize;i++)
         eh+=sqr(exc[i]);

      if (!SUBMODE(innovation_quant)) {/* 1 for spectral folding excitation, 0 for stochastic */
         float g;
         /*speex_bits_pack(bits, 1, 1);*/
         for (i=0;i<st->subframeSize;i++)
            el+=sqr(low_innov[offset+i]);

         /* Gain to use if we want to use the low-band excitation for high-band */
         g=eh/(.01+el);
         g=sqrt(g);

         g *= filter_ratio;

         /* Gain quantization */
         {
            int quant = (int) floor(.5 + 27 + 8.0 * log((g+.0001)));
            if (quant<0)
               quant=0;
            if (quant>31)
               quant=31;
            speex_bits_pack(bits, quant, 5);
            g= .1*exp(quant/9.4);
         }
         /*printf ("folding gain: %f\n", g);*/
         g /= filter_ratio;

      } else {
         float gc, scale, scale_1;

         for (i=0;i<st->subframeSize;i++)
            el+=sqr(low_exc[offset+i]);
         /*speex_bits_pack(bits, 0, 1);*/

         gc = sqrt(1+eh)*filter_ratio/sqrt((1+el)*st->subframeSize);
         {
            int qgc = (int)floor(.5+3.7*(log(gc)+2));
            if (qgc<0)
               qgc=0;
            if (qgc>15)
               qgc=15;
            speex_bits_pack(bits, qgc, 4);
            gc = exp((1/3.7)*qgc-2);
         }

         scale = gc*sqrt(1+el)/filter_ratio;
         scale_1 = 1/scale;
         if (0 && rand()%5==0)
         {
            float sc = 1/sqrt(.1+eh/st->subframeSize);
            if (rand()&1)
               sc=-sc;
            for (i=0;i<st->subframeSize;i++)
            {
               float tmp=exc[i]*sc;
               if (i%8==0)
                  printf ("\nhexc");
               printf (" %f", tmp);
            }
         }

         for (i=0;i<st->subframeSize;i++)
            exc[i]=0;
         exc[0]=1;
         syn_percep_zero(exc, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2, syn_resp, st->subframeSize, st->lpcSize, stack);
         
         /* Reset excitation */
         for (i=0;i<st->subframeSize;i++)
            exc[i]=0;
         
         /* Compute zero response (ringing) of A(z/g1) / ( A(z/g2) * Aq(z) ) */
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
           exc[i]=0;


         for (i=0;i<st->subframeSize;i++)
            target[i]*=scale_1;
         
         /* Reset excitation */
         for (i=0;i<st->subframeSize;i++)
            innov[i]=0;

         /*print_vec(target, st->subframeSize, "\ntarget");*/
         SUBMODE(innovation_quant)(target, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2, 
                                   SUBMODE(innovation_params), st->lpcSize, st->subframeSize, 
                                   innov, syn_resp, bits, stack, (st->complexity+1)>>1);
         /*print_vec(target, st->subframeSize, "after");*/

         for (i=0;i<st->subframeSize;i++)
            exc[i] += innov[i]*scale;

         if (SUBMODE(double_codebook)) {
            void *tmp_stack=stack;
            float *innov2 = PUSH(tmp_stack, st->subframeSize, float);
            for (i=0;i<st->subframeSize;i++)
               innov2[i]=0;
            for (i=0;i<st->subframeSize;i++)
               target[i]*=2.5;
            SUBMODE(innovation_quant)(target, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2, 
                                      SUBMODE(innovation_params), st->lpcSize, st->subframeSize, 
                                      innov2, syn_resp, bits, tmp_stack, (st->complexity+1)>>1);
            for (i=0;i<st->subframeSize;i++)
               innov2[i]*=scale*(1/2.5);
            for (i=0;i<st->subframeSize;i++)
               exc[i] += innov2[i];
         }


         if (0) {
            float en=0;
            for (i=0;i<st->subframeSize;i++)
               en+=exc[i]*exc[i];
            en=sqrt(eh/(1+en));
            for (i=0;i<st->subframeSize;i++)
               exc[i]*=en;
            printf ("correction high: %f\n", en);
         }

      }

         /*Keep the previous memory*/
         for (i=0;i<st->lpcSize;i++)
            mem[i]=st->mem_sp[i];
         /* Final signal synthesis from excitation */
         iir_mem2(exc, st->interp_qlpc, sp, st->subframeSize, st->lpcSize, st->mem_sp);
         
         /* Compute weighted signal again, from synthesized speech (not sure it's the right thing) */
         filter_mem2(sp, st->bw_lpc1, st->bw_lpc2, sw, st->subframeSize, st->lpcSize, st->mem_sw);
   }


#ifndef RELEASE

   /* Reconstruct the original */
   fir_mem_up(st->x0d, h0, st->y0, st->full_frame_size, QMF_ORDER, st->g0_mem, stack);
   fir_mem_up(st->high, h1, st->y1, st->full_frame_size, QMF_ORDER, st->g1_mem, stack);

   for (i=0;i<st->full_frame_size;i++)
      in[i]=2*(st->y0[i]-st->y1[i]);
#endif
   for (i=0;i<st->lpcSize;i++)
      st->old_lsp[i] = st->lsp[i];
   for (i=0;i<st->lpcSize;i++)
      st->old_qlsp[i] = st->qlsp[i];

   st->first=0;
}





void *sb_decoder_init(SpeexMode *m)
{
   SBDecState *st;
   SpeexSBMode *mode;
   st = (SBDecState*)speex_alloc(sizeof(SBDecState));
   st->mode = m;
   mode=(SpeexSBMode*)m->mode;

   st->st_low = speex_decoder_init(mode->nb_mode);
   st->full_frame_size = 2*mode->frameSize;
   st->frame_size = mode->frameSize;
   st->subframeSize = mode->subframeSize;
   st->nbSubframes = mode->frameSize/mode->subframeSize;
   st->lpcSize=8;

   st->submodes=mode->submodes;
   st->submodeID=mode->defaultSubmode;

   st->first=1;
   st->stack = speex_alloc(20000*sizeof(float));

   st->x0d=(float*)speex_alloc(st->frame_size*sizeof(float));
   st->x1d=(float*)speex_alloc(st->frame_size*sizeof(float));
   st->high=(float*)speex_alloc(st->full_frame_size*sizeof(float));
   st->y0=(float*)speex_alloc(st->full_frame_size*sizeof(float));
   st->y1=(float*)speex_alloc(st->full_frame_size*sizeof(float));

   st->h0_mem=(float*)speex_alloc(QMF_ORDER*sizeof(float));
   st->h1_mem=(float*)speex_alloc(QMF_ORDER*sizeof(float));
   st->g0_mem=(float*)speex_alloc(QMF_ORDER*sizeof(float));
   st->g1_mem=(float*)speex_alloc(QMF_ORDER*sizeof(float));

   st->exc=(float*)speex_alloc(st->frame_size*sizeof(float));

   st->qlsp = (float*)speex_alloc(st->lpcSize*sizeof(float));
   st->old_qlsp = (float*)speex_alloc(st->lpcSize*sizeof(float));
   st->interp_qlsp = (float*)speex_alloc(st->lpcSize*sizeof(float));
   st->interp_qlpc = (float*)speex_alloc((st->lpcSize+1)*sizeof(float));

   st->pi_gain = (float*)speex_alloc(st->nbSubframes*sizeof(float));
   st->mem_sp = (float*)speex_alloc(st->lpcSize*sizeof(float));
   return st;
}

void sb_decoder_destroy(void *state)
{
   SBDecState *st;
   st = (SBDecState*)state;
   speex_decoder_destroy(st->st_low);
   speex_free(st->x0d);
   speex_free(st->x1d);
   speex_free(st->high);
   speex_free(st->y0);
   speex_free(st->y1);
   speex_free(st->h0_mem);
   speex_free(st->h1_mem);
   speex_free(st->g0_mem);
   speex_free(st->g1_mem);
   
   speex_free(st->exc);
   speex_free(st->qlsp);
   speex_free(st->old_qlsp);
   speex_free(st->interp_qlsp);
   speex_free(st->interp_qlpc);
   speex_free(st->pi_gain);

   speex_free(st->mem_sp);

   speex_free(st->stack);

   speex_free(state);
}

static void sb_decode_lost(SBDecState *st, float *out, void *stack)
{
   int i;
   for (i=0;i<st->frame_size;i++)
      st->exc[i]*=0.8;
   
   st->first=1;
   
   /* Final signal synthesis from excitation */
   iir_mem2(st->exc, st->interp_qlpc, st->high, st->subframeSize, st->lpcSize, st->mem_sp);
   
   /* Reconstruct the original */
   fir_mem_up(st->x0d, h0, st->y0, st->full_frame_size, QMF_ORDER, st->g0_mem, stack);
   fir_mem_up(st->high, h1, st->y1, st->full_frame_size, QMF_ORDER, st->g1_mem, stack);

   for (i=0;i<st->full_frame_size;i++)
      out[i]=2*(st->y0[i]-st->y1[i]);
   
   return;
}

int sb_decode(void *state, SpeexBits *bits, float *out)
{
   int i, sub;
   SBDecState *st;
   int wideband;
   int ret;
   void *stack;
   float *low_pi_gain, *low_exc, *low_innov;

   st = (SBDecState*)state;
   stack=st->stack;

   /* Decode the low-band */
   ret = speex_decode(st->st_low, bits, st->x0d);

   /* If error decoding the narrowband part, propagate error */
   if (ret!=0)
   {
      return ret;
   }

   if (!bits)
   {
      sb_decode_lost(st, out, stack);
      return 0;
   }

   /*Check "wideband bit"*/
   wideband = speex_bits_peek(bits);
   if (wideband)
   {
      /*Regular wideband frame, read the submode*/
      wideband = speex_bits_unpack_unsigned(bits, 1);
      st->submodeID = speex_bits_unpack_unsigned(bits, SB_SUBMODE_BITS);
   } else
   {
      /*Was a narrowband frame, set "null submode"*/
      st->submodeID = 0;
   }

   for (i=0;i<st->frame_size;i++)
      st->exc[i]=0;

   /* If null mode (no transmission), just set a couple things to zero*/
   if (st->submodes[st->submodeID] == NULL)
   {
      for (i=0;i<st->frame_size;i++)
         st->exc[i]=0;

      st->first=1;

      /* Final signal synthesis from excitation */
      iir_mem2(st->exc, st->interp_qlpc, st->high, st->subframeSize, st->lpcSize, st->mem_sp);

      fir_mem_up(st->x0d, h0, st->y0, st->full_frame_size, QMF_ORDER, st->g0_mem, stack);
      fir_mem_up(st->high, h1, st->y1, st->full_frame_size, QMF_ORDER, st->g1_mem, stack);

      for (i=0;i<st->full_frame_size;i++)
         out[i]=2*(st->y0[i]-st->y1[i]);

      return 0;

   }

   low_pi_gain = PUSH(stack, st->nbSubframes, float);
   low_exc = PUSH(stack, st->frame_size, float);
   low_innov = PUSH(stack, st->frame_size, float);
   speex_decoder_ctl(st->st_low, SPEEX_GET_PI_GAIN, low_pi_gain);
   speex_decoder_ctl(st->st_low, SPEEX_GET_EXC, low_exc);
   speex_decoder_ctl(st->st_low, SPEEX_GET_INNOV, low_innov);

   SUBMODE(lsp_unquant)(st->qlsp, st->lpcSize, bits);
   
   if (st->first)
   {
      for (i=0;i<st->lpcSize;i++)
         st->old_qlsp[i] = st->qlsp[i];
   }
   
   for (sub=0;sub<st->nbSubframes;sub++)
   {
      float *exc, *sp, tmp, filter_ratio, el=0;
      int offset;
      float rl=0,rh=0;
      
      offset = st->subframeSize*sub;
      sp=st->high+offset;
      exc=st->exc+offset;
      
      /* LSP interpolation */
      tmp = (1.0 + sub)/st->nbSubframes;
      for (i=0;i<st->lpcSize;i++)
         st->interp_qlsp[i] = (1-tmp)*st->old_qlsp[i] + tmp*st->qlsp[i];

      /* LSPs to x-domain */
      for (i=0;i<st->lpcSize;i++)
         st->interp_qlsp[i] = cos(st->interp_qlsp[i]);

      lsp_enforce_margin(st->interp_qlsp, st->lpcSize, .002);

      /* LSP to LPC */
      lsp_to_lpc(st->interp_qlsp, st->interp_qlpc, st->lpcSize, stack);

      /* Calculate reponse ratio between the low and high filter in the middle
         of the band (4000 Hz) */
      
         tmp=1;
         st->pi_gain[sub]=0;
         for (i=0;i<=st->lpcSize;i++)
         {
            rh += tmp*st->interp_qlpc[i];
            tmp = -tmp;
            st->pi_gain[sub]+=st->interp_qlpc[i];
         }
         rl = low_pi_gain[sub];
         rl=1/(fabs(rl)+.01);
         rh=1/(fabs(rh)+.01);
         filter_ratio=fabs(.01+rh)/(.01+fabs(rl));
      
      
      for (i=0;i<st->subframeSize;i++)
         exc[i]=0;
      if (!SUBMODE(innovation_unquant))
      {
         float g;
         int quant;

         for (i=0;i<st->subframeSize;i++)
            el+=sqr(low_innov[offset+i]);
         quant = speex_bits_unpack_unsigned(bits, 5);
         g= exp(((float)quant-27)/8.0);
         
         /*printf ("unquant folding gain: %f\n", g);*/
         g /= filter_ratio;
         
         /* High-band excitation using the low-band excitation and a gain */
         for (i=0;i<st->subframeSize;i++)
            exc[i]=.8*g*low_innov[offset+i];
      } else {
         float gc, scale;
         int qgc = speex_bits_unpack_unsigned(bits, 4);
         for (i=0;i<st->subframeSize;i++)
            el+=sqr(low_exc[offset+i]);


         gc = exp((1/3.7)*qgc-2);

         scale = gc*sqrt(1+el)/filter_ratio;


         SUBMODE(innovation_unquant)(exc, SUBMODE(innovation_params), st->subframeSize, 
                                bits, stack);
         for (i=0;i<st->subframeSize;i++)
            exc[i]*=scale;

         if (SUBMODE(double_codebook)) {
            void *tmp_stack=stack;
            float *innov2 = PUSH(tmp_stack, st->subframeSize, float);
            for (i=0;i<st->subframeSize;i++)
               innov2[i]=0;
            SUBMODE(innovation_unquant)(innov2, SUBMODE(innovation_params), st->subframeSize, 
                                bits, tmp_stack);
            for (i=0;i<st->subframeSize;i++)
               innov2[i]*=scale*(1/2.5);
            for (i=0;i<st->subframeSize;i++)
               exc[i] += innov2[i];
         }

      }
      iir_mem2(exc, st->interp_qlpc, sp, st->subframeSize, st->lpcSize, st->mem_sp);

   }

   fir_mem_up(st->x0d, h0, st->y0, st->full_frame_size, QMF_ORDER, st->g0_mem, stack);
   fir_mem_up(st->high, h1, st->y1, st->full_frame_size, QMF_ORDER, st->g1_mem, stack);

   for (i=0;i<st->full_frame_size;i++)
      out[i]=2*(st->y0[i]-st->y1[i]);

   for (i=0;i<st->lpcSize;i++)
      st->old_qlsp[i] = st->qlsp[i];

   st->first=0;

   return 0;
}


void sb_encoder_ctl(void *state, int request, void *ptr)
{
   SBEncState *st;
   st=(SBEncState*)state;
   switch(request)
   {
   case SPEEX_GET_FRAME_SIZE:
      (*(int*)ptr) = st->full_frame_size;
      break;
   case SPEEX_SET_HIGH_MODE:
      st->submodeID = (*(int*)ptr);
      break;
   case SPEEX_SET_LOW_MODE:
      speex_encoder_ctl(st->st_low, SPEEX_SET_MODE, ptr);
      break;
   case SPEEX_SET_VBR:
      speex_encoder_ctl(st->st_low, SPEEX_SET_VBR, ptr);
      break;
   case SPEEX_SET_VBR_QUALITY:
      {
         int q;
         float qual = (*(float*)ptr)+1;
         if (qual>10)
            qual=10;
         q=(int)floor(.5+*(float*)ptr);
         if (q>10)
            q=10;
         speex_encoder_ctl(st->st_low, SPEEX_SET_VBR_QUALITY, &qual);
         speex_encoder_ctl(state, SPEEX_SET_QUALITY, &q);
         break;
      }
   case SPEEX_SET_QUALITY:
      {
         int nb_mode;
         int quality = (*(int*)ptr);
         switch (quality)
         {
         case 0:
            nb_mode=0;
            st->submodeID = 0;
            break;
         case 1:
            nb_mode=1;
            st->submodeID = 1;
            break;
         case 2:
            nb_mode=2;
            st->submodeID = 1;
            break;
         case 3:
            nb_mode=3;
            st->submodeID = 1;
            break;
         case 4:
            nb_mode=4;
            st->submodeID = 1;
            break;
         case 5:
            nb_mode=5;
            st->submodeID = 1;
            break;
         case 6:
            nb_mode=5;
            st->submodeID = 2;
            break;
         case 7:
            nb_mode=6;
            st->submodeID = 2;
            break;
         case 8:
            nb_mode=6;
            st->submodeID = 3;
            break;
         case 9:
            nb_mode=7;
            st->submodeID = 3;
            break;
         case 10:
            nb_mode=7;
            st->submodeID = 4;
            break;
         default:
            fprintf(stderr, "Unknown sb_ctl quality: %d\n", quality);
         }
         speex_encoder_ctl(st->st_low, SPEEX_SET_MODE, &nb_mode);
      }
      break;
   case SPEEX_SET_COMPLEXITY:
      speex_encoder_ctl(st->st_low, SPEEX_SET_COMPLEXITY, ptr);
      st->complexity = (*(int*)ptr);
      if (st->complexity<1)
         st->complexity=1;
      break;
   case SPEEX_GET_COMPLEXITY:
      (*(int*)ptr) = st->complexity;
      break;
   case SPEEX_GET_BITRATE:
      speex_encoder_ctl(st->st_low, request, ptr);
      if (st->submodes[st->submodeID])
         (*(int*)ptr) += 50*SUBMODE(bits_per_frame);
      else
         (*(int*)ptr) += 50*(SB_SUBMODE_BITS+1);
      break;
   case SPEEX_GET_PI_GAIN:
      {
         int i;
         float *g = (float*)ptr;
         for (i=0;i<st->nbSubframes;i++)
            g[i]=st->pi_gain[i];
      }
      break;
   case SPEEX_GET_EXC:
      {
         int i;
         float *e = (float*)ptr;
         for (i=0;i<st->full_frame_size;i++)
            e[i]=0;
         for (i=0;i<st->frame_size;i++)
            e[2*i]=2*st->exc[i];
      }
      break;
   case SPEEX_GET_INNOV:
      {
         int i;
         float *e = (float*)ptr;
         for (i=0;i<st->full_frame_size;i++)
            e[i]=0;
         for (i=0;i<st->frame_size;i++)
            e[2*i]=2*st->exc[i];
      }
      break;
   default:
      fprintf(stderr, "Unknown nb_ctl request: %d\n", request);
   }

}

void sb_decoder_ctl(void *state, int request, void *ptr)
{
   SBDecState *st;
   st=(SBDecState*)state;
   switch(request)
   {
   case SPEEX_GET_FRAME_SIZE:
      (*(int*)ptr) = st->full_frame_size;
      break;
   case SPEEX_SET_ENH:
      speex_decoder_ctl(st->st_low, request, ptr);
      break;
   case SPEEX_GET_BITRATE:
      speex_decoder_ctl(st->st_low, request, ptr);
      if (st->submodes[st->submodeID])
         (*(int*)ptr) += 50*SUBMODE(bits_per_frame);
      else
         (*(int*)ptr) += 50*(SB_SUBMODE_BITS+1);
      break;
   case SPEEX_SET_HANDLER:
      speex_decoder_ctl(st->st_low, SPEEX_SET_HANDLER, ptr);
      break;
   case SPEEX_SET_USER_HANDLER:
      speex_decoder_ctl(st->st_low, SPEEX_SET_USER_HANDLER, ptr);
      break;
   case SPEEX_GET_PI_GAIN:
      {
         int i;
         float *g = (float*)ptr;
         for (i=0;i<st->nbSubframes;i++)
            g[i]=st->pi_gain[i];
      }
      break;
   case SPEEX_GET_EXC:
      {
         int i;
         float *e = (float*)ptr;
         for (i=0;i<st->full_frame_size;i++)
            e[i]=0;
         for (i=0;i<st->frame_size;i++)
            e[2*i]=2*st->exc[i];
      }
      break;
   case SPEEX_GET_INNOV:
      {
         int i;
         float *e = (float*)ptr;
         for (i=0;i<st->full_frame_size;i++)
            e[i]=0;
         for (i=0;i<st->frame_size;i++)
            e[2*i]=2*st->exc[i];
      }
      break;
   default:
      fprintf(stderr, "Unknown sb_ctl request: %d\n", request);
   }

}
