/* Copyright (C) 2002 Jean-Marc Valin 
   File: sb_celp.c

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


#include "speex.h"
#include "sb_celp.h"
#include "stdlib.h"
#include "filters.h"
#include <math.h>
#include "lpc.h"
#include "lsp.h"
#include <stdio.h>
#include "stack_alloc.h"
#include "cb_search.h"
#include "quant_lsp.h"

#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif

extern float stoc[];

#define sqr(x) ((x)*(x))

float quant_high_gain[16]={
   -2.387860,
   -1.504710,
   -0.988013,
   -0.610249,
   -0.310298,
   -0.050495,
   0.188963,
   0.413744,
   0.628971,
   0.840555,
   1.055630,
   1.283410,
   1.544990,
   1.855790,
   2.281910,
   3.002660
};


float quant_high_gain2[8] = {
   -1.51541,
   -0.70324,
   -0.17024,
   0.26748,
   0.67232,
   1.08402,
   1.56110,
   2.25160,
};

#if 0
#define QMF_ORDER 32
static float h0[32] = {
   0.0006910579, -0.001403793,
   -0.001268303, 0.004234195,
   0.001414246, -0.009458318,
   -0.0001303859, 0.01798145,
   -0.004187483, -0.03123862,
   0.01456844, 0.05294745,
   -0.03934878, -0.09980243,
   0.1285579, 0.4664053,
   0.4664053, 0.1285579,
   -0.09980243, -0.03934878,
   0.05294745, 0.01456844,
   -0.03123862, -0.004187483,
   0.01798145, -0.0001303859,
   -0.009458318, 0.001414246,
   0.004234195, -0.001268303,
   -0.001403793, 0.0006910579
};

static float h1[32] = {
   0.0006910579, 0.001403793,
   -0.001268303, -0.004234195,
   0.001414246, 0.009458318,
   -0.0001303859, -0.01798145,
   -0.004187483, 0.03123862,
   0.01456844, -0.05294745,
   -0.03934878, 0.09980243,
   0.1285579, -0.4664053,
   0.4664053, -0.1285579,
   -0.09980243, 0.03934878,
   0.05294745, -0.01456844,
   -0.03123862, 0.004187483,
   0.01798145, 0.0001303859,
   -0.009458318, -0.001414246,
   0.004234195, 0.001268303,
   -0.001403793, -0.0006910579
};
#else 
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
#endif

void sb_encoder_init(SBEncState *st, SpeexMode *mode)
{
   int i;
   encoder_init(&st->st_low, mode);
   st->full_frame_size = 2*st->st_low.frameSize;
   st->frame_size = st->st_low.frameSize;
   st->subframeSize = 40;
   st->nbSubframes = 4;
   st->windowSize = mode->windowSize;
   st->lpcSize=8;

   st->lag_factor = .002;
   st->lpc_floor = 1.0001;
   st->gamma1=.9;
   st->gamma2=.6;
   st->first=1;
   st->stack = calloc(10000, sizeof(float));

   st->x0=calloc(st->full_frame_size, sizeof(float));
   st->x1=calloc(st->full_frame_size, sizeof(float));
   st->x0d=calloc(st->frame_size, sizeof(float));
   st->x1d=calloc(st->frame_size, sizeof(float));
   st->high=calloc(st->full_frame_size, sizeof(float));
   st->y0=calloc(st->full_frame_size, sizeof(float));
   st->y1=calloc(st->full_frame_size, sizeof(float));

   st->h0_mem=calloc(QMF_ORDER, sizeof(float));
   st->h1_mem=calloc(QMF_ORDER, sizeof(float));
   st->g0_mem=calloc(QMF_ORDER, sizeof(float));
   st->g1_mem=calloc(QMF_ORDER, sizeof(float));

   st->buf=calloc(st->windowSize, sizeof(float));
   st->excBuf=calloc(2*st->frame_size, sizeof(float));
   st->exc=st->excBuf+st->frame_size;

   st->res=calloc(st->frame_size, sizeof(float));
   st->sw=calloc(st->frame_size, sizeof(float));
   st->target=calloc(st->frame_size, sizeof(float));
   st->window=calloc(st->windowSize, sizeof(float));
   for (i=0;i<st->windowSize;i++)
      st->window[i]=.5*(1-cos(2*M_PI*i/st->windowSize));

   st->lagWindow = malloc((st->lpcSize+1)*sizeof(float));
   for (i=0;i<st->lpcSize+1;i++)
      st->lagWindow[i]=exp(-.5*sqr(2*M_PI*st->lag_factor*i));

   st->rc = malloc(st->lpcSize*sizeof(float));
   st->autocorr = malloc((st->lpcSize+1)*sizeof(float));
   st->lpc = malloc((st->lpcSize+1)*sizeof(float));
   st->bw_lpc1 = malloc((st->lpcSize+1)*sizeof(float));
   st->bw_lpc2 = malloc((st->lpcSize+1)*sizeof(float));
   st->lsp = malloc(st->lpcSize*sizeof(float));
   st->qlsp = malloc(st->lpcSize*sizeof(float));
   st->old_lsp = malloc(st->lpcSize*sizeof(float));
   st->old_qlsp = malloc(st->lpcSize*sizeof(float));
   st->interp_lsp = malloc(st->lpcSize*sizeof(float));
   st->interp_qlsp = malloc(st->lpcSize*sizeof(float));
   st->interp_lpc = malloc((st->lpcSize+1)*sizeof(float));
   st->interp_qlpc = malloc((st->lpcSize+1)*sizeof(float));

   st->mem_sp = calloc(st->lpcSize, sizeof(float));
   st->mem_sp2 = calloc(st->lpcSize, sizeof(float));
   st->mem_sw = calloc(st->lpcSize, sizeof(float));

}

void sb_encoder_destroy(SBEncState *st)
{
   encoder_destroy(&st->st_low);
   free(st->x0);
   free(st->x0d);
   free(st->x1);
   free(st->x1d);
   free(st->high);
   free(st->y0);
   free(st->y1);
   free(st->h0_mem);
   free(st->h1_mem);
   free(st->g0_mem);
   free(st->g1_mem);
   
   free(st->buf);
   free(st->window);
   free(st->excBuf);
   free(st->sw);
   free(st->res);
   free(st->target);
   free(st->lagWindow);
   free(st->rc);
   free(st->autocorr);
   free(st->lpc);
   free(st->bw_lpc1);
   free(st->bw_lpc2);
   free(st->lsp);
   free(st->qlsp);
   free(st->old_lsp);
   free(st->old_qlsp);
   free(st->interp_lsp);
   free(st->interp_qlsp);
   free(st->interp_lpc);
   free(st->interp_qlpc);

   free(st->mem_sp);
   free(st->mem_sp2);
   free(st->mem_sw);

   free(st->stack);
   
}


void sb_encode(SBEncState *st, float *in, FrameBits *bits)
{
   int i, roots, sub;
   /* Compute the two sub-bands by filtering with h0 and h1*/
   fir_mem(in, h0, st->x0, st->full_frame_size, QMF_ORDER, st->h0_mem);
   fir_mem(in, h1, st->x1, st->full_frame_size, QMF_ORDER, st->h1_mem);
   /* Down-sample x0 and x1 */
   for (i=0;i<st->frame_size;i++)
   {
      st->x0d[i]=st->x0[i<<1];
      st->x1d[i]=st->x1[i<<1];
   }
   /* Encode the narrowband part*/
   encode(&st->st_low, st->x0d, bits);

   /* High-band buffering / sync with low band */
   for (i=0;i<st->frame_size;i++)
   {
      st->excBuf[i]=st->exc[i];
      st->high[i]=st->high[st->frame_size+i];
      st->high[st->frame_size+i]=st->x1d[i];
   }
   
   /* Start encoding the high-band */

   for (i=0;i<st->windowSize;i++)
      st->buf[i] = st->high[i] * st->window[i];

   /* Compute auto-correlation */
   autocorr(st->buf, st->autocorr, st->lpcSize+1, st->windowSize);

   st->autocorr[0] += 1;        /* prevents NANs */
   st->autocorr[0] *= st->lpc_floor; /* Noise floor in auto-correlation domain */
   /* Lag windowing: equivalent to filtering in the power-spectrum domain */
   for (i=0;i<st->lpcSize+1;i++)
      st->autocorr[i] *= st->lagWindow[i];

   /* Levinson-Durbin */
   wld(st->lpc+1, st->autocorr, st->rc, st->lpcSize);
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

   /* LSP quantization */
   lsp_quant_high(st->lsp, st->qlsp, st->lpcSize, bits);
   
   /*printf ("high_lsp:");
   for (i=0;i<st->lpcSize;i++)
      printf (" %f", st->lsp[i]);
      printf ("\n");
   for (i=0;i<st->lpcSize;i++)
   st->qlsp[i]=st->lsp[i];
   */

   if (st->first)
   {
      for (i=0;i<st->lpcSize;i++)
         st->old_lsp[i] = st->lsp[i];
      for (i=0;i<st->lpcSize;i++)
         st->old_qlsp[i] = st->qlsp[i];
   }
   
   for (sub=0;sub<st->nbSubframes;sub++)
   {
      float *exc, *sp, *mem, *res, *target, *sw, tmp, filter_ratio;
      int offset;
      
      offset = st->subframeSize*sub;
      sp=st->high+offset;
      exc=st->excBuf+offset;
      res=st->res+offset;
      target=st->target+offset;
      sw=st->sw+offset;

      mem=PUSH(st->stack, st->lpcSize);
      
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

      bw_lpc(st->gamma1, st->interp_lpc, st->bw_lpc1, st->lpcSize);
      bw_lpc(st->gamma2, st->interp_lpc, st->bw_lpc2, st->lpcSize);

      {
         float tmp, rl=0, rh=0;
         tmp=1;
         for (i=0;i<=st->lpcSize;i++)
         {
            rh += tmp*st->interp_qlpc[i];
            tmp = -tmp;
         }
         rl = st->st_low.pi_gain[sub];
         rl=1/(fabs(rl)+.001);
         rh=1/(fabs(rh)+.001);
         filter_ratio=fabs(.001+rh)/(.001+fabs(rl));
         printf ("filter_ratio: %f\n", filter_ratio);
      }
#if 0 /* 1 for spectral folding excitation, 0 for stochastic */
      for (i=0;i<st->lpcSize;i++)
         mem[i]=st->mem_sp[i];
      residue_mem(sp, st->interp_qlpc, exc, st->subframeSize, st->lpcSize, st->mem_sp);
      {
         float el=0,eh=0,g;
         for (i=0;i<st->subframeSize;i++)
            eh+=sqr(exc[i]);
         for (i=0;i<st->subframeSize;i++)
           el+=sqr(st->st_low.exc[offset+i]);
         g=eh/(.01+el);
         printf ("ratio: %f\n", g);
         g=sqrt(g);
         for (i=0;i<st->subframeSize;i++)
           exc[i]=g*st->st_low.exc[offset+i];
      }
      syn_filt_mem(exc, st->interp_qlpc, sp, st->subframeSize, st->lpcSize, mem);

#else
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

      /* Compute target signal */
      for (i=0;i<st->subframeSize;i++)
         target[i]=sw[i]-res[i];
      {
#if 0
         
         float el=0,eh=0,g;
         residue_mem(sp, st->interp_qlpc, exc, st->subframeSize, st->lpcSize, st->mem_sp2);
         
         for (i=0;i<st->subframeSize;i++)
            eh+=sqr(exc[i]);
         overlap_cb_search(target, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2,
                           &stoc[0], 512, &gain, &ind, st->lpcSize,
                           st->subframeSize);
         for (i=0;i<st->subframeSize;i++)
            exc[i]=gain*stoc[ind+i];
         for (i=0;i<st->subframeSize;i++)
            el+=sqr(exc[i]);
         g=sqrt(eh/(el+.001));
         for (i=0;i<st->subframeSize;i++)
            exc[i]*=g;
         
#else
         int k,N=4;
         float el=0,eh=0,g;
         int *index;
         float *gains;
         gains = PUSH(st->stack, N);
         index = (int*) PUSH(st->stack, N);
         residue_mem(sp, st->interp_qlpc, exc, st->subframeSize, st->lpcSize, st->mem_sp2);
         
         for (i=0;i<st->subframeSize;i++)
            eh+=sqr(exc[i]);

         for (i=0;i<st->subframeSize;i++)
            exc[i]=0;
         for (k=0;k<N;k++)
         {
            int of=k*st->subframeSize/N;
         overlap_cb_search(target+of, st->interp_qlpc, st->bw_lpc1, st->bw_lpc2,
                           &stoc[0], 64, &gains[k], &index[k], st->lpcSize,
                           st->subframeSize/N);

         frame_bits_pack(bits,index[k],6);

         for (i=0;i<st->subframeSize;i++)
            res[i]=0;
         for (i=0;i<st->subframeSize/N;i++)
            res[of+i]=gains[k]*stoc[index[k]+i];
         residue_zero(res, st->bw_lpc1, res, st->subframeSize, st->lpcSize);
         syn_filt_zero(res, st->interp_qlpc, res, st->subframeSize, st->lpcSize);
         syn_filt_zero(res, st->bw_lpc2, res, st->subframeSize, st->lpcSize);
         for (i=0;i<st->subframeSize;i++)
            target[i]-=res[i];
         for (i=0;i<st->subframeSize/N;i++)
            exc[of+i]+=gains[k]*stoc[index[k]+i];
         }
         for (i=0;i<st->subframeSize;i++)
            el+=sqr(exc[i]);
         g=sqrt(eh/(el+.001));

         for (i=0;i<st->subframeSize;i++)
            el+=sqr(st->st_low.exc[offset+i]);

         for (k=0;k<N;k++)
         {
            int sign=0;
            float quant, dist, min_dist;
            int best_ind;
            int of=k*st->subframeSize/N;
            gains[k]*=g;

            if (gains[k]<0)
            {
               sign=1;
               gains[k] = -gains[k];
            }
            quant = log((1+gains[k])*filter_ratio/(1+sqrt(el/st->subframeSize)));
            min_dist = sqr(quant-quant_high_gain2[0]);
            best_ind=0;
            for (i=1;i<8;i++)
            {
               dist = sqr(quant-quant_high_gain2[i]);
               if (dist<min_dist)
               {
                  best_ind=i;
                  min_dist=dist;
               }
            }
            quant=quant_high_gain2[best_ind];

            frame_bits_pack(bits,sign,1);
            frame_bits_pack(bits,best_ind,3);

            gains[k]=exp(quant)*(1+sqrt(el/st->subframeSize))/filter_ratio;
            if (sign)
               gains[k] = -gains[k];

            for (i=0;i<st->subframeSize/N;i++)
               exc[of+i]=gains[k]*stoc[index[k]+i];
         }
         /*for (i=0;i<st->subframeSize;i++)
           exc[i]*=g;*/
         POP(st->stack);
         POP(st->stack);
#endif
      }

      /*Keep the previous memory*/
      for (i=0;i<st->lpcSize;i++)
         mem[i]=st->mem_sp[i];
      /* Final signal synthesis from excitation */
      syn_filt_mem(exc, st->interp_qlpc, sp, st->subframeSize, st->lpcSize, st->mem_sp);
       
      /* Compute weighted signal again, from synthesized speech (not sure it's the right thing) */
      residue_mem(sp, st->bw_lpc1, sw, st->subframeSize, st->lpcSize, mem);
      syn_filt_mem(sw, st->bw_lpc2, sw, st->subframeSize, st->lpcSize, st->mem_sw);
#endif

      POP(st->stack);
   }

   /* Up-sample coded low-band and high-band*/
   for (i=0;i<st->frame_size;i++)
   {
      st->x0[(i<<1)]=st->x0d[i];
      st->x1[(i<<1)]=st->high[i];
      st->x0[(i<<1)+1]=0;
      st->x1[(i<<1)+1]=0;
   }
   /* Reconstruct the original */
   fir_mem(st->x0, h0, st->y0, st->full_frame_size, QMF_ORDER, st->g0_mem);
   fir_mem(st->x1, h1, st->y1, st->full_frame_size, QMF_ORDER, st->g1_mem);
   for (i=0;i<st->full_frame_size;i++)
      in[i]=2*(st->y0[i]-st->y1[i]);

   for (i=0;i<st->lpcSize;i++)
      st->old_lsp[i] = st->lsp[i];
   for (i=0;i<st->lpcSize;i++)
      st->old_qlsp[i] = st->qlsp[i];

   st->first=0;
}





void sb_decoder_init(SBDecState *st, SpeexMode *mode)
{
   decoder_init(&st->st_low, mode);
   st->full_frame_size = 2*st->st_low.frameSize;
   st->frame_size = st->st_low.frameSize;
   st->subframeSize = 40;
   st->nbSubframes = 4;
   st->lpcSize=8;

   st->first=1;
   st->stack = calloc(10000, sizeof(float));

   st->x0=calloc(st->full_frame_size, sizeof(float));
   st->x1=calloc(st->full_frame_size, sizeof(float));
   st->x0d=calloc(st->frame_size, sizeof(float));
   st->x1d=calloc(st->frame_size, sizeof(float));
   st->high=calloc(st->full_frame_size, sizeof(float));
   st->y0=calloc(st->full_frame_size, sizeof(float));
   st->y1=calloc(st->full_frame_size, sizeof(float));

   st->h0_mem=calloc(QMF_ORDER, sizeof(float));
   st->h1_mem=calloc(QMF_ORDER, sizeof(float));
   st->g0_mem=calloc(QMF_ORDER, sizeof(float));
   st->g1_mem=calloc(QMF_ORDER, sizeof(float));

   st->exc=calloc(st->frame_size, sizeof(float));

   st->qlsp = malloc(st->lpcSize*sizeof(float));
   st->old_qlsp = malloc(st->lpcSize*sizeof(float));
   st->interp_qlsp = malloc(st->lpcSize*sizeof(float));
   st->interp_qlpc = malloc((st->lpcSize+1)*sizeof(float));

   st->mem_sp = calloc(st->lpcSize, sizeof(float));
   st->mem_sw = calloc(st->lpcSize, sizeof(float));

}

void sb_decoder_destroy(SBDecState *st)
{
   decoder_destroy(&st->st_low);
   free(st->x0);
   free(st->x0d);
   free(st->x1);
   free(st->x1d);
   free(st->high);
   free(st->y0);
   free(st->y1);
   free(st->h0_mem);
   free(st->h1_mem);
   free(st->g0_mem);
   free(st->g1_mem);
   
   free(st->exc);
   free(st->qlsp);
   free(st->old_qlsp);
   free(st->interp_qlsp);
   free(st->interp_qlpc);

   free(st->mem_sp);
   free(st->mem_sw);

   free(st->stack);
   
}


void sb_decode(SBDecState *st, FrameBits *bits, float *out)
{
   int i, sub;

   /* Decode the low-band */
   decode(&st->st_low, bits, st->x0d);

   for (i=0;i<st->frame_size;i++)
      st->exc[i]=0;

   lsp_unquant_high(st->qlsp, st->lpcSize, bits);
   
   if (st->first)
   {
      for (i=0;i<st->lpcSize;i++)
         st->old_qlsp[i] = st->qlsp[i];
   }
   
   for (sub=0;sub<st->nbSubframes;sub++)
   {
      float *exc, *sp, tmp, filter_ratio;
      int offset;
      
      offset = st->subframeSize*sub;
      sp=st->high+offset;
      exc=st->exc+offset;
      
      /* LSP interpolation */
      tmp = (.5 + sub)/st->nbSubframes;
      for (i=0;i<st->lpcSize;i++)
         st->interp_qlsp[i] = (1-tmp)*st->old_qlsp[i] + tmp*st->qlsp[i];

      /* LSPs to x-domain */
      for (i=0;i<st->lpcSize;i++)
         st->interp_qlsp[i] = cos(st->interp_qlsp[i]);

      /* LSP to LPC */
      lsp_to_lpc(st->interp_qlsp, st->interp_qlpc, st->lpcSize, st->stack);

      {
         float rl=0, rh=0;
         tmp=1;
         for (i=0;i<=st->lpcSize;i++)
         {
            rh += tmp*st->interp_qlpc[i];
            tmp = -tmp;
         }
         rl = st->st_low.pi_gain[sub];
         rl=1/(fabs(rl)+.001);
         rh=1/(fabs(rh)+.001);
         filter_ratio=fabs(.001+rh)/(.001+fabs(rl));
         printf ("filter_ratio: %f\n", filter_ratio);
      }

      {
         int k,N=4,el=0;
         int *index;
         float *gains;
         int of;
         gains = PUSH(st->stack, N);
         index = (int*) PUSH(st->stack, N);

         for (i=0;i<st->subframeSize;i++)
           el+=sqr(st->st_low.exc[offset+i]);
         
         for (k=0;k<N;k++)
            index[k] = frame_bits_unpack_unsigned(bits,6);

         for (k=0;k<N;k++)
         {
            int sign, gain_ind;
            sign = frame_bits_unpack_unsigned(bits,1);
            gain_ind = frame_bits_unpack_unsigned(bits,3);
            gains[k]=exp(quant_high_gain2[gain_ind])*(1+sqrt(el/st->subframeSize))/filter_ratio;
            if (sign)
               gains[k] =- gains[k];
         }

         for (k=0;k<N;k++)
         {
            of=k*st->subframeSize/N;
            for (i=0;i<st->subframeSize/N;i++)
               exc[of+i]=gains[k]*stoc[index[k]+i];
         }

         POP(st->stack);
         POP(st->stack);
      }

      syn_filt_mem(exc, st->interp_qlpc, sp, st->subframeSize, st->lpcSize, st->mem_sp);

   }

   /* Up-sample coded low-band and high-band*/
   for (i=0;i<st->frame_size;i++)
   {
      st->x0[(i<<1)]=st->x0d[i];
      st->x1[(i<<1)]=st->high[i];
      st->x0[(i<<1)+1]=0;
      st->x1[(i<<1)+1]=0;
   }
   /* Reconstruct the original */
   fir_mem(st->x0, h0, st->y0, st->full_frame_size, QMF_ORDER, st->g0_mem);
   fir_mem(st->x1, h1, st->y1, st->full_frame_size, QMF_ORDER, st->g1_mem);
   for (i=0;i<st->full_frame_size;i++)
      out[i]=2*(st->y0[i]-st->y1[i]);

   for (i=0;i<st->lpcSize;i++)
      st->old_qlsp[i] = st->qlsp[i];

   st->first=0;

}
