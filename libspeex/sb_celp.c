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


#include "speex.h"
#include "sb_celp.h"
#include "stdlib.h"
#include "filters.h"
#include <math.h>
#include "lpc.h"
#include "lsp.h"
#include <stdio.h>


#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif

#define sqr(x) ((x)*(x))

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

   st->lag_factor = .01;
   st->lpc_floor = 1.001;
   st->first=1;
   st->stack = calloc(2000, sizeof(float));

   st->x0=calloc(st->full_frame_size, sizeof(float));
   st->x1=calloc(st->full_frame_size, sizeof(float));
   st->x0d=calloc(st->frame_size, sizeof(float));
   st->x1d=calloc(st->frame_size, sizeof(float));
   st->high=calloc(st->full_frame_size, sizeof(float));
   st->y0=calloc(st->full_frame_size, sizeof(float));
   st->y1=calloc(st->full_frame_size, sizeof(float));

   st->h0_mem=calloc(32, sizeof(float));
   st->h1_mem=calloc(32, sizeof(float));
   st->g0_mem=calloc(32, sizeof(float));
   st->g1_mem=calloc(32, sizeof(float));

   st->buf=calloc(st->windowSize, sizeof(float));
   st->window=calloc(st->windowSize, sizeof(float));
   for (i=0;i<st->windowSize;i++)
      st->window[i]=.5*(1-cos(2*M_PI*i/st->windowSize));

   st->lagWindow = malloc((st->lpcSize+1)*sizeof(float));
   for (i=0;i<st->lpcSize+1;i++)
      st->lagWindow[i]=exp(-.5*sqr(2*M_PI*st->lag_factor*i));

   st->rc = malloc(st->lpcSize*sizeof(float));
   st->autocorr = malloc((st->lpcSize+1)*sizeof(float));
   st->lpc = malloc((st->lpcSize+1)*sizeof(float));
   st->lsp = malloc(st->lpcSize*sizeof(float));
   st->qlsp = malloc(st->lpcSize*sizeof(float));
   st->old_lsp = malloc(st->lpcSize*sizeof(float));
   st->old_qlsp = malloc(st->lpcSize*sizeof(float));
   st->interp_lsp = malloc(st->lpcSize*sizeof(float));
   st->interp_qlsp = malloc(st->lpcSize*sizeof(float));
   st->interp_lpc = malloc(st->lpcSize*sizeof(float));
   st->interp_qlpc = malloc(st->lpcSize*sizeof(float));


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
   free(st->lagWindow);
   free(st->rc);
   free(st->autocorr);
   free(st->lpc);
   free(st->lsp);
   free(st->qlsp);
   free(st->old_lsp);
   free(st->old_qlsp);
   free(st->interp_lsp);
   free(st->interp_qlsp);
   free(st->interp_lpc);
   free(st->interp_qlpc);

   free(st->stack);
}


void sb_encode(SBEncState *st, float *in, FrameBits *bits)
{
   int i, roots, sub;
   /* Compute the two sub-bands by filtering with h0 and h1*/
   fir_mem(in, h0, st->x0, st->full_frame_size, 32, st->h0_mem);
   fir_mem(in, h1, st->x1, st->full_frame_size, 32, st->h1_mem);
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
      st->high[i]=st->high[(st->frame_size)+i];
      st->high[(st->frame_size)+i]=st->x1d[i];
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
   
   /* FIXME: Need to really quantize the LSPs*/
   for (i=0;i<st->lpcSize;i++)
      st->qlsp[i]=st->lsp[i];
   if (st->first)
   {
      for (i=0;i<st->lpcSize;i++)
         st->old_lsp[i] = st->lsp[i];
      for (i=0;i<st->lpcSize;i++)
         st->old_qlsp[i] = st->qlsp[i];
   }
   
   for (sub=0;sub<st->nbSubframes;sub++)
   {
      float *exc, *sp, tmp;
      int offset;
      
      offset = st->subframeSize*sub;
      sp=st->high+offset;
      exc=st->exc+offset;

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
   fir_mem(st->x0, h0, st->y0, st->full_frame_size, 32, st->g0_mem);
   fir_mem(st->x1, h1, st->y1, st->full_frame_size, 32, st->g1_mem);
   for (i=0;i<st->full_frame_size;i++)
      in[i]=2*(st->y0[i]-st->y1[i]);

   st->first=0;
}
