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
   encoder_init(&st->st_low, mode);
   st->frame_size = 2*st->st_low.frameSize;
   st->x0=calloc(st->frame_size, sizeof(float));
   st->x1=calloc(st->frame_size, sizeof(float));
   st->x0d=calloc(st->frame_size>>1, sizeof(float));
   st->x1d=calloc(st->frame_size>>1, sizeof(float));
   st->high=calloc(st->frame_size, sizeof(float));
   st->y0=calloc(st->frame_size, sizeof(float));
   st->y1=calloc(st->frame_size, sizeof(float));

   st->h0_mem=calloc(32, sizeof(float));
   st->h1_mem=calloc(32, sizeof(float));
   st->g0_mem=calloc(32, sizeof(float));
   st->g1_mem=calloc(32, sizeof(float));
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
}


void sb_encode(SBEncState *st, float *in, FrameBits *bits)
{
   int i;
   fir_mem(in, h0, st->x0, st->frame_size, 32, st->h0_mem);
   fir_mem(in, h1, st->x1, st->frame_size, 32, st->h1_mem);
   for (i=0;i<st->frame_size>>1;i++)
   {
      st->x0d[i]=st->x0[i<<1];
      st->x1d[i]=st->x1[i<<1];
   }
   encode(&st->st_low, st->x0d, bits);
   for (i=0;i<st->frame_size>>1;i++)
   {
      st->high[i]=st->high[(st->frame_size>>1)+i];
      st->high[(st->frame_size>>1)+i]=st->x1d[i];
   }
   for (i=0;i<st->frame_size>>1;i++)
   {
      st->x0[(i<<1)]=st->x0d[i];
      st->x1[(i<<1)]=st->high[i];
      st->x0[(i<<1)+1]=0;
      st->x1[(i<<1)+1]=0;
   }
   fir_mem(st->x0, h0, st->y0, st->frame_size, 32, st->g0_mem);
   fir_mem(st->x1, h1, st->y1, st->frame_size, 32, st->g1_mem);
   for (i=0;i<st->frame_size;i++)
      in[i]=2*(st->y0[i]-st->y1[i]);
}
