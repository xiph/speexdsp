/* Copyright (C) 2005 Jean-Marc Valin, CSIRO, Christopher Montgomery
   File: vorbis_psy.c

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

/*#define VORBIS_PSYCHO*/

#ifdef VORBIS_PSYCHO

#include "misc.h"
#include "smallft.h"
#include "lpc.h"
#include "vorbis_psy.h"

VorbisPsy *vorbis_psy_init(int rate, int size)
{
   VorbisPsy *psy = speex_alloc(sizeof(VorbisPsy));
   psy->size = size;
   spx_drft_init(&psy->lookup, size);
   return psy;
}

void vorbis_psy_destroy(VorbisPsy *psy)
{
   spx_drft_clear(&psy->lookup);
   speex_free(psy);
}

void compute_curve(VorbisPsy *psy, float *audio, float *curve)
{
   int len = psy->size >> 1;
   int i;
   for (i=0;i<len;i++)
      curve[i] = 1.;
}

/* Transform a masking curve (power spectrum) into a pole-zero filter */
void curve_to_lpc(VorbisPsy *psy, float *curve, float *awk1, float *awk2, int ord)
{
   int i;
   float ac[psy->size];
   int len = psy->size >> 1;
   for (i=0;i<2*len;i++)
      ac[i] = 0;
   for (i=1;i<len;i++)
      ac[2*i-1] = curve[i];
   ac[0] = curve[0];
   ac[2*len-1] = curve[len-1];
   
   spx_drft_backward(&psy->lookup, ac);
   _spx_lpc(awk1, ac, ord);
#if 0
   for (i=0;i<ord;i++)
      awk2[i] = 0;
#else
   /* Use the second (awk2) filter to correct the first one */
   for (i=0;i<2*len;i++)
      ac[i] = 0;   
   for (i=0;i<ord;i++)
      ac[i+1] = awk1[i];
   ac[0] = 1;
   spx_drft_forward(&psy->lookup, ac);
   /* Compute (power) response of awk1 (all zero) */
   ac[0] *= ac[0];
   for (i=1;i<len;i++)
      ac[i] = ac[2*i-1]*ac[2*i-1] + ac[2*i]*ac[2*i];
   ac[len] = ac[2*len-1]*ac[2*len-1];
   /* Compute correction required */
   for (i=0;i<len;i++)
      curve[i] = 1. / (1e-6f+curve[i]*ac[i]);

   for (i=0;i<2*len;i++)
      ac[i] = 0;
   for (i=1;i<len;i++)
      ac[2*i-1] = curve[i];
   ac[0] = curve[0];
   ac[2*len-1] = curve[len-1];
   
   spx_drft_backward(&psy->lookup, ac);
   _spx_lpc(awk2, ac, ord);
   
#endif
}

#if 0
#include <stdio.h>
#include <math.h>

#define ORDER 10
#define CURVE_SIZE 24

int main()
{
   int i;
   float curve[CURVE_SIZE];
   float awk1[ORDER], awk2[ORDER];
   for (i=0;i<CURVE_SIZE;i++)
      scanf("%f ", &curve[i]);
   for (i=0;i<CURVE_SIZE;i++)
      curve[i] = pow(10.f, .1*curve[i]);
   curve_to_lpc(curve, CURVE_SIZE, awk1, awk2, ORDER);
   for (i=0;i<ORDER;i++)
      printf("%f ", awk1[i]);
   printf ("\n");
   for (i=0;i<ORDER;i++)
      printf("%f ", awk2[i]);
   printf ("\n");
   return 0;
}
#endif

#endif
