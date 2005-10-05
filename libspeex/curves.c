/* Copyright (C) 2005 Jean-Marc Valin / CSIRO
   File: curves

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


#include "smallft.h"
#include "lpc.h"

struct drft_lookup lookup;

/* FIXME: This is an horrible kludge */
static void fft_init(int size)
{
   static initialized = -1;
   if (size != initialized)
   {
      if (initialized != -1)
         spx_drft_clear(&lookup);
      spx_drft_init(&lookup, size);
      initialized = size;
   }
}
/* */
void curve_to_lpc(float *curve, int len, float *awk1, float *awk2, int ord)
{
   int i;
   float ac[len*2];
   for (i=0;i<2*len;i++)
      ac[i] = 0;
   for (i=1;i<len;i++)
      ac[2*i-1] = curve[i];
   ac[0] = curve[0];
   ac[2*len-1] = curve[len-1];
   
   fft_init(2*len);
   spx_drft_backward(&lookup, ac);
   _spx_lpc(awk1, ac, ord);
   for (i=0;i<ord;i++)
      awk2[i] = 0;
}


#endif
