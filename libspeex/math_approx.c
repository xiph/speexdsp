/* Copyright (C) 2002 Jean-Marc Valin 
   File: math_approx.c
   Various math approximation functions for Speex

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

#include "math_approx.h"
#include "misc.h"

#ifdef FIXED_POINT

#define C0 3634
#define C1 21173
#define C2 -12627
#define C3 4215

spx_word16_t sqroot(spx_word32_t x)
{
   int k=0;
   spx_word32_t rt;

#if 1
   if (x>16777216)
   {
      x>>=10;
      k+=5;
   }
   if (x>1048576)
   {
      x>>=6;
      k+=3;
   }
   if (x>262144)
   {
      x>>=4;
      k+=2;
   }
   if (x>32768)
   {
      x>>=2;
      k+=1;
   }
   if (x>16384)
   {
      x>>=2;
      k+=1;
   }
#else
   while (x>16384)
   {
      x>>=2;
      k++;
      }
#endif
   rt = ADD16(C0, MULT16_16_Q14(x, ADD16(C1, MULT16_16_Q14(x, ADD16(C2, MULT16_16_Q14(x, (C3)))))));
   rt <<= k;
   rt >>=7;
   return rt;
}

#endif
