/* Copyright (C) 2002 Jean-Marc Valin 
   File: ltp.c
   Lont-Term Prediction functions

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


static float inner_prod(float *a, float *b, int len)
{
  float sum;
  __asm__ __volatile__ (
  "
  push %%eax
  push %%edi
  push %%ecx
  xorps %%xmm3, %%xmm3
  xorps %%xmm4, %%xmm4

  sub $20, %%ecx

.mul20_loop%=:

  movups (%%eax), %%xmm0
  movups (%%edi), %%xmm1
  mulps %%xmm0, %%xmm1

  movups 16(%%eax), %%xmm5
  movups 16(%%edi), %%xmm6
  mulps %%xmm5, %%xmm6
  addps %%xmm1, %%xmm3

  movups 32(%%eax), %%xmm0
  movups 32(%%edi), %%xmm1
  mulps %%xmm0, %%xmm1
  addps %%xmm6, %%xmm4

  movups 48(%%eax), %%xmm5
  movups 48(%%edi), %%xmm6
  mulps %%xmm5, %%xmm6
  addps %%xmm1, %%xmm3

  movups 64(%%eax), %%xmm0
  movups 64(%%edi), %%xmm1
  mulps %%xmm0, %%xmm1
  addps %%xmm6, %%xmm4
  addps %%xmm1, %%xmm3


  add $80, %%eax
  add $80, %%edi

  sub $20,  %%ecx

  jae .mul20_loop%=

  addps %%xmm4, %%xmm3

  movhlps %%xmm3, %%xmm4
  addps %%xmm4, %%xmm3
  movaps %%xmm3, %%xmm4
  shufps $0x55, %%xmm4, %%xmm4
  addss %%xmm4, %%xmm3
  movss %%xmm3, (%%edx)
  
  pop %%ecx
  pop %%edi
  pop %%eax
  "
  : : "a" (a), "D" (b), "c" (len), "d" (&sum) : "memory");
  return sum;
}
