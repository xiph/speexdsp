/* Copyright (C) 2002 Jean-Marc Valin 
   File: filters.c
   Various analysis/synthesis filters

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

void filter_mem2(float *x, float *_num, float *_den, float *y, int N, int ord, float *_mem)
{
   float __num[20], __den[20], __mem[20];
   float *num, *den, *mem;
   int i;

   num = (float*)(((int)(__num+4))&0xfffffff0)-1;
   den = (float*)(((int)(__den+4))&0xfffffff0)-1;
   mem = (float*)(((int)(__mem+4))&0xfffffff0)-1;
   for (i=0;i<=10;i++)
      num[i]=den[i]=0;
   for (i=0;i<10;i++)
      mem[i]=0;

   for (i=0;i<ord+1;i++)
   {
      num[i]=_num[i];
      den[i]=_den[i];
   }
   for (i=0;i<ord;i++)
      mem[i]=_mem[i];
   for (i=0;i<N;i+=4)
   {

      __asm__ __volatile__ 
      ("
       movss (%1), %%xmm0
       movss (%0), %%xmm1
       addss %%xmm0, %%xmm1
       movss %%xmm1, (%2)
       shufps $0x00, %%xmm0, %%xmm0
       shufps $0x00, %%xmm1, %%xmm1

       movaps 4(%3),  %%xmm2
       movaps 4(%4),  %%xmm3
       mulps  %%xmm0, %%xmm2
       mulps  %%xmm1, %%xmm3
       movaps 20(%3), %%xmm4
       mulps  %%xmm0, %%xmm4
       addps  4(%0),  %%xmm2
       movaps 20(%4), %%xmm5
       mulps  %%xmm1, %%xmm5
       addps  20(%0), %%xmm4
       subps  %%xmm3, %%xmm2
       movups %%xmm2, (%0)
       subps  %%xmm5, %%xmm4
       movups %%xmm4, 16(%0)

       movss  36(%3), %%xmm2
       mulss  %%xmm0, %%xmm2
       movss  36(%4), %%xmm3
       mulss  %%xmm1, %%xmm3
       addss  36(%0), %%xmm2
       movss  40(%3), %%xmm4
       mulss  %%xmm0, %%xmm4
       movss  40(%4), %%xmm5
       mulss  %%xmm1, %%xmm5
       subss  %%xmm3, %%xmm2
       movss  %%xmm2, 32(%0)       
       subss  %%xmm5, %%xmm4
       movss  %%xmm4, 36(%0)



       movss 4(%1), %%xmm0
       movss (%0), %%xmm1
       addss %%xmm0, %%xmm1
       movss %%xmm1, 4(%2)
       shufps $0x00, %%xmm0, %%xmm0
       shufps $0x00, %%xmm1, %%xmm1

       movaps 4(%3),  %%xmm2
       movaps 4(%4),  %%xmm3
       mulps  %%xmm0, %%xmm2
       mulps  %%xmm1, %%xmm3
       movaps 20(%3), %%xmm4
       mulps  %%xmm0, %%xmm4
       addps  4(%0),  %%xmm2
       movaps 20(%4), %%xmm5
       mulps  %%xmm1, %%xmm5
       addps  20(%0), %%xmm4
       subps  %%xmm3, %%xmm2
       movups %%xmm2, (%0)
       subps  %%xmm5, %%xmm4
       movups %%xmm4, 16(%0)

       movss  36(%3), %%xmm2
       mulss  %%xmm0, %%xmm2
       movss  36(%4), %%xmm3
       mulss  %%xmm1, %%xmm3
       addss  36(%0), %%xmm2
       movss  40(%3), %%xmm4
       mulss  %%xmm0, %%xmm4
       movss  40(%4), %%xmm5
       mulss  %%xmm1, %%xmm5
       subss  %%xmm3, %%xmm2
       movss  %%xmm2, 32(%0)       
       subss  %%xmm5, %%xmm4
       movss  %%xmm4, 36(%0)



       movss 8(%1), %%xmm0
       movss (%0), %%xmm1
       addss %%xmm0, %%xmm1
       movss %%xmm1, 8(%2)
       shufps $0x00, %%xmm0, %%xmm0
       shufps $0x00, %%xmm1, %%xmm1

       movaps 4(%3),  %%xmm2
       movaps 4(%4),  %%xmm3
       mulps  %%xmm0, %%xmm2
       mulps  %%xmm1, %%xmm3
       movaps 20(%3), %%xmm4
       mulps  %%xmm0, %%xmm4
       addps  4(%0),  %%xmm2
       movaps 20(%4), %%xmm5
       mulps  %%xmm1, %%xmm5
       addps  20(%0), %%xmm4
       subps  %%xmm3, %%xmm2
       movups %%xmm2, (%0)
       subps  %%xmm5, %%xmm4
       movups %%xmm4, 16(%0)

       movss  36(%3), %%xmm2
       mulss  %%xmm0, %%xmm2
       movss  36(%4), %%xmm3
       mulss  %%xmm1, %%xmm3
       addss  36(%0), %%xmm2
       movss  40(%3), %%xmm4
       mulss  %%xmm0, %%xmm4
       movss  40(%4), %%xmm5
       mulss  %%xmm1, %%xmm5
       subss  %%xmm3, %%xmm2
       movss  %%xmm2, 32(%0)       
       subss  %%xmm5, %%xmm4
       movss  %%xmm4, 36(%0)



       movss 12(%1), %%xmm0
       movss (%0), %%xmm1
       addss %%xmm0, %%xmm1
       movss %%xmm1, 12(%2)
       shufps $0x00, %%xmm0, %%xmm0
       shufps $0x00, %%xmm1, %%xmm1

       movaps 4(%3),  %%xmm2
       movaps 4(%4),  %%xmm3
       mulps  %%xmm0, %%xmm2
       mulps  %%xmm1, %%xmm3
       movaps 20(%3), %%xmm4
       mulps  %%xmm0, %%xmm4
       addps  4(%0),  %%xmm2
       movaps 20(%4), %%xmm5
       mulps  %%xmm1, %%xmm5
       addps  20(%0), %%xmm4
       subps  %%xmm3, %%xmm2
       movups %%xmm2, (%0)
       subps  %%xmm5, %%xmm4
       movups %%xmm4, 16(%0)

       movss  36(%3), %%xmm2
       mulss  %%xmm0, %%xmm2
       movss  36(%4), %%xmm3
       mulss  %%xmm1, %%xmm3
       addss  36(%0), %%xmm2
       movss  40(%3), %%xmm4
       mulss  %%xmm0, %%xmm4
       movss  40(%4), %%xmm5
       mulss  %%xmm1, %%xmm5
       subss  %%xmm3, %%xmm2
       movss  %%xmm2, 32(%0)       
       subss  %%xmm5, %%xmm4
       movss  %%xmm4, 36(%0)

       "
       : : "r" (mem), "r" (x+i), "r" (y+i), "r" (num), "r" (den)
       : "memory" );

   }
   for (i=0;i<ord;i++)
      _mem[i]=mem[i];

}


void iir_mem2(float *x, float *_den, float *y, int N, int ord, float *_mem)
{
   float  __den[20], __mem[20];
   float *den, *mem;
   int i;

   den = (float*)(((int)(__den+4))&0xfffffff0)-1;
   mem = (float*)(((int)(__mem+4))&0xfffffff0)-1;
   for (i=0;i<=10;i++)
      den[i]=0;
   for (i=0;i<10;i++)
      mem[i]=0;
   for (i=0;i<ord+1;i++)
   {
      den[i]=_den[i];
   }
   for (i=0;i<ord;i++)
      mem[i]=_mem[i];

   for (i=0;i<N;i++)
   {
#if 0
      y[i] = x[i] + mem[0];
      for (j=0;j<ord-1;j++)
      {
         mem[j] = mem[j+1] - den[j+1]*y[i];
      }
      mem[ord-1] = - den[ord]*y[i];
#else
      __asm__ __volatile__ 
      ("
       movss (%1), %%xmm0
       movss (%0), %%xmm1
       addss %%xmm0, %%xmm1
       movss %%xmm1, (%2)
       shufps $0x00, %%xmm0, %%xmm0
       shufps $0x00, %%xmm1, %%xmm1

       
       movaps 4(%3),  %%xmm2
       movaps 20(%3), %%xmm3
       mulps  %%xmm1, %%xmm2
       mulps  %%xmm1, %%xmm3
       movss  36(%3), %%xmm4
       movss  40(%3), %%xmm5
       mulss  %%xmm1, %%xmm4
       mulss  %%xmm1, %%xmm5
       movaps 4(%0),  %%xmm6
       subps  %%xmm2, %%xmm6
       movups %%xmm6, (%0)
       movaps 20(%0), %%xmm7
       subps  %%xmm3, %%xmm7
       movups %%xmm7, 16(%0)


       movss  36(%0), %%xmm7
       subss  %%xmm4, %%xmm7
       movss  %%xmm7, 32(%0)       
       xorps  %%xmm2, %%xmm2
       subss  %%xmm5, %%xmm2
       movss  %%xmm2, 36(%0)

       "
       : : "r" (mem), "r" (x+i), "r" (y+i), "r" (den)
       : "memory" );
#endif
   }
   for (i=0;i<ord;i++)
      _mem[i]=mem[i];

}

