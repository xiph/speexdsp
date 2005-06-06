/* Copyright (C) 2005 Analog Devices
   Author: Jean-Marc Valin 
   File: filters_bfin.h
   Various analysis/synthesis filters (Blackfin version)

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

#include <stdio.h>

#define OVERRIDE_NORMALIZE16
int normalize16(const spx_sig_t *x, spx_word16_t *y, spx_sig_t max_scale, int len)
{
   spx_sig_t max_val=1;
   int sig_shift;

   __asm__ 
   (
   "%0 = 0;\n\t"
   "I0 = %1;\n\t"
   "L0 = 0;\n\t"
   "LOOP norm_max%= LC0 = %2;\n\t"
   "LOOP_BEGIN norm_max%=;\n\t"
      "R1 = [I0++];\n\t"
      "R1 = ABS R1;\n\t"
      "%0 = MAX(%0, R1);\n\t"
   "LOOP_END norm_max%=;\n\t"
   : "=&d" (max_val)
   : "a" (x), "a" (len)
   : "R1"
   );

   sig_shift=0;
   while (max_val>max_scale)
   {
      sig_shift++;
      max_val >>= 1;
   }

   __asm__ __volatile__ 
   (
   "I0 = %0;\n\t"
   "L0 = 0;\n\t"
   "I1 = %1;\n\t"
   "L1 = 0;\n\t"
   "R0 = [I0++];\n\t"
   "LOOP norm_shift%= LC0 = %3 >> 1;\n\t"
   "LOOP_BEGIN norm_shift%=;\n\t"
      "R1 = ASHIFT R0 by %2.L || R2 = [I0++];\n\t"
      "R3 = ASHIFT R2 by %2.L || R0 = [I0++];\n\t"
      "R3 = PACK(R3.L, R1.L);\n\t"
      "[I1++] = R3;\n\t"
   "LOOP_END norm_shift%=;\n\t"
   : : "a" (x), "a" (y), "d" (-sig_shift), "a" (len)
   : "I0", "L0", "I1", "L1", "R0", "R1", "R2", "R3", "memory"
   );
   return sig_shift;
}

#define OVERRIDE_FILTER_MEM2
void filter_mem2(const spx_sig_t *_x, const spx_coef_t *num, const spx_coef_t *den, spx_sig_t *_y, int N, int ord, spx_mem_t *mem)
{
   spx_word16_t x[N],y[N];
   spx_word16_t *xx, *yy;
   xx = x;
   yy = y;
   
   __asm__ __volatile__
   (
   /* Register setup */
   "R0 = %7;\n\t"      /*ord */
   
   "P0 = %4;\n\t"
   "I0 = P0;\n\t"
   "B0 = P0;\n\t"
   "L0 = 0;\n\t"
   
   "P1 = %5;\n\t"
   "I1 = P1;\n\t"
   "B1 = P1;\n\t"
   "L1 = 0;\n\t"
   
   "P2 = %0;\n\t"
   "P3 = %1;\n\t"
   "I2 = P2;\n\t"
   "I3 = P3;\n\t"
   "L2 = 0;\n\t"
   "L3 = 0;\n\t"
   
   "P4 = %8;\n\t"
   "P0 = %2;\n\t"
   "P1 = %3;\n\t"
   
   /* First sample */
   "R1 = [P4++];\n\t"
   "R1 <<= 1;\n\t"
   "R2 = [P0++];\n\t"
   "R1 = R1 + R2;\n\t"
   "[P1++] = R1;\n\t"
   "R1 <<= 2;\n\t"
   "W[P3] = R1.H;\n\t"
   "R2 <<= 2;\n\t"
   "W[P2] = R2.H;\n\t"

   /* Samples 1 to ord-1 (using memory) */
   "R0 += -1;\n\t"
   "R3 = 0;\n\t"
   "LC0 = R0;\n\t"
   "LOOP filter_start%= LC0;\n\t"
   "LOOP_BEGIN filter_start%=;\n\t"
      "R3 += 1;\n\t"
      "LC1 = R3;\n\t"
      
      "R1 = [P4++];\n\t"
      "A1 = R1;\n\t"
      "I0 = B0;\n\t"
      "I1 = B1;\n\t"
      "I2 = P2;\n\t"
      "I3 = P3;\n\t"
      "P2 += 2;\n\t"
      "P3 += 2;\n\t"
      "LOOP filter_start_inner%= LC1;\n\t"
      "LOOP_BEGIN filter_start_inner%=;\n\t"
         "R4.L = W[I0++];\n\t"
         "R5.L = W[I2--];\n\t"
         "A1 += R4.L*R5.L (IS);\n\t"
         "R4.L = W[I1++];\n\t"
         "R5.L = W[I3--];\n\t"
         "A1 -= R4.L*R5.L (IS);\n\t"
      "LOOP_END filter_start_inner%=;\n\t"
   
      "R1 = A1;\n\t"
      "R1 <<= 1;\n\t"
      "R2 = [P0++];\n\t"
      "R1 = R1 + R2;\n\t"
      "[P1++] = R1;\n\t"
      "R1 <<= 2;\n\t"
      "W[P3] = R1.H;\n\t"
      "R2 <<= 2;\n\t"
      "W[P2] = R2.H;\n\t"
   "LOOP_END filter_start%=;\n\t"

   /* Samples ord to N*/   
   "R0 = %7;\n\t"
   "R0 <<= 1;\n\t"
   "I0 = B0;\n\t"
   "I1 = B1;\n\t"
   "L0 = R0;\n\t"
   "L1 = R0;\n\t"
   
   "R0 = %7;\n\t"
   "R2 = %6;\n\t"
   "R2 = R2 - R0;\n\t"
   "R4.L = W[I0++];\n\t"
   "R6.L = W[I1++];\n\t"
   "LC0 = R2;\n\t"
   "LOOP filter_mid%= LC0;\n\t"
   "LOOP_BEGIN filter_mid%=;\n\t"
      "LC1 = R0;\n\t"
      "A1 = 0;\n\t"

      "I2 = P2;\n\t"
      "I3 = P3;\n\t"
      "P2 += 2;\n\t"
      "P3 += 2;\n\t"
      "R5.L = W[I2--];\n\t"
      "R7.L = W[I3--];\n\t"
      "LOOP filter_mid_inner%= LC1;\n\t"
      "LOOP_BEGIN filter_mid_inner%=;\n\t"
         "A1 += R4.L*R5.L (IS) || R4.L = W[I0++] || R5.L = W[I2--];\n\t"
         "A1 -= R6.L*R7.L (IS) || R6.L = W[I1++] || R7.L = W[I3--];\n\t"
      "LOOP_END filter_mid_inner%=;\n\t"
      "R1 = A1;\n\t"
      "R1 <<= 1;\n\t"
      "R2 = [P0++];\n\t"
      "R1 = R1 + R2;\n\t"
      "[P1++] = R1;\n\t"
      "R1 <<= 2;\n\t"
      "W[P3] = R1.H;\n\t"
      "R2 <<= 2;\n\t"
      "W[P2] = R2.H;\n\t"
   "LOOP_END filter_mid%=;\n\t"
     
   /* Update memory */
   "P4 = %8;\n\t"
   "R0 = %7;\n\t"
   "LC0 = R0;\n\t"
   "P0 = B0;\n\t"
   "P1 = B1;\n\t"
   "LOOP mem_update%= LC0;\n\t"
   "LOOP_BEGIN mem_update%=;\n\t"
      "A0 = 0;\n\t"
      "I2 = P2;\n\t"
      "I3 = P3;\n\t"
      "I0 = P0;\n\t"
      "I1 = P1;\n\t"
      "P0 += 2;\n\t"
      "P1 += 2;\n\t"
      "R0 = LC0;LC1=R0;\n\t"
      "R5.L = W[I2--];\n\t"
      "R7.L = W[I3--];\n\t"
      "R4.L = W[I0++];\n\t"
      "R6.L = W[I1++];\n\t"
      "LOOP mem_accum%= LC1;\n\t"
      "LOOP_BEGIN mem_accum%=;\n\t"
         "A0 += R4.L*R5.L (IS) || R4.L = W[I0++] || R5.L = W[I2--];\n\t"
         "A0 -= R6.L*R7.L (IS) || R6.L = W[I1++] || R7.L = W[I3--];\n\t"
      "LOOP_END mem_accum%=;\n\t"
      "R0 = A0;\n\t"
      "[P4++] = R0;\n\t"
   "LOOP_END mem_update%=;\n\t"

   : : "m" (xx), "m" (yy), "m" (_x), "m" (_y), "m" (num), "m" (den), "m" (N), "m" (ord), "m" (mem)
   : "R0", "R1", "R2", "R3", "R4", "R5", "R7", "P0", "P1", "P2", "P3", "P4", "B0", "B1", "I0", "I1", "I2", "I3", "L0", "L1", "L2", "L3", "memory"
   );

}




#define OVERRIDE_IIR_MEM2
void iir_mem2(const spx_sig_t *_x, const spx_coef_t *den, spx_sig_t *_y, int N, int ord, spx_mem_t *mem)
{
   spx_word16_t x[N],y[N];
   spx_word16_t *xx, *yy;
   xx = x;
   yy = y;
   __asm__ __volatile__
   (
   /* Register setup */
   "R0 = %5;\n\t"      /*ord */
   
   "P1 = %3;\n\t"
   "I1 = P1;\n\t"
   "B1 = P1;\n\t"
   "L1 = 0;\n\t"
   
   "P3 = %0;\n\t"
   "I3 = P3;\n\t"
   "L3 = 0;\n\t"
   
   "P4 = %6;\n\t"
   "P0 = %1;\n\t"
   "P1 = %2;\n\t"
   
   /* First sample */
   "R1 = [P4++];\n\t"
   "R1 <<= 1;\n\t"
   "R2 = [P0++];\n\t"
   "R1 = R1 + R2;\n\t"
   "[P1++] = R1;\n\t"
   "R1 <<= 2;\n\t"
   "W[P3] = R1.H;\n\t"
   "R2 <<= 2;\n\t"

   /* Samples 1 to ord-1 (using memory) */
   "R0 += -1;\n\t"
   "R3 = 0;\n\t"
   "LC0 = R0;\n\t"
   "LOOP filter_start%= LC0;\n\t"
   "LOOP_BEGIN filter_start%=;\n\t"
      "R3 += 1;\n\t"
      "LC1 = R3;\n\t"
      
      "R1 = [P4++];\n\t"
      "A1 = R1;\n\t"
      "I1 = B1;\n\t"
      "I3 = P3;\n\t"
      "P3 += 2;\n\t"
      "LOOP filter_start_inner%= LC1;\n\t"
      "LOOP_BEGIN filter_start_inner%=;\n\t"
         "R4.L = W[I1++];\n\t"
         "R5.L = W[I3--];\n\t"
         "A1 -= R4.L*R5.L (IS);\n\t"
      "LOOP_END filter_start_inner%=;\n\t"
   
      "R1 = A1;\n\t"
      "R1 <<= 1;\n\t"
      "R2 = [P0++];\n\t"
      "R1 = R1 + R2;\n\t"
      "[P1++] = R1;\n\t"
      "R1 <<= 2;\n\t"
      "W[P3] = R1.H;\n\t"
      "R2 <<= 2;\n\t"
   "LOOP_END filter_start%=;\n\t"

   /* Samples ord to N*/   
   "R0 = %5;\n\t"
   "R0 <<= 1;\n\t"
   "I1 = B1;\n\t"
   "L1 = R0;\n\t"
   
   "R0 = %5;\n\t"
   "R2 = %4;\n\t"
   "R2 = R2 - R0;\n\t"
   "R6.L = W[I1++];\n\t"
   "LC0 = R2;\n\t"
   "LOOP filter_mid%= LC0;\n\t"
   "LOOP_BEGIN filter_mid%=;\n\t"
      "LC1 = R0;\n\t"
      "A1 = 0;\n\t"
      "I3 = P3;\n\t"
      "P3 += 2;\n\t"
      "R7.L = W[I3--];\n\t"
      "LOOP filter_mid_inner%= LC1;\n\t"
      "LOOP_BEGIN filter_mid_inner%=;\n\t"
         "A1 -= R6.L*R7.L (IS) || R6.L = W[I1++] || R7.L = W[I3--];\n\t"
      "LOOP_END filter_mid_inner%=;\n\t"
      "R1 = A1;\n\t"
      "R1 <<= 1;\n\t"
      "R2 = [P0++];\n\t"
      "R1 = R1 + R2;\n\t"
      "[P1++] = R1;\n\t"
      "R1 <<= 2;\n\t"
      "W[P3] = R1.H;\n\t"
      "R2 <<= 2;\n\t"
   "LOOP_END filter_mid%=;\n\t"
     
   /* Update memory */
   "P4 = %6;\n\t"
   "R0 = %5;\n\t"
   "LC0 = R0;\n\t"
   "P0 = B0;\n\t"
   "P1 = B1;\n\t"
   "LOOP mem_update%= LC0;\n\t"
   "LOOP_BEGIN mem_update%=;\n\t"
      "A0 = 0;\n\t"
      "I3 = P3;\n\t"
      "I1 = P1;\n\t"
      "P0 += 2;\n\t"
      "P1 += 2;\n\t"
      "R0 = LC0;LC1=R0;\n\t"
      "R7.L = W[I3--];\n\t"
      "R6.L = W[I1++];\n\t"
      "LOOP mem_accum%= LC1;\n\t"
      "LOOP_BEGIN mem_accum%=;\n\t"
         "A0 -= R6.L*R7.L (IS) || R6.L = W[I1++] || R7.L = W[I3--];\n\t"
      "LOOP_END mem_accum%=;\n\t"
      "R0 = A0;\n\t"
      "[P4++] = R0;\n\t"
   "LOOP_END mem_update%=;\n\t"

   : : "m" (yy), "m" (_x), "m" (_y), "m" (den), "m" (N), "m" (ord), "m" (mem)
   : "R0", "R1", "R2", "R3", "R4", "R5", "R7", "P0", "P1", "P2", "P3", "P4", "B0", "B1", "I0", "I1", "I2", "I3", "L0", "L1", "L2", "L3", "memory"
   );

}

#define OVERRIDE_FIR_MEM2
void fir_mem2(const spx_sig_t *x, const spx_coef_t *num, spx_sig_t *y, int N, int ord, spx_mem_t *mem)
{
   int i;
   spx_coef_t den2[12];
   spx_coef_t *den;
   den = (spx_coef_t*)((((int)den2)+4)&0xfffffffc);
   for (i=0;i<10;i++)
      den[i] = 0;
   filter_mem2(x, num, den, y, N, ord, mem);
}

#define min(a,b) ((a)<(b) ? (a):(b))

#define OVERRIDE_COMPUTE_IMPULSE_RESPONSE
void compute_impulse_response(const spx_coef_t *ak, const spx_coef_t *awk1, const spx_coef_t *awk2, spx_word16_t *y, int N, int ord, char *stack)
{
   int i,j;
   VARDECL(spx_word16_t *ytmp);
   ALLOC(ytmp, N, spx_word16_t);
   
   y[0] = LPC_SCALING;
   for (i=0;i<ord;i++)
      y[i+1] = awk1[i];
   i++;
   for (;i<N;i++)
      y[i] = 0;

   for (i=0;i<N;i++)
   {
      spx_word32_t yi = SHL32(EXTEND32(y[i]),LPC_SHIFT);
      spx_word32_t yi2 = 0;
      for (j=0;j<min(i,ord);j++)
      {
         yi = MAC16_16(yi, awk2[j], -ytmp[i-j-1]);
         yi2 = MAC16_16(yi2, ak[j], -y[i-j-1]);
      }
      ytmp[i] = EXTRACT16(SHR32(yi,LPC_SHIFT));
      yi2 = ADD32(yi2,SHL32(yi,1));
      y[i] = EXTRACT16(SHR32(yi2,LPC_SHIFT));
   }

}



#if 0 /* Equivalent C function for filter_mem2 */
void filter_mem2(const spx_sig_t *_x, const spx_coef_t *num, const spx_coef_t *den, spx_sig_t *_y, int N, int ord, spx_mem_t *mem)
{
   int i,j;
   spx_word16_t xi,yi,nyi;
   spx_word16_t x[N],y[N];
   spx_word16_t *xx, *yy;
   xx = x;
   yy = y;
   for (i=0;i<N;i++)
   {
      x[i] = EXTRACT16(SHR32(_x[i],SIG_SHIFT));
   }
   
   for (i=0;i<ord;i++)
   {
      spx_word32_t yi = mem[i];
      for (j=0;j<i;j++)
      {
         yi = MAC16_16(yi, num[j], x[i-j-1]);
         yi = MAC16_16(yi, den[j], -y[i-j-1]);
      }
      _y[i] = ADD32(_x[i],SHL32(yi,1));
      y[i] = EXTRACT16(SHR32(_y[i],SIG_SHIFT));
   }
   for (i=ord;i<N;i++)
   {
      spx_word32_t yi = 0;
      for (j=0;j<ord;j++)
      {
         yi = MAC16_16(yi, num[j], x[i-j-1]);
         yi = MAC16_16(yi, den[j], -y[i-j-1]);
      }
      _y[i] = ADD32(_x[i],SHL32(yi,1));
      y[i] = EXTRACT16(SHR32(_y[i],SIG_SHIFT));
   }

   for (i=0;i<ord;i++)
   {
      spx_mem_t m = 0;
      for (j=0;j<ord-i;j++)
      {
         m = MAC16_16(m, x[N-1-j], num[j+i]);
         m = MAC16_16(m, -y[N-1-j], den[j+i]);
      }
      mem[i] = m;
   }
}
#endif
