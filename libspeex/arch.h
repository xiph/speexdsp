/* Copyright (C) 2003 Jean-Marc Valin */
/**
   @file arch.h
   @brief Various architecture definitions Speex
*/
/*
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

#ifndef ARCH_H
#define ARCH_H

#define ABS(x) ((x) < 0 ? (-(x)) : (x))

#ifdef FIXED_POINT

typedef signed short spx_word16_t;
typedef signed int   spx_word32_t;
typedef long long    spx_word64_t;
typedef spx_word32_t spx_mem_t;
typedef spx_word16_t spx_coef_t;
typedef spx_word16_t spx_lsp_t;
typedef spx_word32_t spx_sig_t;

#define LPC_SCALING  8192
#define SIG_SCALING  16384
#define LSP_SCALING  8192.
#define GAMMA_SCALING 32768.
#define GAIN_SCALING 64
#define GAIN_SCALING_1 0.015625

#define LPC_SHIFT    13
#define SIG_SHIFT    14

#define VERY_SMALL 0

#ifdef COUNT_MIPS
extern long long spx_mips;
#define MIPS_INC spx_mips++,
#else
#define MIPS_INC
#endif

#define SHR(a,shift) ((a) >> (shift))
#define SHL(a,shift) ((a) << (shift))


#define ADD16(a,b) (MIPS_INC(short)((short)(a)+(short)(b)))
#define SUB16(a,b) (MIPS_INC(a)-(b))
#define ADD32(a,b) (MIPS_INC(a)+(b))
#define SUB32(a,b) (MIPS_INC(a)-(b))
#define ADD64(a,b) (MIPS_INC(a)+(b))

#define PSHR(a,shift) (SHR((a)+(1<<((shift)-1)),shift))

/* result fits in 16 bits */
#define MULT16_16_16(a,b)     (MIPS_INC(short)(((short)(a))*((short)(b))))
/* Kludge: just making sure results are on 32 bits */
#ifdef ARM_ASM
static inline spx_word32_t MULT16_16(spx_word16_t x, spx_word16_t y) {
  int res;
  asm volatile("smulbb  %0,%1,%2;\n"
              : "=&r"(res)
              : "%r"(x),"r"(y));
  return(MIPS_INC res);
}
#else
#define MULT16_16(a,b)     (MIPS_INC((short)(a))*((short)(b)))
#endif



#ifdef ARM_ASM

static inline spx_word32_t MAC16_16(spx_word32_t a, spx_word16_t x, spx_word32_t y) {
  int res;
  asm volatile("smlabb  %0,%1,%2,%3;\n"
              : "=&r"(res)
               : "%r"(x),"r"(y),"r"(a));
  return(MIPS_INC res);
}

#else

#define MAC16_16(c,a,b)     (ADD32((c),MULT16_16((a),(b))))

#endif

#define MULT16_32_Q12(a,b) ADD32(MULT16_16((a),SHR((b),12)), SHR(MULT16_16((a),((b)&0x00000fff)),12))
#define MULT16_32_Q13(a,b) ADD32(MULT16_16((a),SHR((b),13)), SHR(MULT16_16((a),((b)&0x00001fff)),13))
#define MULT16_32_Q14(a,b) ADD32(MULT16_16((a),SHR((b),14)), SHR(MULT16_16((a),((b)&0x00003fff)),14))

#ifdef ARM_ASM
static inline spx_word32_t MULT16_32_Q15(spx_word16_t x, spx_word32_t y) {
  int res;
  asm volatile("smulwb  %0,%1,%2;\n"
              : "=&r"(res)
               : "%r"(y<<1),"r"(x));
  return(MIPS_INC res);
}
static inline spx_word32_t MAC16_32_Q15(spx_word32_t a, spx_word16_t x, spx_word32_t y) {
  int res;
  asm volatile("smlawb  %0,%1,%2,%3;\n"
              : "=&r"(res)
               : "%r"(y<<1),"r"(x),"r"(a));
  return(MIPS_INC res);
}
static inline spx_word32_t MULT16_32_Q11(spx_word16_t x, spx_word32_t y) {
  int res;
  asm volatile("smulwb  %0,%1,%2;\n"
              : "=&r"(res)
               : "%r"(y<<5),"r"(x));
  return(MIPS_INC res);
}
static inline spx_word32_t MAC16_32_Q11(spx_word32_t a, spx_word16_t x, spx_word32_t y) {
  int res;
  asm volatile("smlawb  %0,%1,%2,%3;\n"
              : "=&r"(res)
               : "%r"(y<<5),"r"(x),"r"(a));
  return(MIPS_INC res);
}

#else
#define MULT16_32_Q11(a,b) ADD32(MULT16_16((a),SHR((b),11)), SHR(MULT16_16((a),((b)&0x000007ff)),11))
#define MAC16_32_Q11(c,a,b) ADD32(c,ADD32(MULT16_16((a),SHR((b),11)), SHR(MULT16_16((a),((b)&0x000007ff)),11)))

#define MULT16_32_Q15(a,b) ADD32(MULT16_16((a),SHR((b),15)), SHR(MULT16_16((a),((b)&0x00007fff)),15))
#define MAC16_32_Q15(c,a,b) ADD32(c,ADD32(MULT16_16((a),SHR((b),15)), SHR(MULT16_16((a),((b)&0x00007fff)),15)))
#endif


#define MULT16_16_Q13(a,b) (SHR(MULT16_16((a),(b)),13))
#define MULT16_16_Q14(a,b) (SHR(MULT16_16((a),(b)),14))
#define MULT16_16_Q15(a,b) (SHR(MULT16_16((a),(b)),15))

#define MULT16_16_P13(a,b) (SHR(ADD32(4096,MULT16_16((a),(b))),13))
#define MULT16_16_P14(a,b) (SHR(ADD32(8192,MULT16_16((a),(b))),14))
#define MULT16_16_P15(a,b) (SHR(ADD32(16384,MULT16_16((a),(b))),15))

#define MUL_16_32_R15(a,bh,bl) ADD32(MULT16_16((a),(bh)), SHR(MULT16_16((a),(bl)),15))



#define DIV32_16(a,b) ((short)(((signed int)(a))/((short)(b))))
#define DIV32(a,b) (((signed int)(a))/((signed int)(b)))

#else

typedef float spx_mem_t;
typedef float spx_coef_t;
typedef float spx_lsp_t;
typedef float spx_sig_t;
typedef float spx_word16_t;
typedef float spx_word32_t;
typedef float spx_word64_t;

#define LPC_SCALING  1.
#define SIG_SCALING  1.
#define LSP_SCALING  1.
#define GAMMA_SCALING 1.
#define GAIN_SCALING 1.
#define GAIN_SCALING_1 1.

#define LPC_SHIFT    0
#define SIG_SHIFT    0

#define VERY_SMALL 1e-30

#define PSHR(a,shift)       (a)
#define SHR(a,shift)       (a)
#define SHL(a,shift)       (a)
#define ADD16(a,b) ((a)+(b))
#define SUB16(a,b) ((a)-(b))
#define ADD32(a,b) ((a)+(b))
#define SUB32(a,b) ((a)-(b))
#define ADD64(a,b) ((a)+(b))
#define MULT16_16_16(a,b)     ((a)*(b))
#define MULT16_16(a,b)     ((a)*(b))
#define MAC16_16(c,a,b)     ((c)+(a)*(b))

#define MULT16_32_Q11(a,b)     ((a)*(b))
#define MULT16_32_Q13(a,b)     ((a)*(b))
#define MULT16_32_Q14(a,b)     ((a)*(b))
#define MULT16_32_Q15(a,b)     ((a)*(b))

#define MAC16_32_Q11(c,a,b)     ((c)+(a)*(b))
#define MAC16_32_Q15(c,a,b)     ((c)+(a)*(b))

#define MULT16_16_Q13(a,b)     ((a)*(b))
#define MULT16_16_Q14(a,b)     ((a)*(b))
#define MULT16_16_Q15(a,b)     ((a)*(b))

#define DIV32_16(a,b)     ((a)/(b))
#define DIV32(a,b)     ((a)/(b))


#endif



#endif
