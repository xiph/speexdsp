/* Copyright (C) 2002 Jean-Marc Valin */
/**
   @file misc.h
   @brief Various compatibility routines for Speex
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

#ifndef MISC_H
#define MISC_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifndef SPEEX_VERSION
#define SPEEX_MAJOR_VERSION 1         /**< Major Speex version. */
#define SPEEX_MINOR_VERSION 1         /**< Minor Speex version. */
#define SPEEX_MICRO_VERSION 14        /**< Micro Speex version. */
#define SPEEX_EXTRA_VERSION ""        /**< Extra Speex version. */
#define SPEEX_VERSION "speex-1.2beta2"  /**< Speex version string. */
#endif

/* A couple test to catch stupid option combinations */
#ifdef FIXED_POINT

#ifdef _USE_SSE
#error SSE is only for floating-point
#endif
#if ((defined (ARM4_ASM)||defined (ARM4_ASM)) && defined(BFIN_ASM)) || (defined (ARM4_ASM)&&defined(ARM5E_ASM))
#error Make up your mind. What CPU do you have?
#endif
#ifdef VORBIS_PSYCHO
#error Vorbis-psy model currently not implemented in fixed-point
#endif

#else

#if defined (ARM4_ASM) || defined(ARM5E_ASM) || defined(BFIN_ASM)
#error I suppose you can have a [ARM4/ARM5E/Blackfin] that has float instructions?
#endif
#ifdef FIXED_POINT_DEBUG
#error "Don't you think enabling fixed-point is a good thing to do if you want to debug that?"
#endif


#endif

#include "arch.h"

/** Convert little endian */
static inline spx_int32_t le_int(spx_int32_t i)
{
#if !defined(__LITTLE_ENDIAN__) && ( defined(WORDS_BIGENDIAN) || defined(__BIG_ENDIAN__) )
   spx_uint32_t ui, ret;
   ui = i;
   ret =  ui>>24;
   ret |= (ui>>8)&0x0000ff00;
   ret |= (ui<<8)&0x00ff0000;
   ret |= (ui<<24);
   return ret;
#else
   return i;
#endif
}

#define speex_fatal(str) _speex_fatal(str, __FILE__, __LINE__);
#define speex_assert(cond) {if (!(cond)) {speex_fatal("assertion failed: " #cond);}}

#ifndef RELEASE
static inline void print_vec(float *vec, int len, char *name)
{
   int i;
   printf ("%s ", name);
   for (i=0;i<len;i++)
      printf (" %f", vec[i]);
   printf ("\n");
}
#endif

#ifdef FIXED_DEBUG
long long spx_mips=0;
#endif


/** Speex wrapper for calloc. To do your own dynamic allocation, all you need to do is replace this function, speex_realloc and speex_free */
#ifndef OVERRIDE_SPEEX_ALLOC
static inline void *speex_alloc (int size)
{
   return calloc(size,1);
}
#endif

/** Same as speex_alloc, except that the area is only needed inside a Speex call (might cause problem with wideband though) */
#ifndef OVERRIDE_SPEEX_ALLOC_SCRATCH
static inline void *speex_alloc_scratch (int size)
{
   return calloc(size,1);
}
#endif

/** Speex wrapper for realloc. To do your own dynamic allocation, all you need to do is replace this function, speex_alloc and speex_free */
#ifndef OVERRIDE_SPEEX_REALLOC
static inline void *speex_realloc (void *ptr, int size)
{
   return realloc(ptr, size);
}
#endif

/** Speex wrapper for calloc. To do your own dynamic allocation, all you need to do is replace this function, speex_realloc and speex_alloc */
#ifndef OVERRIDE_SPEEX_FREE
static inline void speex_free (void *ptr)
{
   free(ptr);
}
#endif

/** Same as speex_free, except that the area is only needed inside a Speex call (might cause problem with wideband though) */
#ifndef OVERRIDE_SPEEX_FREE_SCRATCH
static inline void speex_free_scratch (void *ptr)
{
   free(ptr);
}
#endif

/** Print warning message with integer argument to stderr */
#ifndef OVERRIDE_SPEEX_MOVE
static inline void *speex_move (void *dest, void *src, int n)
{
   return memmove(dest,src,n);
}
#endif

#ifndef OVERRIDE_SPEEX_FATAL
static inline void _speex_fatal(const char *str, const char *file, int line)
{
   fprintf (stderr, "Fatal (internal) error in %s, line %d: %s\n", file, line, str);
   exit(1);
}
#endif

#ifndef OVERRIDE_SPEEX_WARNING
static inline void speex_warning(const char *str)
{
#ifndef DISABLE_WARNINGS
   fprintf (stderr, "warning: %s\n", str);
#endif
}
#endif

#ifndef OVERRIDE_SPEEX_WARNING_INT
static inline void speex_warning_int(const char *str, int val)
{
#ifndef DISABLE_WARNINGS
   fprintf (stderr, "warning: %s %d\n", str, val);
#endif
}
#endif

#ifndef OVERRIDE_SPEEX_NOTIFY
static inline void speex_notify(const char *str)
{
#ifndef DISABLE_NOTIFICATIONS
   fprintf (stderr, "notification: %s\n", str);
#endif
}
#endif

#ifdef FIXED_POINT
/** Generate a pseudo-random number */
static inline spx_word16_t speex_rand(spx_word16_t std, spx_int32_t *seed)
{
   spx_word32_t res;
   *seed = 1664525 * *seed + 1013904223;
   res = MULT16_16(EXTRACT16(SHR32(*seed,16)),std);
   return EXTRACT16(PSHR32(SUB32(res, SHR32(res, 3)),14));
}
#else
/** Generate a pseudo-random number */
static inline spx_word16_t speex_rand(spx_word16_t std, spx_int32_t *seed)
{
   const unsigned int jflone = 0x3f800000;
   const unsigned int jflmsk = 0x007fffff;
   union {int i; float f;} ran;
   *seed = 1664525 * *seed + 1013904223;
   ran.i = jflone | (jflmsk & *seed);
   ran.f -= 1.5;
   return 3.4642*std*ran.f;
}
#endif

#ifndef OVERRIDE_SPEEX_PUTC
/** Speex wrapper for putc */
static inline void _speex_putc(int ch, void *file)
{
   FILE *f = (FILE *)file;
   fprintf(f, "%c", ch);
}
#endif

#endif
