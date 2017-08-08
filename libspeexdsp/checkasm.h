/* Copyright (C) 2017 Tristan Matthews */
/**
   @file checkasm.h
   @brief ASM validation functions
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

#ifndef RESAMPLE_CHECKASM_H
#define RESAMPLE_CHECKASM_H

#include <math.h>
#include <stdio.h>
#include "os_support.h"

#define RESAMPLE_CHECK_INNER_PRODUCT(a, b, len, actual, expected, fmt, epsilon) \
   do {                                                \
     int k;                                            \
     double sum_magnitudes = 0.0;                      \
     for (k = 0; k < len; k++) {                       \
       sum_magnitudes += a[k]*a[k] + b[k]*b[k];        \
       expected += a[k]*b[k];                          \
     }                                                 \
     double normalized_error = fabs(expected - actual)/sum_magnitudes; \
     if (normalized_error > epsilon) {                 \
       fprintf(stderr, "ASM mismatch: Error:"fmt" exceeds:"fmt"\n", normalized_error, epsilon); \
       speex_assert(0);                                \
     }                                                 \
   } while(0)

#endif
