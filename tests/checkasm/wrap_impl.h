/* Copyright (C) 2026 Tristan Matthews */
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

/* #includes a libspeexdsp module privately into a wrap_<mod>_<variant>.c TU.
 * The module's wrap_<mod>_rename.h sets CKA_MODULE_SRC to its file name
 * (found through the libspeexdsp include path). Include this AFTER that
 * rename header + wrap.h and any per-variant #undef USE_* / *_FORCE_ON.
 *
 * #undef HAVE_CONFIG_H stops the module (or os_support.h) from re-including
 * config.h, which would re-define the USE_* macros the caller just cleared;
 * HAVE_CONFIG_H is set on the command line, but the #undef still takes
 * effect for the rest of this TU. The full module pulls in its whole
 * public API plus internal helpers, many of which a given TU never calls,
 * so the pragma quiets the resulting -Wunused-function noise. */

#ifndef CKA_MODULE_SRC
#  error "CKA_MODULE_SRC not set; include the module's wrap_*_rename.h first"
#endif

#undef HAVE_CONFIG_H

#if defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunused-function"
#endif

#include CKA_MODULE_SRC

#if defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif
