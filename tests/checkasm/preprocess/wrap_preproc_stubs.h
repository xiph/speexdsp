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

/* Link stubs for the fftwrap/filterbank/mdf symbols preprocess.c
 * references (renamed per TU by wrap_preproc_rename.h). The kernel tests
 * never run the preprocessor itself, so they are never called. Include
 * right after ../wrap_impl.h. */

void *spx_fft_init(int size) { (void) size; return 0; }
void spx_fft_destroy(void *table) { (void) table; }
void spx_fft(void *table, spx_word16_t *in, spx_word16_t *out)
{ (void) table; (void) in; (void) out; }
void spx_ifft(void *table, spx_word16_t *in, spx_word16_t *out)
{ (void) table; (void) in; (void) out; }
FilterBank *filterbank_new(int banks, spx_word32_t sampling, int len, int type)
{ (void) banks; (void) sampling; (void) len; (void) type; return 0; }
void filterbank_destroy(FilterBank *bank) { (void) bank; }
void filterbank_compute_bank32(FilterBank *bank, spx_word32_t *ps, spx_word32_t *mel)
{ (void) bank; (void) ps; (void) mel; }
void filterbank_compute_psd16(FilterBank *bank, spx_word16_t *mel, spx_word16_t *psd)
{ (void) bank; (void) mel; (void) psd; }
void speex_echo_get_residual(SpeexEchoState *st, spx_word32_t *Yout, int len)
{ (void) st; (void) Yout; (void) len; }
