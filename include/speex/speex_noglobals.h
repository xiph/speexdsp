/* Copyright (C) 2004 CSIRO Australia */
/* Copyright (C) 2002 Jean-Marc Valin*/
/**
  @file speex_noglobals.h
  @brief Dynamically allocates the different modes of the codec
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

#ifndef SPEEX_NOGLOBALS_H
#define SPEEX_NOGLOBALS_H

#ifdef __cplusplus
extern "C" {
#endif

/** Default narrowband mode */
const SpeexMode * speex_nb_mode_new (void);
void speex_nb_mode_free (const SpeexMode * mode);

/** Default wideband mode */
const SpeexMode * speex_wb_mode_new (void);
void speex_wb_mode_free (const SpeexMode * mode);

/** Default "ultra-wideband" mode */
const SpeexMode * speex_uwb_mode_new (void);
void speex_uwb_mode_free (const SpeexMode * mode);

/** Query modes available */
const SpeexMode * speex_mode_new_byID (int id);

/** Free a mode, using its ID */
void speex_mode_free_byID (SpeexMode * mode, int id);

#ifdef __cplusplus
}
#endif


#endif
