/* Copyright (C) 2006 Jean-Marc Valin */
/**
   @file filterbank.c
   @brief Converting between psd and filterbank
 */
/*
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:

   1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

   3. The name of the author may not be used to endorse or promote products
   derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
   IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
   OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
   INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
   HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
   STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "filterbank.h"
#include "misc.h"
#include <math.h>

#define toBARK(n)   (13.1f*atan(.00074f*(n))+2.24f*atan((n)*(n)*1.85e-8f)+1e-4f*(n))
#define toMEL(n)    (2595.f*log10(1.f+(n)/700.f))
      
#ifdef FIXED_POINT
#define Q15(x) (floor(.5+32767.*(x)))
#else
#define Q15(x) (x)
#endif
      
FilterBank *filterbank_new(int banks, float sampling, int len, int type)
{
   FilterBank *bank;
   float df, max_mel, mel_interval;
   int i;
   int id1;
   int id2;
   df = .5*sampling/len;
   max_mel = toBARK(.5*sampling);
   mel_interval = max_mel/(banks-1);
   
   bank = (FilterBank*)speex_alloc(sizeof(FilterBank));
   bank->nb_banks = banks;
   bank->len = len;
   bank->bank_left = (int*)speex_alloc(len*sizeof(int));
   bank->bank_right = (int*)speex_alloc(len*sizeof(int));
   bank->filter_left = (spx_word16_t*)speex_alloc(len*sizeof(spx_word16_t));
   bank->filter_right = (spx_word16_t*)speex_alloc(len*sizeof(spx_word16_t));
   /* Think I can safely disable normalisation that for fixed-point (and probably float as well) */
#ifndef FIXED_POINT
   bank->scaling = (float*)speex_alloc(banks*sizeof(float));
#endif
   for (i=0;i<len;i++)
   {
      float curr_freq;
      float mel;
      float val;
      curr_freq = i*df;
      mel = toBARK(curr_freq);
      if (mel > max_mel)
         break;
      id1 = (int)(floor(mel/mel_interval));
      if (id1>banks-2)
      {
         id1 = banks-2;
         val = 1;
      } else {
         val = (mel - id1*mel_interval)/mel_interval;
      }
      id2 = id1+1;
      bank->bank_left[i] = id1;
      bank->filter_left[i] = Q15(1-val);
      bank->bank_right[i] = id2;
      bank->filter_right[i] = Q15(val);
   }
   
   /* Think I can safely disable normalisation that for fixed-point (and probably float as well) */
#ifndef FIXED_POINT
   for (i=0;i<bank->nb_banks;i++)
      bank->scaling[i] = 0;
   for (i=0;i<bank->len;i++)
   {
      int id = bank->bank_left[i];
      bank->scaling[id] += bank->filter_left[i];
      id = bank->bank_right[i];
      bank->scaling[id] += bank->filter_right[i];
   }
   for (i=0;i<bank->nb_banks;i++)
      bank->scaling[i] = Q15_ONE/(bank->scaling[i]);
#endif
   return bank;
}

void filterbank_destroy(FilterBank *bank)
{
   speex_free(bank->bank_left);
   speex_free(bank->bank_right);
   speex_free(bank->filter_left);
   speex_free(bank->filter_right);
#ifndef FIXED_POINT
   speex_free(bank->scaling);
#endif
   speex_free(bank);
}

void filterbank_compute_bank32(FilterBank *bank, spx_word32_t *ps, spx_word32_t *mel)
{
   int i;
   for (i=0;i<bank->nb_banks;i++)
      mel[i] = 0;

   for (i=0;i<bank->len;i++)
   {
      int id;
      id = bank->bank_left[i];
      mel[id] += MULT16_32_P15(bank->filter_left[i],ps[i]);
      id = bank->bank_right[i];
      mel[id] += MULT16_32_P15(bank->filter_right[i],ps[i]);
   }
   /* Think I can safely disable normalisation that for fixed-point (and probably float as well) */
#ifndef FIXED_POINT
   /*for (i=0;i<bank->nb_banks;i++)
      mel[i] = MULT16_32_P15(Q15(bank->scaling[i]),mel[i]);
   */
#endif
}

void filterbank_compute_bank(FilterBank *bank, float *ps, float *mel)
{
   int i;
   for (i=0;i<bank->nb_banks;i++)
      mel[i] = 0;

   for (i=0;i<bank->len;i++)
   {
      int id = bank->bank_left[i];
      mel[id] += bank->filter_left[i]*ps[i];
      id = bank->bank_right[i];
      mel[id] += bank->filter_right[i]*ps[i];
   }
#ifndef FIXED_POINT
   for (i=0;i<bank->nb_banks;i++)
      mel[i] *= bank->scaling[i];
#endif
}

void filterbank_compute_psd(FilterBank *bank, float *mel, float *ps)
{
   int i;
   for (i=0;i<bank->len;i++)
   {
      int id = bank->bank_left[i];
      ps[i] = mel[id]*bank->filter_left[i];
      id = bank->bank_right[i];
      ps[i] += mel[id]*bank->filter_right[i];
   }
}

void filterbank_compute_psd16(FilterBank *bank, spx_word16_t *mel, spx_word16_t *ps)
{
   int i;
   for (i=0;i<bank->len;i++)
   {
      spx_word32_t tmp;
      int id1, id2;
      id1 = bank->bank_left[i];
      id2 = bank->bank_right[i];
      tmp = MULT16_16(mel[id1],bank->filter_left[i]);
      tmp += MULT16_16(mel[id2],bank->filter_right[i]);
      ps[i] = EXTRACT16(PSHR32(tmp,15));
   }
}


