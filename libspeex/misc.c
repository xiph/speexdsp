/* Copyright (C) 2002 Jean-Marc Valin 
   File: mics.c
   Various utility routines for Speex

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

#include "misc.h"
#include <stdlib.h>
#include <string.h>

unsigned int be_int(unsigned int i)
{
   unsigned int ret=i;
#ifndef WORDS_BIGENDIAN
   ret =  i>>24;
   ret += (i>>8)&0x0000ff00;
   ret += (i<<8)&0x00ff0000;
   ret += (i<<24);
#endif
   return ret;
}

unsigned int le_int(unsigned int i)
{
   unsigned int ret=i;
#ifdef WORDS_BIGENDIAN
   ret =  i>>24;
   ret += (i>>8)&0x0000ff00;
   ret += (i<<8)&0x00ff0000;
   ret += (i<<24);
#endif
   return ret;
}

unsigned short be_short(unsigned short s)
{
   unsigned short ret=s;
#ifndef WORDS_BIGENDIAN
   ret =  s>>8;
   ret += s<<8;
#endif
   return ret;
}

unsigned short le_short(unsigned short s)
{
   unsigned short ret=s;
#ifdef WORDS_BIGENDIAN
   ret =  s>>8;
   ret += s<<8;
#endif
   return ret;
}

void *speex_alloc (int size)
{
   return calloc(size,1);
}

void speex_free (void *ptr)
{
   free(ptr);
}

void *speex_move (void *dest, void *src, int n)
{
   return memmove(dest,src,n);
}
