/* Copyright (C) 2002 Jean-Marc Valin 
   File: speex_bits.c

   Handles bit packing/unpacking

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

#include "speex_bits.h"
#include "misc.h"

void speex_bits_init(SpeexBits *bits)
{
   int i;
   bits->bytes = (char*)speex_alloc(MAX_BYTES_PER_FRAME);

   for (i=0;i<MAX_BYTES_PER_FRAME;i++)
      bits->bytes[i]=0;
   bits->nbBits=0;
   bits->bytePtr=0;
   bits->bitPtr=0;
   bits->owner=1;
}

void speex_bits_init_buffer(SpeexBits *bits, void *buff)
{
   int i;
   bits->bytes = (char*)buff;

   for (i=0;i<MAX_BYTES_PER_FRAME;i++)
      bits->bytes[i]=0;
   bits->nbBits=0;
   bits->bytePtr=0;
   bits->bitPtr=0;
   bits->owner=0;
}

void speex_bits_destroy(SpeexBits *bits)
{
   if (bits->owner)
      speex_free(bits->bytes);
   /* Will do something once the allocation is dynamic */
}

void speex_bits_reset(SpeexBits *bits)
{
   int i;
   for (i=0;i<MAX_BYTES_PER_FRAME;i++)
      bits->bytes[i]=0;
   bits->nbBits=0;
   bits->bytePtr=0;
   bits->bitPtr=0;
}

void speex_bits_rewind(SpeexBits *bits)
{
   bits->bytePtr=0;
   bits->bitPtr=0;
}

void speex_bits_read_from(SpeexBits *bits, char *bytes, int len)
{
   int i;
   if (len > MAX_BYTES_PER_FRAME)
   {
      speex_error ("Trying to init frame with too many bits");
   }
   for (i=0;i<len;i++)
      bits->bytes[i]=bytes[i];
   bits->nbBits=len<<3;
   bits->bytePtr=0;
   bits->bitPtr=0;
}

void speex_bits_flush(SpeexBits *bits)
{
   int i;
   if (bits->bytePtr>0)
   {
      for (i=bits->bytePtr;i<((bits->nbBits+7)>>3);i++)
         bits->bytes[i-bits->bytePtr]=bits->bytes[i];
   }
   bits->nbBits -= bits->bytePtr<<3;
   bits->bytePtr=0;
}

void speex_bits_read_whole_bytes(SpeexBits *bits, char *bytes, int len)
{
   int i,pos;
   speex_bits_flush(bits);
   pos=bits->nbBits>>3;
   for (i=0;i<len;i++)
      bits->bytes[pos+i]=bytes[i];
   bits->nbBits+=len<<3;
}

int speex_bits_write(SpeexBits *bits, char *bytes, int max_len)
{
   int i;
   if (max_len > ((bits->nbBits+7)>>3))
      max_len = ((bits->nbBits+7)>>3);
   for (i=0;i<max_len;i++)
      bytes[i]=bits->bytes[i];
   return max_len;
}

int speex_bits_write_whole_bytes(SpeexBits *bits, char *bytes, int max_len)
{
   int i;
   if (max_len > ((bits->nbBits)>>3))
      max_len = ((bits->nbBits)>>3);
   for (i=0;i<max_len;i++)
      bytes[i]=bits->bytes[i];
   
   if (bits->bitPtr>0)
      bits->bytes[0]=bits->bytes[max_len];
   else
      bits->bytes[0]=0;
   for (i=1;i<((bits->nbBits)>>3)+1;i++)
      bits->bytes[i]=0;
   bits->bytePtr=0;
   bits->nbBits &= 7;
   return max_len;
}


void speex_bits_pack(SpeexBits *bits, int data, int nbBits)
{
   unsigned int d=data;
   while(nbBits)
   {
      int bit;
      bit = (d>>(nbBits-1))&1;
      bits->bytes[bits->bytePtr] |= bit<<(7-bits->bitPtr);
      bits->bitPtr++;

      if (bits->bitPtr==8)
      {
         bits->bitPtr=0;
         bits->bytePtr++;
      }
      bits->nbBits++;
      nbBits--;
   }
}

int speex_bits_unpack_signed(SpeexBits *bits, int nbBits)
{
   unsigned int d=speex_bits_unpack_unsigned(bits,nbBits);
   /* If number is negative */
   if (d>>(nbBits-1))
   {
      d |= (-1)<<nbBits;
   }
   return d;
}

unsigned int speex_bits_unpack_unsigned(SpeexBits *bits, int nbBits)
{
   unsigned int d=0;
   while(nbBits)
   {
      d<<=1;
      d |= (bits->bytes[bits->bytePtr]>>(7-bits->bitPtr))&1;
      bits->bitPtr++;
      if (bits->bitPtr==8)
      {
         bits->bitPtr=0;
         bits->bytePtr++;
      }
      nbBits--;
   }
   return d;
}

unsigned int speex_bits_peek_unsigned(SpeexBits *bits, int nbBits)
{
   unsigned int d=0;
   int bitPtr, bytePtr;
   char *bytes;
   bitPtr=bits->bitPtr;
   bytePtr=bits->bytePtr;
   bytes = bits->bytes;
   while(nbBits)
   {
      d<<=1;
      d |= (bytes[bytePtr]>>(7-bitPtr))&1;
      bitPtr++;
      if (bitPtr==8)
      {
         bitPtr=0;
         bytePtr++;
      }
      nbBits--;
   }
   return d;
}

int speex_bits_peek(SpeexBits *bits)
{
   return (bits->bytes[bits->bytePtr]>>(7-bits->bitPtr))&1;
}

void speex_bits_advance(SpeexBits *bits, int n)
{
   int nbytes, nbits;
   nbytes = n >> 3;
   nbits = n & 7;
   
   bits->bytePtr += nbytes;
   bits->bitPtr += nbits;
   
   if (bits->bitPtr>7)
   {
      bits->bitPtr-=8;
      bits->bytePtr++;
   }
}

int speex_bits_nbytes(SpeexBits *bits)
{
   return ((bits->nbBits+7)>>3);
}
