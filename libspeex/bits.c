/* Copyright (C) 2002 Jean-Marc Valin 
   File: bits.c

   Handles bit packing/unpacking

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.
   
   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

#include "bits.h"
#include <stdio.h>

void frame_bits_init(FrameBits *bits)
{
   int i;
   for (i=0;i<MAX_BYTES_PER_FRAME;i++)
      bits->bytes[i]=0;
   bits->nbBits=0;
   bits->bytePtr=0;
   bits->bitPtr=0;
}

void frame_bits_destroy(FrameBits *bits)
{
   /* Will do something once the allocation is dynamic */
}

void frame_bits_reset(FrameBits *bits)
{
   int i;
   for (i=0;i<MAX_BYTES_PER_FRAME;i++)
      bits->bytes[i]=0;
   bits->nbBits=0;
   bits->bytePtr=0;
   bits->bitPtr=0;
}

void frame_bits_init_from(FrameBits *bits, char *bytes, int len)
{
   int i;
   if (len > MAX_BYTES_PER_FRAME)
   {
      fprintf (stderr, "Trying to init frame with too many bits");
      exit(1);
   }
   for (i=0;i<len;i++)
      bits->bytes[i]=bytes[i];
   bits->nbBits=len<<3;
   bits->bytePtr=0;
   bits->bitPtr=0;
}

void frame_bits_write(FrameBits *bits, char *bytes, int max_len)
{
   int i;
   if (max_len > ((bits->nbBits+7)>>3))
      max_len = ((bits->nbBits+7)>>3);
   for (i=0;i<max_len;i++)
      bytes[i]=bits->bytes[i];
}

void frame_bits_pack(FrameBits *bits, int data, int nbBits)
{
   int i;
   unsigned int d=data;
   while(nbBits)
   {
      int bit;
      bit = (d>>(nbBits-1))&1;
      bits->bytes[bits->bytePtr] |= bit<<(7-bits->bitPtr);
      bits->bitPtr++;
      fprintf(stderr, "%d %d\n", nbBits, bit);
      if (bits->bitPtr==8)
      {
         bits->bitPtr=0;
         bits->bytePtr++;
      }
      bits->nbBits++;
      nbBits--;
   }
}

int frame_bits_unpack_signed(FrameBits *bits, int nbBits)
{
   unsigned int d=frame_bits_unpack_unsigned(bits,nbBits);
   /* If number is negative */
   if (d>>(nbBits-1))
   {
      d |= (-1)<<nbBits;
   }
   return d;
}

unsigned int frame_bits_unpack_unsigned(FrameBits *bits, int nbBits)
{
   int i;
   unsigned int d=0;
   while(nbBits)
   {
      int bit;
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
