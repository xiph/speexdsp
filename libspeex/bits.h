/* Copyright (C) 2002 Jean-Marc Valin 
   File: bits.h

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

#ifndef BITS_H
#define BITS_H

#define MAX_BYTES_PER_FRAME 1000

typedef struct FrameBits {
   char bytes[MAX_BYTES_PER_FRAME];
   int  nbBits;
   int  bytePtr;
   int  bitPtr;
} FrameBits;

void frame_bits_init(FrameBits *bits);

void frame_bits_destroy(FrameBits *bits);

void frame_bits_reset(FrameBits *bits);

void frame_bits_init_from(FrameBits *bits, char *bytes, int len);

void frame_bits_write(FrameBits *bits, char *bytes, int max_len);

void frame_bits_pack(FrameBits *bits, int data, int nbBits);

int frame_bits_unpack_signed(FrameBits *bits, int nbBits);

unsigned int frame_bits_unpack_unsigned(FrameBits *bits, int nbBits);


#endif
