/* Copyright (C) 2002 Jean-Marc Valin 
   File: speex_bits.h

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

typedef struct SpeexBits {
   char bytes[MAX_BYTES_PER_FRAME];
   int  nbBits;
   int  bytePtr;
   int  bitPtr;
} SpeexBits;

void speex_bits_init(SpeexBits *bits);

void speex_bits_destroy(SpeexBits *bits);

void speex_bits_reset(SpeexBits *bits);

void speex_bits_rewind(SpeexBits *bits);

void speex_bits_init_from(SpeexBits *bits, char *bytes, int len);

void speex_bits_read_whole_bytes(SpeexBits *bits, char *bytes, int len);

int speex_bits_write(SpeexBits *bits, char *bytes, int max_len);

int speex_bits_write_whole_bytes(SpeexBits *bits, char *bytes, int max_len);

void speex_bits_pack(SpeexBits *bits, int data, int nbBits);

int speex_bits_unpack_signed(SpeexBits *bits, int nbBits);

unsigned int speex_bits_unpack_unsigned(SpeexBits *bits, int nbBits);


#endif
