/* Copyright (C) 2002 Jean-Marc Valin 
   File: speex_header.h

   Describes the different modes of the codec

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


#ifndef SPEEX_HEADER_H
#define SPEEX_HEADER_H

typedef struct SpeexHeader {
   char speex_version[20];
   int speex_header_version;
   int rate;
   int mode;
   int nb_channels;
   int byte_rate;
   int frame_size;
} SpeexHeader;

char *speex_header_to_packet(SpeexHeader *header, int *size);

SpeexHeader *speex_packet_to_header(char *packet, int size);

#endif
