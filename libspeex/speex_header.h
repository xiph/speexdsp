/* Copyright (C) 2002 Jean-Marc Valin */
/**
   @file speex_header.h
   @brief Describes the Speex header
*/
/*
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

#ifdef __cplusplus
extern "C" {
#endif

struct SpeexMode;

/** Maximum number of characters for encoding the Speex version number in the header */
#define SPEEX_HEADER_VERSION_LENGTH 20

/** Current version of the Speex header */
#define SPEEX_HEADER_VERSION -1

/** Speex header info for file-based formats */
typedef struct SpeexHeader {
   char speex_string[8]; /**< Identifies a Speex bit-stream, always set to "Speex   " */
   char speex_version[SPEEX_HEADER_VERSION_LENGTH]; /**< Speex version */
   int speex_header_version; /**< Version number for the header */
   int header_size; /**< Total size of the header ( sizeof(SpeexHeader) ) */
   int rate; /**< Sampling rate used */
   int mode; /**< Mode used (0 for narrowband, 1 for wideband) */
   int mode_bitstream_version; /**< Version ID of the bit-stream */
   int nb_channels; /**< Number of channels encoded */
   int bitrate; /**< Bit-rate used */
   int frame_size; /**< Size of frames */
   int vbr; /**< 1 for a VBR encoding, 0 otherwise */
   int frames_per_packet; /**< Number of frames stored per Ogg packet */
   int reserved1; /**< Reserved for future use */
   int reserved2; /**< Reserved for future use */
   int reserved3; /**< Reserved for future use */
} SpeexHeader;

/** Initializes a SpeexHeader using basic information */
void speex_init_header(SpeexHeader *header, int rate, int nb_channels, struct SpeexMode *m);

/** Creates the header packet from the header itself (mostly involves endianness conversion) */
char *speex_header_to_packet(SpeexHeader *header, int *size);

/** Creates a SpeexHeader from a packet */
SpeexHeader *speex_packet_to_header(char *packet, int size);

#ifdef __cplusplus
}
#endif


#endif
