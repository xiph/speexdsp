/* Copyright (C) 2002 Jean-Marc Valin 
   File: modes.h

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

#ifndef MODES_H
#define MODES_H

typedef struct SpeexMode {
   int     frameSize;
   int     subframeSize;
   int     windowSize;
   int     lpcSize;
   int     bufSize;
   int     pitchStart;
   int     pitchEnd;
   float   gamma1;
   float   gamma2;
   /* Should add info about LSP quantization, pitch gain quantization 
      and other codebooks */
} SpeexMode;

extern SpeexMode nb_mode;

#endif
