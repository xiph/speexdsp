/* Copyright (C) 2002 Jean-Marc Valin
   File speex_callbacks.c
   Callback handling and in-band signalling


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

#include "speex_callbacks.h"

int speex_inband_handler(SpeexBits *bits, SpeexCallback *callback_list, void *state)
{
   int id;
   SpeexCallback *callback;

   id=speex_bits_unpack_unsigned(bits, 4);
   callback = callback_list+id;

   if (callback->func)
   {
      return callback->func(bits, state, callback->data);
   } else
   {
      int adv;
      if (id<4)
         adv = 4;
      else if (id<8)
         adv = 8;
      else if (id<12)
         adv = 16;
      else if (id<14)
         adv = 32;
      else 
         adv = 64;
      speex_bits_advance(bits, adv);
   }
   return 0;
}

int speex_std_mode_request_handler(SpeexBits *bits, void *state, void *data)
{
   return 0;
}

int speex_std_high_mode_request_handler(SpeexBits *bits, void *state, void *data)
{
   return 0;
}

int speex_std_vbr_quality_request_handler(SpeexBits *bits, void *state, void *data)
{
   return 0;
}


int speex_std_char_handler(SpeexBits *bits, void *state, void *data)
{
   return 0;
}



int speex_default_user_handler(SpeexBits *bits, void *state, void *data)
{
   int req_size = speex_bits_unpack_unsigned(bits, 6);
   speex_bits_advance(bits, 5+8*req_size);
   return 0;
}
