/* Copyright (C) 2002 Jean-Marc Valin 
   File: speex_jitter.h

   Adaptive jitter buffer for Speex

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

#include <stdlib.h>
#include "speex.h"
#include "speex_bits.h"
#include "speex_jitter.h"

void speex_jitter_init(SpeexJitter *jitter, void *decoder, int sampling_rate)
{
   int i;
   for (i=0;i<SPEEX_JITTER_MAX_BUFFER_SIZE;i++)
   {
      jitter->len[i]=-1;
      jitter->timestamp[i]=-1;
   }

   jitter->dec = decoder;
   speex_decoder_ctl(decoder, SPEEX_GET_FRAME_SIZE, &jitter->frame_size);
   jitter->frame_time = 1000*jitter->frame_size / sampling_rate;

   speex_bits_init(&jitter->current_packet);
   jitter->valid_bits = 0;

   jitter->buffer_size = 10;

   jitter->pointer_timestamp = -jitter->frame_time * jitter->buffer_size;
   
   jitter->underflow_count = 0;
   jitter->drop_frame = 0;
   jitter->interp_frame = 0;
}

void speex_jitter_destroy(SpeexJitter *jitter)
{
}


void speex_jitter_put(SpeexJitter *jitter, char *packet, int len, int timestamp)
{
   int i,j;

   /* Cleanup buffer (remove old packets that weren't played) */
   for (i=0;i<SPEEX_JITTER_MAX_BUFFER_SIZE;i++)
   {
      if (jitter->timestamp[i]<jitter->pointer_timestamp)
         jitter->len[i]=-1;
   }

   /*Find an empty slot in the buffer*/
   for (i=0;i<SPEEX_JITTER_MAX_BUFFER_SIZE;i++)
   {
      if (jitter->len[i]==-1)
         break;
   }

   if (i==SPEEX_JITTER_MAX_BUFFER_SIZE)
   {
      /*No place left in the buffer*/
      
      /*skip some frame(s) */
      return;
   }
   
   /* Copy packet in buffer */
   if (len>SPEEX_JITTER_MAX_PACKET_SIZE)
      len=SPEEX_JITTER_MAX_PACKET_SIZE;
   for (j=0;j<len;j++)
      jitter->buf[i][j]=packet[j];
   jitter->timestamp[i]=timestamp;
   jitter->len[i]=len;
   
   /* Detect when it's time to drop frames (too much stuff in buffer) */
   if (timestamp-jitter->pointer_timestamp > (jitter->buffer_size+3)*jitter->frame_time)
      jitter->drop_frame = 1;

   /*Detect when it's time to interpolate a frame (not wnough stuff in buffer) */
   if (timestamp-jitter->pointer_timestamp < (jitter->buffer_size-3)*jitter->frame_time)
      jitter->underflow_count++;
   else
      jitter->underflow_count=0;
   if (jitter->underflow_count > 10)
   {
      jitter->underflow_count=0;
      jitter->interp_frame = 1;
   }

   /* Adjust the buffer size depending on network conditions */
}

void speex_jitter_get(SpeexJitter *jitter, float *out)
{
   int i;
   int ret;
   
   /* Handle frame interpolation (receiving too fast) */
   if (jitter->interp_frame)
   {
      speex_decode(jitter->dec, NULL, out);
      jitter->interp_frame = 0;
      return;
   }

   /* Increment timestamp */
   jitter->pointer_timestamp += jitter->frame_time;

   /* Handle frame dropping (receiving too fast) */
   if (jitter->drop_frame)
   {
      jitter->pointer_timestamp += jitter->frame_time;
      jitter->drop_frame = 0;
   }

   /* Send zeros while we fill in the buffer */
   if (jitter->pointer_timestamp<0)
   {
      for (i=0;i<jitter->frame_size;i++)
         out[i]=0;
      return;
   }
   
   /* Search the buffer for a packet with the right timestamp */
   for (i=0;i<SPEEX_JITTER_MAX_BUFFER_SIZE;i++)
   {
      if (jitter->len[i]!=-1 && jitter->timestamp[i]==jitter->pointer_timestamp)
         break;
   }
   
   if (i==SPEEX_JITTER_MAX_BUFFER_SIZE)
   {
      /* No packet found */
      if (jitter->valid_bits)
      {
         /* Try decoding last received packet */
         ret = speex_decode(jitter->dec, &jitter->current_packet, out);
         if (ret == 0)
            return;
         else
            jitter->valid_bits = 0;
      }

      /*Packet is late or lost*/
      speex_decode(jitter->dec, NULL, out);
   } else {
      /* Found the right packet */
      speex_bits_read_from(&jitter->current_packet, jitter->buf[i], jitter->len[i]);
      jitter->len[i]=-1;
      /* Decode packet */
      ret = speex_decode(jitter->dec, &jitter->current_packet, out);
      if (ret == 0)
      {
         jitter->valid_bits = 1;
      } else {
         /* Error while decoding */
         for (i=0;i<jitter->frame_size;i++)
            out[i]=0;
      }
   }


}

