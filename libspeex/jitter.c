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

#ifndef NULL
#define NULL 0
#endif

#include "speex.h"
#include "speex_bits.h"
#include "speex_jitter.h"
#include <stdio.h>

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

   jitter->buffer_size = 4;

   jitter->pointer_timestamp = -jitter->frame_time * jitter->buffer_size;
   jitter->reset_state = 1;
   jitter->lost_count = 0;
   jitter->loss_rate = 0;
}

void speex_jitter_destroy(SpeexJitter *jitter)
{
}


void speex_jitter_put(SpeexJitter *jitter, char *packet, int len, int timestamp)
{
   int i,j;
   int timestamp_offset;
   int arrival_margin;
   
   if (jitter->reset_state)
   {
      jitter->reset_state=0;
      jitter->pointer_timestamp = timestamp-jitter->frame_time * jitter->buffer_size;
      jitter->drift_average = 0;
      jitter->late_ratio = 0;
      jitter->early_ratio = 0;
      jitter->ontime_ratio = 0;
   }
   
   /* Cleanup buffer (remove old packets that weren't played) */
   for (i=0;i<SPEEX_JITTER_MAX_BUFFER_SIZE;i++)
   {
      if (jitter->timestamp[i]<jitter->pointer_timestamp)
      {
         jitter->len[i]=-1;
         /*if (jitter->timestamp[i] != -1)
            fprintf (stderr, "discarding %d %d\n", jitter->timestamp[i], jitter->pointer_timestamp);*/
      }
   }

   /*Find an empty slot in the buffer*/
   for (i=0;i<SPEEX_JITTER_MAX_BUFFER_SIZE;i++)
   {
      if (jitter->len[i]==-1)
         break;
   }

   /*fprintf(stderr, "%d %d %f\n", timestamp, jitter->pointer_timestamp, jitter->drift_average);*/
   if (i==SPEEX_JITTER_MAX_BUFFER_SIZE)
   {
      fprintf (stderr, "Buffer is full, discarding frame\n");
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
   
   /* Don't count late packets when adjusting the synchro (we're taking care of them elsewhere) */
   if (timestamp > jitter->pointer_timestamp)
   {
      timestamp_offset = timestamp-jitter->pointer_timestamp - (jitter->buffer_size+1)*jitter->frame_time;
      jitter->drift_average = .99*jitter->drift_average + .01*timestamp_offset;
   } else {
      fprintf (stderr, "frame for timestamp %d arrived too late\n", timestamp);
   }

   /* Adjust the buffer size depending on network conditions */
   arrival_margin = (timestamp - jitter->pointer_timestamp - jitter->frame_time) - (int)jitter->drift_average;
   if (arrival_margin >= -2*jitter->frame_time)
   {
      jitter->late_ratio *= .99;
      jitter->early_ratio *= .99;
      jitter->ontime_ratio *= .99;
   }
   
   if (arrival_margin < 0 && arrival_margin >= -3*jitter->frame_time)
   {
      jitter->late_ratio += .01;
   } else if (arrival_margin <= jitter->frame_time && arrival_margin >= 0) 
   {
      jitter->ontime_ratio += .01;
   } else if (arrival_margin > jitter->frame_time)
   {
      jitter->early_ratio += .01;
   }
   
   if (jitter->late_ratio + jitter->ontime_ratio < .01 && jitter->early_ratio > .8)
   {
      jitter->buffer_size--;
      jitter->early_ratio = 0;
      jitter->drift_average += jitter->frame_time;
      fprintf (stderr, "buffer %d -> %d\n", jitter->buffer_size+1, jitter->buffer_size);
   }
   if (jitter->late_ratio > .03)
   {
      jitter->buffer_size++;
      jitter->late_ratio = 0;
      jitter->drift_average -= jitter->frame_time;
      fprintf (stderr, "buffer %d -> %d\n", jitter->buffer_size-1, jitter->buffer_size);
   }
   /*fprintf (stderr, "margin : %d %d %f %f %f %f\n", arrival_margin, jitter->buffer_size, 100*jitter->loss_rate, 100*jitter->late_ratio, 100*jitter->ontime_ratio, 100*jitter->early_ratio);*/
}

void speex_jitter_get(SpeexJitter *jitter, short *out)
{
   int i;
   int ret;
   int drop, interp;
   int drop_margin, interp_margin;
   
   drop_margin = (1+jitter->buffer_size) * jitter->frame_time/2;
   interp_margin = jitter->frame_time/2;
   /*drop_margin = 15;
   interp_margin = 15;*/
   
   /* Handle frame interpolation (playing too fast) */
   if (-jitter->drift_average > interp_margin)
   {
      fprintf (stderr, "interpolate frame\n");
      speex_decode_int(jitter->dec, NULL, out);
      jitter->drift_average += jitter->frame_time;
      return;
   }

   /* Increment timestamp */
   jitter->pointer_timestamp += jitter->frame_time;

   /* Handle frame dropping (receiving too fast) */
   if (jitter->drift_average > drop_margin)
   {
      fprintf (stderr, "drop frame\n");
      jitter->pointer_timestamp += jitter->frame_time;
      jitter->drift_average -= jitter->frame_time;
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
         ret = speex_decode_int(jitter->dec, &jitter->current_packet, out);
         if (ret == 0)
            return;
         else
            jitter->valid_bits = 0;
      }

      fprintf (stderr, "lost/late frame %d\n", jitter->pointer_timestamp);
      /*Packet is late or lost*/
      speex_decode_int(jitter->dec, NULL, out);
      jitter->lost_count++;
      if (jitter->lost_count>=25)
      {
         jitter->lost_count = 0;
         jitter->reset_state = 1;
         speex_decoder_ctl(jitter->dec, SPEEX_RESET_STATE, NULL);
      }
      jitter->loss_rate = .999*jitter->loss_rate + .001;
   } else {
      /* Found the right packet */
      speex_bits_read_from(&jitter->current_packet, jitter->buf[i], jitter->len[i]);
      jitter->len[i]=-1;
      /* Decode packet */
      ret = speex_decode_int(jitter->dec, &jitter->current_packet, out);
      if (ret == 0)
      {
         jitter->valid_bits = 1;
      } else {
         /* Error while decoding */
         for (i=0;i<jitter->frame_size;i++)
            out[i]=0;
      }
      jitter->loss_rate = .999*jitter->loss_rate;
   }


}

int speex_jitter_get_pointer_timestamp(SpeexJitter *jitter)
{
   return jitter->pointer_timestamp;
}

int speex_jitter_set_buffersize(SpeexJitter *jitter, int size)
{
   jitter->buffer_size = size;
}
