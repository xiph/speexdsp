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

/* TODO:
 - Implement JITTER_BUFFER_SET_DESTROY_CALLBACK
 - Multiple get / get_multiple()
 - User data?
 - Use JitterBufferPacket internally
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include "misc.h"
#include <speex/speex.h>
#include <speex/speex_bits.h>
#include <speex/speex_jitter.h>
#include "os_support.h"

#ifndef NULL
#define NULL 0
#endif

#define LATE_BINS 15
#define MAX_MARGIN 30                     /**< Number of bins in margin histogram */

#define SPEEX_JITTER_MAX_BUFFER_SIZE 200   /**< Maximum number of packets in jitter buffer */



#define GT32(a,b) (((spx_int32_t)((a)-(b)))>0)
#define GE32(a,b) (((spx_int32_t)((a)-(b)))>=0)
#define LT32(a,b) (((spx_int32_t)((a)-(b)))<0)
#define LE32(a,b) (((spx_int32_t)((a)-(b)))<=0)

/** Jitter buffer structure */
struct JitterBuffer_ {
   spx_uint32_t pointer_timestamp;                                        /**< Timestamp of what we will *get* next */
   spx_uint32_t current_timestamp;                                        /**< Timestamp of the local clock (what we will *play* next) */

   char *buf[SPEEX_JITTER_MAX_BUFFER_SIZE];                               /**< Buffer of packets (NULL if slot is free) */
   spx_uint32_t timestamp[SPEEX_JITTER_MAX_BUFFER_SIZE];                  /**< Timestamp of packet                 */
   int span[SPEEX_JITTER_MAX_BUFFER_SIZE];                                /**< Timestamp of packet                 */
   int len[SPEEX_JITTER_MAX_BUFFER_SIZE];                                 /**< Number of bytes in packet           */
   
   void (*destroy) (void *);                                           /**< Callback for destroying a packet */

   spx_int32_t sub_clock;                                                 /** Time since last get() call (incremented by tick()) */
   int resolution;                                                        /**< Time resolution for histogram (timestamp units) */
   int delay_step;                                                        /**< Size of the steps when adjusting buffering (timestamp units) */
   int reset_state;                                                       /**< True if state was just reset        */
   int buffer_margin;                                                     /**< How many frames we want to keep in the buffer (lower bound) */
   int late_cutoff;                                                       /**< How late must a packet be for it not to be considered at all */
   int interp_requested;                                                  /**< An interpolation is requested by speex_jitter_update_delay() */

   int lost_count;                                                        /**< Number of consecutive lost packets  */
   float shortterm_margin[MAX_MARGIN];                                    /**< Short term margin histogram         */
   float longterm_margin[MAX_MARGIN];                                     /**< Long term margin histogram          */
   float loss_rate;                                                       /**< Average loss rate                   */
};

/** Initialise jitter buffer */
JitterBuffer *jitter_buffer_init(int resolution)
{
   JitterBuffer *jitter = (JitterBuffer*)speex_alloc(sizeof(JitterBuffer));
   if (jitter)
   {
      int i;
      for (i=0;i<SPEEX_JITTER_MAX_BUFFER_SIZE;i++)
         jitter->buf[i]=NULL;
      jitter->resolution = resolution;
      jitter->delay_step = resolution;
      jitter->buffer_margin = 1;
      jitter->late_cutoff = 50;
      jitter->destroy = NULL;
      jitter_buffer_reset(jitter);
   }
   return jitter;
}

/** Reset jitter buffer */
void jitter_buffer_reset(JitterBuffer *jitter)
{
   int i;
   for (i=0;i<SPEEX_JITTER_MAX_BUFFER_SIZE;i++)
   {
      if (jitter->buf[i])
      {
         if (jitter->destroy)
            jitter->destroy(jitter->buf[i]);
         else
            speex_free(jitter->buf[i]);
         jitter->buf[i] = NULL;
      }
   }
   /* Timestamp is actually undefined at this point */
   jitter->pointer_timestamp = 0;
   jitter->current_timestamp = 0;
   jitter->reset_state = 1;
   jitter->lost_count = 0;
   jitter->loss_rate = 0;
   for (i=0;i<MAX_MARGIN;i++)
   {
      jitter->shortterm_margin[i] = 0;
      jitter->longterm_margin[i] = 0;
   }
   /*fprintf (stderr, "reset\n");*/
}

/** Destroy jitter buffer */
void jitter_buffer_destroy(JitterBuffer *jitter)
{
   jitter_buffer_reset(jitter);
   speex_free(jitter);
}

/** Put one packet into the jitter buffer */
void jitter_buffer_put(JitterBuffer *jitter, const JitterBufferPacket *packet)
{
   int i,j;
   spx_int32_t arrival_margin;
   spx_int32_t arrival_time;
   /*fprintf (stderr, "put packet %d %d\n", timestamp, span);*/
   if (jitter->reset_state)
   {
      jitter->reset_state=0;
      jitter->pointer_timestamp = packet->timestamp;
      jitter->current_timestamp = packet->timestamp;
      /*fprintf(stderr, "reset to %d\n", timestamp);*/
   }
   
   /* Cleanup buffer (remove old packets that weren't played) */
   for (i=0;i<SPEEX_JITTER_MAX_BUFFER_SIZE;i++)
   {
      /* Make sure we don't discard a "just-late" packet in case we want to play it next (if we interpolate). */
      if (jitter->buf[i] && LE32(jitter->timestamp[i] + jitter->span[i], jitter->pointer_timestamp))
      {
         /*fprintf (stderr, "cleaned (not played)\n");*/
         if (jitter->destroy)
            jitter->destroy(jitter->buf[i]);
         else
            speex_free(jitter->buf[i]);
         jitter->buf[i] = NULL;
      }
   }

   /*Find an empty slot in the buffer*/
   for (i=0;i<SPEEX_JITTER_MAX_BUFFER_SIZE;i++)
   {
      if (jitter->buf[i]==NULL)
         break;
   }

   /*fprintf(stderr, "%d %d %f\n", timestamp, jitter->pointer_timestamp, jitter->drift_average);*/
   /*No place left in the buffer*/
   if (i==SPEEX_JITTER_MAX_BUFFER_SIZE)
   {
      int earliest=jitter->timestamp[0];
      i=0;
      for (j=1;j<SPEEX_JITTER_MAX_BUFFER_SIZE;j++)
      {
         if (!jitter->buf[i] || LT32(jitter->timestamp[j],earliest))
         {
            earliest = jitter->timestamp[j];
            i=j;
         }
      }
      if (jitter->destroy)
         jitter->destroy(jitter->buf[i]);
      else
         speex_free(jitter->buf[i]);
      jitter->buf[i]=NULL;
      if (jitter->lost_count>20)
      {
         jitter_buffer_reset(jitter);
      }
      /*fprintf (stderr, "Buffer is full, discarding earliest frame %d (currently at %d)\n", timestamp, jitter->pointer_timestamp);*/      
   }
   
   /* Copy packet in buffer */
   if (jitter->destroy)
   {
      jitter->buf[i] = packet->data;
   } else {
      jitter->buf[i]=(char*)speex_alloc(packet->len);
      for (j=0;j<packet->len;j++)
         jitter->buf[i][j]=packet->data[j];
   }
   jitter->timestamp[i]=packet->timestamp;
   jitter->span[i]=packet->span;
   jitter->len[i]=packet->len;
   
   if (jitter->sub_clock == -1)
      arrival_time = jitter->pointer_timestamp;
   else
      arrival_time = jitter->current_timestamp + jitter->sub_clock;
   
   if (jitter->current_timestamp + jitter->sub_clock > jitter->pointer_timestamp)
      speex_warning("something's wrong with the time");
   
   /* Adjust the buffer size depending on network conditions.
      The arrival margin is how much in advance (or late) the packet it */
   /* FIXME: We might be one tick off in considering what's late */
   arrival_margin = (((spx_int32_t)packet->timestamp) - ((spx_int32_t)arrival_time))/jitter->resolution - jitter->buffer_margin;
   
   if (arrival_margin >= -jitter->late_cutoff)
   {
      /* Here we compute the histogram based on the time of arrival of the packet.
         This is based on a (first-order) recursive average. We keep both a short-term
         histogram and a long-term histogram */
      spx_int32_t int_margin;
      /* First, apply the "damping" of the recursive average to all bins */
      for (i=0;i<MAX_MARGIN;i++)
      {
         jitter->shortterm_margin[i] *= .98;
         jitter->longterm_margin[i] *= .995;
      }
      /* What histogram bin the packet should be counted in */
      int_margin = LATE_BINS + arrival_margin;
      if (int_margin>MAX_MARGIN-1)
         int_margin = MAX_MARGIN-1;
      if (int_margin<0)
         int_margin = 0;
      /* Add the packet to the right bin */
      jitter->shortterm_margin[int_margin] += .02;
      jitter->longterm_margin[int_margin] += .005;
   } else {
      /* Packet has arrived *way* too late, we pretty much consider it lost and not take it into account in the histogram */
      /*fprintf (stderr, "way too late = %d\n", arrival_margin);*/
      if (jitter->lost_count>20)
      {
         jitter_buffer_reset(jitter);
      }
   }
#if 0 /* Enable to check how much is being buffered */
   if (rand()%1000==0)
   {
      int count = 0;
      for (j=0;j<SPEEX_JITTER_MAX_BUFFER_SIZE;j++)
      {
         if (jitter->buf[j])
            count++;
      }
      fprintf (stderr, "buffer_size = %d\n", count);
   }
#endif
}

/** Get one packet from the jitter buffer */
int jitter_buffer_get(JitterBuffer *jitter, JitterBufferPacket *packet, spx_int32_t desired_span, spx_int32_t *start_offset)
{
   int i;
   unsigned int j;
   float late_ratio_short;
   float late_ratio_long;
   float ontime_ratio_short;
   float ontime_ratio_long;
   float early_ratio_short;
   float early_ratio_long;
   int incomplete = 0;
   
   jitter->sub_clock = -1;
   jitter->current_timestamp = jitter->pointer_timestamp;
   
   if (jitter->interp_requested)
   {
      jitter->interp_requested = 0;
      if (start_offset)
         *start_offset = 0;
      packet->timestamp = jitter->pointer_timestamp;
      packet->span = jitter->delay_step;
      jitter->pointer_timestamp += jitter->delay_step;
      packet->len = 0;
      return JITTER_BUFFER_MISSING;
   }
   /*fprintf (stderr, "get packet %d %d\n", jitter->pointer_timestamp, jitter->current_timestamp);*/
   
   /* Compiling arrival statistics */
   
   late_ratio_short = 0;
   late_ratio_long = 0;
   /* Count the proportion of packets that are late */
   for (i=0;i<LATE_BINS;i++)
   {
      late_ratio_short += jitter->shortterm_margin[i];
      late_ratio_long += jitter->longterm_margin[i];
   }
   /* Count the proportion of packets that are just on time */
   ontime_ratio_short = jitter->shortterm_margin[LATE_BINS];
   ontime_ratio_long = jitter->longterm_margin[LATE_BINS];
   early_ratio_short = early_ratio_long = 0;
   /* Count the proportion of packets that are early */
   for (i=LATE_BINS+1;i<MAX_MARGIN;i++)
   {
      early_ratio_short += jitter->shortterm_margin[i];
      early_ratio_long += jitter->longterm_margin[i];
   }
   if (0&&jitter->pointer_timestamp%1000==0)
   {
      /*fprintf (stderr, "%f %f %f %f %f %f\n", early_ratio_short, early_ratio_long, ontime_ratio_short, ontime_ratio_long, late_ratio_short, late_ratio_long);*/
      /*fprintf (stderr, "%f %f\n", early_ratio_short + ontime_ratio_short + late_ratio_short, early_ratio_long + ontime_ratio_long + late_ratio_long);*/
   }
   
   
   /* Searching for the packet that fits best */
   
   /* Search the buffer for a packet with the right timestamp and spanning the whole current chunk */
   for (i=0;i<SPEEX_JITTER_MAX_BUFFER_SIZE;i++)
   {
      if (jitter->buf[i] && jitter->timestamp[i]==jitter->pointer_timestamp && GE32(jitter->timestamp[i]+jitter->span[i],jitter->pointer_timestamp+desired_span))
         break;
   }
   
   /* If no match, try for an "older" packet that still spans (fully) the current chunk */
   if (i==SPEEX_JITTER_MAX_BUFFER_SIZE)
   {
      for (i=0;i<SPEEX_JITTER_MAX_BUFFER_SIZE;i++)
      {
         if (jitter->buf[i] && LE32(jitter->timestamp[i], jitter->pointer_timestamp) && GE32(jitter->timestamp[i]+jitter->span[i],jitter->pointer_timestamp+desired_span))
            break;
      }
   }
   
   /* If still no match, try for an "older" packet that spans part of the current chunk */
   if (i==SPEEX_JITTER_MAX_BUFFER_SIZE)
   {
      for (i=0;i<SPEEX_JITTER_MAX_BUFFER_SIZE;i++)
      {
         if (jitter->buf[i] && LE32(jitter->timestamp[i], jitter->pointer_timestamp) && GT32(jitter->timestamp[i]+jitter->span[i],jitter->pointer_timestamp))
            break;
      }
   }
   
   /* If still no match, try for earliest packet possible */
   if (i==SPEEX_JITTER_MAX_BUFFER_SIZE)
   {
      int found = 0;
      spx_uint32_t best_time=0;
      int best_span=0;
      int besti=0;
      for (i=0;i<SPEEX_JITTER_MAX_BUFFER_SIZE;i++)
      {
         /* check if packet starts within current chunk */
         if (jitter->buf[i] && LT32(jitter->timestamp[i],jitter->pointer_timestamp+desired_span) && GE32(jitter->timestamp[i],jitter->pointer_timestamp))
         {
            if (!found || LT32(jitter->timestamp[i],best_time) || (jitter->timestamp[i]==best_time && GT32(jitter->span[i],best_span)))
            {
               best_time = jitter->timestamp[i];
               best_span = jitter->span[i];
               besti = i;
               found = 1;
            }
         }
      }
      if (found)
      {
         i=besti;
         incomplete = 1;
         /*fprintf (stderr, "incomplete: %d %d %d %d\n", jitter->timestamp[i], jitter->pointer_timestamp, chunk_size, jitter->span[i]);*/
      }
   }

   /* If we find something */
   if (i!=SPEEX_JITTER_MAX_BUFFER_SIZE)
   {
      /* We (obviously) haven't lost this packet */
      jitter->lost_count = 0;
      jitter->loss_rate = .999*jitter->loss_rate;
      /* Check for potential overflow */
      packet->len = jitter->len[i];
      /* Copy packet */
      if (jitter->destroy)
      {
         packet->data = jitter->buf[i];
      } else {
         for (j=0;j<packet->len;j++)
            packet->data[j] = jitter->buf[i][j];
         /* Remove packet */
         speex_free(jitter->buf[i]);
      }
      jitter->buf[i] = NULL;
      /* Set timestamp and span (if requested) */
      if (start_offset)
         *start_offset = (spx_int32_t)jitter->timestamp[i]-(spx_int32_t)jitter->pointer_timestamp;
      packet->timestamp = jitter->timestamp[i];
      packet->span = jitter->span[i];
      /* Point to the end of the current packet */
      jitter->pointer_timestamp = jitter->timestamp[i]+jitter->span[i];
      if (incomplete)
         return JITTER_BUFFER_INCOMPLETE;
      else
         return JITTER_BUFFER_OK;
   }
   
   
   /* If we haven't found anything worth returning */
   
   /*fprintf (stderr, "not found\n");*/
   jitter->lost_count++;
   /*fprintf (stderr, "m");*/
   /*fprintf (stderr, "lost_count = %d\n", jitter->lost_count);*/
   jitter->loss_rate = .999*jitter->loss_rate + .001;
   if (start_offset)
      *start_offset = 0;
   packet->timestamp = jitter->pointer_timestamp;
   packet->span = desired_span;
   jitter->pointer_timestamp += desired_span;
   packet->len = 0;
   
   /* Adjusting the buffering bssed on the amount of packets that are early/on time/late */   
   if (late_ratio_short > .1 || late_ratio_long > .03)
   {
      /* If too many packets are arriving late */
      jitter->shortterm_margin[MAX_MARGIN-1] += jitter->shortterm_margin[MAX_MARGIN-2];
      jitter->longterm_margin[MAX_MARGIN-1] += jitter->longterm_margin[MAX_MARGIN-2];
      for (i=MAX_MARGIN-3;i>=0;i--)
      {
         jitter->shortterm_margin[i+1] = jitter->shortterm_margin[i];
         jitter->longterm_margin[i+1] = jitter->longterm_margin[i];         
      }
      jitter->shortterm_margin[0] = 0;
      jitter->longterm_margin[0] = 0;            
      jitter->pointer_timestamp -= jitter->delay_step;
      jitter->current_timestamp -= jitter->delay_step;
      /*fprintf (stderr, "i");*/
      /*fprintf (stderr, "interpolate (getting some slack)\n");*/
   }

   return JITTER_BUFFER_MISSING;

}

/** Get pointer timestamp of jitter buffer */
int jitter_buffer_get_pointer_timestamp(JitterBuffer *jitter)
{
   return jitter->pointer_timestamp;
}

void jitter_buffer_tick(JitterBuffer *jitter)
{
   if (jitter->sub_clock == -1)
      jitter->sub_clock = 0;
   else
      jitter->sub_clock += jitter->resolution;
}

/* Let the jitter buffer know it's the right time to adjust the buffering delay to the network conditions */
int jitter_buffer_update_delay(JitterBuffer *jitter, JitterBufferPacket *packet, spx_int32_t *start_offset)
{
   int i;
   float late_ratio_short;
   float late_ratio_long;
   float ontime_ratio_short;
   float ontime_ratio_long;
   float early_ratio_short;
   float early_ratio_long;
   
   /*fprintf (stderr, "get packet %d %d\n", jitter->pointer_timestamp, jitter->current_timestamp);*/

   /* FIXME: This should be only what remaining of the current tick */
   late_ratio_short = 0;
   late_ratio_long = 0;
   /* Count the proportion of packets that are late */
   for (i=0;i<LATE_BINS;i++)
   {
      late_ratio_short += jitter->shortterm_margin[i];
      late_ratio_long += jitter->longterm_margin[i];
   }
   /* Count the proportion of packets that are just on time */
   ontime_ratio_short = jitter->shortterm_margin[LATE_BINS];
   ontime_ratio_long = jitter->longterm_margin[LATE_BINS];
   early_ratio_short = early_ratio_long = 0;
   /* Count the proportion of packets that are early */
   for (i=LATE_BINS+1;i<MAX_MARGIN;i++)
   {
      early_ratio_short += jitter->shortterm_margin[i];
      early_ratio_long += jitter->longterm_margin[i];
   }
   
   /* Adjusting the buffering bssed on the amount of packets that are early/on time/late */   
   if (late_ratio_short > .1 || late_ratio_long > .03)
   {
      /* If too many packets are arriving late */
      jitter->shortterm_margin[MAX_MARGIN-1] += jitter->shortterm_margin[MAX_MARGIN-2];
      jitter->longterm_margin[MAX_MARGIN-1] += jitter->longterm_margin[MAX_MARGIN-2];
      for (i=MAX_MARGIN-3;i>=0;i--)
      {
         jitter->shortterm_margin[i+1] = jitter->shortterm_margin[i];
         jitter->longterm_margin[i+1] = jitter->longterm_margin[i];         
      }
      jitter->shortterm_margin[0] = 0;
      jitter->longterm_margin[0] = 0;            
      jitter->pointer_timestamp -= jitter->delay_step;
      jitter->current_timestamp -= jitter->delay_step;
      jitter->interp_requested = 1;
      return JITTER_BUFFER_ADJUST_INTERPOLATE;
   
   } else if (late_ratio_short + ontime_ratio_short < .005 && late_ratio_long + ontime_ratio_long < .01 && early_ratio_short > .8)
   {
      /* Many frames arriving early */
      jitter->shortterm_margin[0] += jitter->shortterm_margin[1];
      jitter->longterm_margin[0] += jitter->longterm_margin[1];
      for (i=1;i<MAX_MARGIN-1;i++)
      {
         jitter->shortterm_margin[i] = jitter->shortterm_margin[i+1];
         jitter->longterm_margin[i] = jitter->longterm_margin[i+1];         
      }
      jitter->shortterm_margin[MAX_MARGIN-1] = 0;
      jitter->longterm_margin[MAX_MARGIN-1] = 0;      
      /*fprintf (stderr, "drop frame\n");*/
      /*fprintf (stderr, "d");*/
      jitter->pointer_timestamp += jitter->delay_step;
      jitter->current_timestamp += jitter->delay_step;
      return JITTER_BUFFER_ADJUST_DROP;
   }
   
   return JITTER_BUFFER_ADJUST_OK;
}

/* Used like the ioctl function to control the jitter buffer parameters */
int jitter_buffer_ctl(JitterBuffer *jitter, int request, void *ptr)
{
   int count, i;
   switch(request)
   {
      case JITTER_BUFFER_SET_MARGIN:
         jitter->buffer_margin = *(spx_int32_t*)ptr;
         break;
      case JITTER_BUFFER_GET_MARGIN:
         *(spx_int32_t*)ptr = jitter->buffer_margin;
         break;
      case JITTER_BUFFER_GET_AVALIABLE_COUNT:
         count = 0;
         for (i=0;i<SPEEX_JITTER_MAX_BUFFER_SIZE;i++)
         {
            if (jitter->buf[i] && LE32(jitter->pointer_timestamp, jitter->timestamp[i]))
            {
               count++;
            }
         }
         *(spx_int32_t*)ptr = count;
         break;
      case JITTER_BUFFER_SET_DESTROY_CALLBACK:
         jitter->destroy = (void (*) (void *))ptr;
         break;
      case JITTER_BUFFER_GET_DESTROY_CALLBACK:
         *(void (**) (void *))ptr = jitter->destroy;
         break;
      case JITTER_BUFFER_SET_DELAY_STEP:
         jitter->delay_step = *(spx_int32_t*)ptr;
         break;
      case JITTER_BUFFER_GET_DELAY_STEP:
         *(spx_int32_t*)ptr = jitter->delay_step;
         break;
      default:
         speex_warning_int("Unknown jitter_buffer_ctl request: ", request);
         return -1;
   }
   return 0;
}

