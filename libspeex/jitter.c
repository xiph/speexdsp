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

#define TSUB(a,b) ((spx_int32_t)((a)-(b)))

#define GT32(a,b) (((spx_int32_t)((a)-(b)))>0)
#define GE32(a,b) (((spx_int32_t)((a)-(b)))>=0)
#define LT32(a,b) (((spx_int32_t)((a)-(b)))<0)
#define LE32(a,b) (((spx_int32_t)((a)-(b)))<=0)

/** Jitter buffer structure */
struct JitterBuffer_ {
   spx_uint32_t pointer_timestamp;                                        /**< Timestamp of what we will *get* next */
   spx_uint32_t current_timestamp;                                        /**< Timestamp of the local clock (what we will *play* next) */
   spx_uint32_t last_returned_timestamp;
   spx_uint32_t next_stop;
   
   JitterBufferPacket packets[SPEEX_JITTER_MAX_BUFFER_SIZE];              /**< Packets stored in the buffer */
   spx_uint32_t arrival[SPEEX_JITTER_MAX_BUFFER_SIZE];              /**< Packet arrival time (0 means it was late, even though it's a valid timestamp) */
   
   void (*destroy) (void *);                                              /**< Callback for destroying a packet */

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
         jitter->packets[i].data=NULL;
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
      if (jitter->packets[i].data)
      {
         if (jitter->destroy)
            jitter->destroy(jitter->packets[i].data);
         else
            speex_free(jitter->packets[i].data);
         jitter->packets[i].data = NULL;
      }
   }
   /* Timestamp is actually undefined at this point */
   jitter->pointer_timestamp = 0;
   jitter->current_timestamp = 0;
   jitter->next_stop = 0;
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

void update_histogram(JitterBuffer *jitter, spx_int32_t arrival_margin)
{
   int i;
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
}

/** Put one packet into the jitter buffer */
void jitter_buffer_put(JitterBuffer *jitter, const JitterBufferPacket *packet)
{
   int i,j;
   int late;
   spx_int32_t arrival_margin;
   spx_int32_t arrival_time;
   /*fprintf (stderr, "put packet %d %d\n", timestamp, span);*/
   
   /* Syncing on the first packet to arrive */
   if (jitter->reset_state)
   {
      jitter->reset_state=0;
      jitter->pointer_timestamp = packet->timestamp;
      jitter->current_timestamp = packet->timestamp;
      jitter->next_stop = packet->timestamp;
      /*fprintf(stderr, "reset to %d\n", timestamp);*/
   }
   
   /* Cleanup buffer (remove old packets that weren't played) */
   for (i=0;i<SPEEX_JITTER_MAX_BUFFER_SIZE;i++)
   {
      /* Make sure we don't discard a "just-late" packet in case we want to play it next (if we interpolate). */
      if (jitter->packets[i].data && LE32(jitter->packets[i].timestamp + jitter->packets[i].span, jitter->pointer_timestamp))
      {
         /*fprintf (stderr, "cleaned (not played)\n");*/
         if (jitter->destroy)
            jitter->destroy(jitter->packets[i].data);
         else
            speex_free(jitter->packets[i].data);
         jitter->packets[i].data = NULL;
      }
   }

   fprintf(stderr, "arrival: %d %d %d\n", packet->timestamp, jitter->next_stop, jitter->pointer_timestamp);
   /* Check if packet is late (could still be useful though) */
   if (LT32(packet->timestamp, jitter->next_stop))
   {
      fprintf(stderr, "late by %d\n", jitter->next_stop - packet->timestamp);
      
      /* The arrival margin is how much in advance (or in this case late) the packet it (in resolution units) */
      arrival_margin = (((spx_int32_t)packet->timestamp) - ((spx_int32_t)jitter->next_stop))/jitter->resolution - jitter->buffer_margin;
   
      fprintf(stderr, "put arrival_margin = %d\n", arrival_margin);

      update_histogram(jitter, arrival_margin);
      late = 1;
   } else {
      late = 0;
   }
   
   /* Only insert the packet if it's not hopelessly late (i.e. totally useless) */
   if (GE32(packet->timestamp+packet->span, jitter->pointer_timestamp))
   {

      /*Find an empty slot in the buffer*/
      for (i=0;i<SPEEX_JITTER_MAX_BUFFER_SIZE;i++)
      {
         if (jitter->packets[i].data==NULL)
            break;
      }
      
      /*No place left in the buffer, need to make room for it by discarding the oldest packet */
      if (i==SPEEX_JITTER_MAX_BUFFER_SIZE)
      {
         int earliest=jitter->packets[0].timestamp;
         i=0;
         for (j=1;j<SPEEX_JITTER_MAX_BUFFER_SIZE;j++)
         {
            if (!jitter->packets[i].data || LT32(jitter->packets[j].timestamp,earliest))
            {
               earliest = jitter->packets[j].timestamp;
               i=j;
            }
         }
         if (jitter->destroy)
            jitter->destroy(jitter->packets[i].data);
         else
            speex_free(jitter->packets[i].data);
         jitter->packets[i].data=NULL;
         if (jitter->lost_count>20)
         {
            jitter_buffer_reset(jitter);
         }
         /*fprintf (stderr, "Buffer is full, discarding earliest frame %d (currently at %d)\n", timestamp, jitter->pointer_timestamp);*/      
      }
   
      /* Copy packet in buffer */
      if (jitter->destroy)
      {
         jitter->packets[i].data = packet->data;
      } else {
         jitter->packets[i].data=(char*)speex_alloc(packet->len);
         for (j=0;j<packet->len;j++)
            jitter->packets[i].data[j]=packet->data[j];
      }
      jitter->packets[i].timestamp=packet->timestamp;
      jitter->packets[i].span=packet->span;
      jitter->packets[i].len=packet->len;
      jitter->packets[i].user_data=packet->user_data;
      if (late)
         jitter->arrival[i] = 0;
      else
         jitter->arrival[i] = jitter->next_stop;
   }
   
   
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
   
   if (jitter->current_timestamp + jitter->sub_clock > jitter->pointer_timestamp)
      speex_warning("something's wrong with the time");

   jitter->sub_clock = -1;
   jitter->current_timestamp = jitter->pointer_timestamp;
   jitter->last_returned_timestamp = jitter->pointer_timestamp;
         
   if (jitter->interp_requested)
   {
      jitter->interp_requested = 0;
      if (start_offset)
         *start_offset = 0;
      packet->timestamp = jitter->pointer_timestamp;
      packet->span = jitter->delay_step;
      jitter->pointer_timestamp += jitter->delay_step;
      packet->len = 0;
      fprintf (stderr, "Deferred interpolate\n");
      
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
      if (jitter->packets[i].data && jitter->packets[i].timestamp==jitter->pointer_timestamp && GE32(jitter->packets[i].timestamp+jitter->packets[i].span,jitter->pointer_timestamp+desired_span))
         break;
   }
   
   /* If no match, try for an "older" packet that still spans (fully) the current chunk */
   if (i==SPEEX_JITTER_MAX_BUFFER_SIZE)
   {
      for (i=0;i<SPEEX_JITTER_MAX_BUFFER_SIZE;i++)
      {
         if (jitter->packets[i].data && LE32(jitter->packets[i].timestamp, jitter->pointer_timestamp) && GE32(jitter->packets[i].timestamp+jitter->packets[i].span,jitter->pointer_timestamp+desired_span))
            break;
      }
   }
   
   /* If still no match, try for an "older" packet that spans part of the current chunk */
   if (i==SPEEX_JITTER_MAX_BUFFER_SIZE)
   {
      for (i=0;i<SPEEX_JITTER_MAX_BUFFER_SIZE;i++)
      {
         if (jitter->packets[i].data && LE32(jitter->packets[i].timestamp, jitter->pointer_timestamp) && GT32(jitter->packets[i].timestamp+jitter->packets[i].span,jitter->pointer_timestamp))
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
         if (jitter->packets[i].data && LT32(jitter->packets[i].timestamp,jitter->pointer_timestamp+desired_span) && GE32(jitter->packets[i].timestamp,jitter->pointer_timestamp))
         {
            if (!found || LT32(jitter->packets[i].timestamp,best_time) || (jitter->packets[i].timestamp==best_time && GT32(jitter->packets[i].span,best_span)))
            {
               best_time = jitter->packets[i].timestamp;
               best_span = jitter->packets[i].span;
               besti = i;
               found = 1;
            }
         }
      }
      if (found)
      {
         i=besti;
         incomplete = 1;
         /*fprintf (stderr, "incomplete: %d %d %d %d\n", jitter->packets[i].timestamp, jitter->pointer_timestamp, chunk_size, jitter->packets[i].span);*/
      }
   }

   /* If we find something */
   if (i!=SPEEX_JITTER_MAX_BUFFER_SIZE)
   {
      
      
      /* We (obviously) haven't lost this packet */
      jitter->lost_count = 0;
      jitter->loss_rate = .999*jitter->loss_rate;
      
      /* In this case, 0 isn't as a valid timestamp */
      if (jitter->arrival[i] != 0)
      {
         spx_int32_t arrival_margin;
         fprintf(stderr, "early by %d\n", jitter->packets[i].timestamp - jitter->arrival[i]);
         
         /* The arrival margin is how much in advance (or in this case late) the packet it (in resolution units) */
         arrival_margin = (((spx_int32_t)jitter->packets[i].timestamp) - ((spx_int32_t)jitter->arrival[i]))/jitter->resolution - jitter->buffer_margin;
   
         fprintf(stderr, "get arrival_margin = %d\n", arrival_margin);

         update_histogram(jitter, arrival_margin);

      }
      
      
      /* FIXME: Check for potential overflow */
      packet->len = jitter->packets[i].len;
      /* Copy packet */
      if (jitter->destroy)
      {
         packet->data = jitter->packets[i].data;
      } else {
         for (j=0;j<packet->len;j++)
            packet->data[j] = jitter->packets[i].data[j];
         /* Remove packet */
         speex_free(jitter->packets[i].data);
      }
      jitter->packets[i].data = NULL;
      /* Set timestamp and span (if requested) */
      if (start_offset)
         *start_offset = (spx_int32_t)jitter->packets[i].timestamp-(spx_int32_t)jitter->pointer_timestamp;
      
      packet->timestamp = jitter->packets[i].timestamp;
      jitter->last_returned_timestamp = packet->timestamp;
      
      packet->span = jitter->packets[i].span;
      packet->user_data = jitter->packets[i].user_data;
      /* Point to the end of the current packet */
      jitter->pointer_timestamp = jitter->packets[i].timestamp+jitter->packets[i].span;

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
      fprintf (stderr, "Forced to interpolate\n");

      /*fprintf (stderr, "i");*/
      /*fprintf (stderr, "interpolate (getting some slack)\n");*/
   }

   return JITTER_BUFFER_MISSING;

}

int jitter_buffer_get_another(JitterBuffer *jitter, JitterBufferPacket *packet)
{
   int i, j;
   for (i=0;i<SPEEX_JITTER_MAX_BUFFER_SIZE;i++)
   {
      if (jitter->packets[i].data && jitter->packets[i].timestamp==jitter->last_returned_timestamp)
         break;
   }
   if (i!=SPEEX_JITTER_MAX_BUFFER_SIZE)
   {
      /* Copy packet */
      packet->len = jitter->packets[i].len;
      if (jitter->destroy)
      {
         packet->data = jitter->packets[i].data;
      } else {
         for (j=0;j<packet->len;j++)
            packet->data[j] = jitter->packets[i].data[j];
         /* Remove packet */
         speex_free(jitter->packets[i].data);
      }
      jitter->packets[i].data = NULL;
      packet->timestamp = jitter->packets[i].timestamp;
      packet->span = jitter->packets[i].span;
      packet->user_data = jitter->packets[i].user_data;
      return JITTER_BUFFER_OK;
   } else {
      packet->data = NULL;
      packet->len = 0;
      packet->span = 0;
      return JITTER_BUFFER_MISSING;
   }
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
   jitter->sub_clock += jitter->resolution;
   jitter->next_stop = jitter->pointer_timestamp;
}

void jitter_buffer_remaining_span(JitterBuffer *jitter, spx_uint32_t rem)
{
   jitter->next_stop = jitter->pointer_timestamp - rem;
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
      fprintf (stderr, "Decision to interpolate\n");
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
      fprintf (stderr, "Decision to drop\n");
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
            if (jitter->packets[i].data && LE32(jitter->pointer_timestamp, jitter->packets[i].timestamp))
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

