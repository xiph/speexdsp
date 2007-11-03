/* Copyright (C) 2002 Jean-Marc Valin */
/**
   @file speex_jitter.h
   @brief Adaptive jitter buffer for Speex
*/
/*
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

#ifndef SPEEX_JITTER_H
#define SPEEX_JITTER_H
/** @defgroup JitterBuffer JitterBuffer: Adaptive jitter buffer
 *  This is the jitter buffer that reorders UDP/RTP packets and adjusts the buffer size
 * to maintain good quality and low latency.
 *  @{
 */

#include "speex/speex_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Generic adaptive jitter buffer state */
struct JitterBuffer_;

/** Generic adaptive jitter buffer state */
typedef struct JitterBuffer_ JitterBuffer;

/** Definition of an incoming packet */
typedef struct _JitterBufferPacket JitterBufferPacket;

/** Definition of an incoming packet */
struct _JitterBufferPacket {
   char        *data;       /**< Data bytes contained in the packet */
   spx_uint32_t len;        /**< Length of the packet in bytes */
   spx_uint32_t timestamp;  /**< Timestamp for the packet */
   spx_uint32_t span;       /**< Time covered by the packet (same units as timestamp) */
   spx_uint32_t user_data;  /**< Put whatever data you like here (it's ignored by the jitter buffer) */
};

/** Packet has been retrieved */
#define JITTER_BUFFER_OK 0
/** Packet is missing */
#define JITTER_BUFFER_MISSING 1
/** Packet is incomplete (does not cover the entire tick */
#define JITTER_BUFFER_INCOMPLETE 2
/** There was an error in the jitter buffer */
#define JITTER_BUFFER_INTERNAL_ERROR -1
/** Invalid argument */
#define JITTER_BUFFER_BAD_ARGUMENT -2


/** Set minimum amount of extra buffering required (margin) */
#define JITTER_BUFFER_SET_MARGIN 0
/** Get minimum amount of extra buffering required (margin) */
#define JITTER_BUFFER_GET_MARGIN 1
/* JITTER_BUFFER_SET_AVAILABLE_COUNT wouldn't make sense */

/** Get the amount of available packets currently buffered */
#define JITTER_BUFFER_GET_AVAILABLE_COUNT 3
/** Included because of an early misspelling (will remove in next release) */
#define JITTER_BUFFER_GET_AVALIABLE_COUNT 3

/**  */
#define JITTER_BUFFER_SET_DESTROY_CALLBACK 4
/**  */
#define JITTER_BUFFER_GET_DESTROY_CALLBACK 5

/**  */
#define JITTER_BUFFER_SET_DELAY_STEP 6
/**  */
#define JITTER_BUFFER_GET_DELAY_STEP 7


#define JITTER_BUFFER_ADJUST_INTERPOLATE -1
#define JITTER_BUFFER_ADJUST_OK 0
#define JITTER_BUFFER_ADJUST_DROP 1

/** Initialises jitter buffer 
 * 
 * @param tick Number of samples per "tick", i.e. the time period of the elements that will be retrieved
 * @return Newly created jitter buffer state
 */
JitterBuffer *jitter_buffer_init(int tick);

/** Restores jitter buffer to its original state 
 * 
 * @param jitter Jitter buffer state
 */
void jitter_buffer_reset(JitterBuffer *jitter);

/** Destroys jitter buffer 
 * 
 * @param jitter Jitter buffer state
 */
void jitter_buffer_destroy(JitterBuffer *jitter);

/** Put one packet into the jitter buffer
 * 
 * @param jitter Jitter buffer state
 * @param packet Incoming packet
*/
void jitter_buffer_put(JitterBuffer *jitter, const JitterBufferPacket *packet);

/** Get one packet from the jitter buffer
 * 
 * @param jitter Jitter buffer state
 * @param packet Returned packet
 * @param desired_span Number of samples (or units) we wish to get from the buffer (no guarantee)
 * @param current_timestamp Timestamp for the returned packet 
*/
int jitter_buffer_get(JitterBuffer *jitter, JitterBufferPacket *packet, spx_int32_t desired_span, spx_int32_t *start_offset);

/** Used right after jitter_buffer_get() to obtain another packet that would have the same timestamp.
 * This is mainly useful for media where a single "frame" can be split into several packets.
 * 
 * @param jitter Jitter buffer state
 * @param packet Returned packet
 */
int jitter_buffer_get_another(JitterBuffer *jitter, JitterBufferPacket *packet);

/** Get pointer timestamp of jitter buffer
 * 
 * @param jitter Jitter buffer state
*/
int jitter_buffer_get_pointer_timestamp(JitterBuffer *jitter);

/** Advance by one tick
 * 
 * @param jitter Jitter buffer state
*/
void jitter_buffer_tick(JitterBuffer *jitter);

/** Telling the jitter buffer about the remaining data in the application buffer
 * @param jitter Jitter buffer state
 * @param rem Amount of data buffered by the application (timestamp units)
 */
void jitter_buffer_remaining_span(JitterBuffer *jitter, spx_uint32_t rem);

/** Used like the ioctl function to control the jitter buffer parameters
 * 
 * @param jitter Jitter buffer state
 * @param request ioctl-type request (one of the JITTER_BUFFER_* macros)
 * @param ptr Data exchanged to-from function
 * @return 0 if no error, -1 if request in unknown
*/
int jitter_buffer_ctl(JitterBuffer *jitter, int request, void *ptr);

int jitter_buffer_update_delay(JitterBuffer *jitter, JitterBufferPacket *packet, spx_int32_t *start_offset);

/* @} */

#ifdef __cplusplus
}
#endif

#endif
