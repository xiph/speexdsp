#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "speex/speex_jitter.h"
#include <stdio.h>
#include <stdlib.h>

union jbpdata {
  unsigned int idx;
  unsigned char data[4];
};

void synthIn(JitterBufferPacket *in, union jbpdata *d, int idx, int span) {
  d->idx = idx;

  in->data = (char*)d->data;
  in->len = sizeof(*d);
  in->timestamp = idx * 10;
  in->span = span * 10;
  in->sequence = idx;
  in->user_data = 0;
}

static int errors;

void jitterFill(JitterBuffer *jb) {
   char buffer[65536];
   JitterBufferPacket in, out;
   union jbpdata d;
   int i;

   out.data = buffer;

   jitter_buffer_reset(jb);

   for(i=0;i<100;++i) {
     synthIn(&in, &d, i, 1);
     jitter_buffer_put(jb, &in);

     out.len = 65536;
     if (jitter_buffer_get(jb, &out, 10, NULL) != JITTER_BUFFER_OK) {
       printf("Fill test failed iteration %d\n", i);
       errors++;
     }
     if (out.timestamp != i * 10) {
       printf("Fill test expected %d got %d\n", i*10, out.timestamp);
       errors++;
     }
     jitter_buffer_tick(jb);
   }
}

/* Packets for the destroy callback test. Nothing is really freed until the
   end, so that a packet released twice is reported instead of crashing. */
#define MAX_OWNED 16
static char *owned_data[MAX_OWNED];
static int owned_released[MAX_OWNED];
static int nb_owned;
static int nb_destroyed;
static int owned_errors;

static void ownedFail(const char *msg) {
   printf("Destroy callback test: %s\n", msg);
   owned_errors++;
}

static void ownedRelease(void *data) {
   int i;
   for (i=0;i<nb_owned;i++) {
      if (owned_data[i] == data)
         break;
   }
   if (i == nb_owned)
      ownedFail("released an unknown packet");
   else if (owned_released[i])
      ownedFail("packet released twice");
   else
      owned_released[i] = 1;
}

static void ownedDestroy(void *data) {
   nb_destroyed++;
   ownedRelease(data);
}

static char *ownedPut(JitterBuffer *jb, unsigned int timestamp) {
   JitterBufferPacket in;
   char *data = malloc(4);

   owned_data[nb_owned] = data;
   owned_released[nb_owned++] = 0;
   in.data = data;
   in.len = 4;
   in.timestamp = timestamp;
   in.span = 10;
   in.sequence = 0;
   in.user_data = 0;
   jitter_buffer_put(jb, &in);
   return data;
}

/* With a destroy callback, the buffer keeps the pointers it is given: packets
   that are returned belong to the caller, and all the others must go through
   the callback, exactly once in both cases. */
static int jitterDestroyCallback(void) {
   JitterBufferPacket out;
   char *first;
   int i, destroyed;
   JitterBuffer *jb = jitter_buffer_init(10);

   jitter_buffer_ctl(jb, JITTER_BUFFER_SET_DESTROY_CALLBACK, (void*)ownedDestroy);

   first = ownedPut(jb, 0);
   ownedPut(jb, 10);
   ownedPut(jb, 20);

   if (jitter_buffer_get(jb, &out, 10, NULL) != JITTER_BUFFER_OK || out.data != first)
      ownedFail("did not get the first packet back");
   else
      ownedRelease(out.data);
   if (jitter_buffer_get_another(jb, &out) == JITTER_BUFFER_OK)
      ownedFail("got another packet when there was only one");
   jitter_buffer_tick(jb);

   /* Move well past the packets that were put */
   for (i=0;i<9;i++) {
      if (jitter_buffer_get(jb, &out, 10, NULL) == JITTER_BUFFER_OK)
         ownedRelease(out.data);
      jitter_buffer_tick(jb);
   }

   /* Too late to be inserted: there is no way to tell, so it must be destroyed */
   destroyed = nb_destroyed;
   ownedPut(jb, 0);
   if (nb_destroyed != destroyed+1 || !owned_released[nb_owned-1])
      ownedFail("late packet was not destroyed");

   /* Still queued when the buffer goes away */
   ownedPut(jb, jitter_buffer_get_pointer_timestamp(jb));
   ownedPut(jb, jitter_buffer_get_pointer_timestamp(jb)+10);
   jitter_buffer_destroy(jb);

   for (i=0;i<nb_owned;i++) {
      if (!owned_released[i])
         ownedFail("packet leaked");
      free(owned_data[i]);
   }
   if (!owned_errors)
      printf("Destroy callback: %d packets, %d destroyed\n", nb_owned, nb_destroyed);
   return owned_errors;
}

int main()
{
   char buffer[65536];
   JitterBufferPacket in, out;
   union jbpdata d;
   int i;

   JitterBuffer *jb = jitter_buffer_init(10);

   out.data = buffer;

   /* Frozen sender case */
   jitterFill(jb);
   for(i=0;i<100;++i) {
     out.len = 65536;
     jitter_buffer_get(jb, &out, 10, NULL);
     jitter_buffer_tick(jb);
   }
   synthIn(&in, &d, 100, 1);
   jitter_buffer_put(jb, &in);
   out.len = 65536;
   if (jitter_buffer_get(jb, &out, 10, NULL) != JITTER_BUFFER_OK) {
     printf("Failed frozen sender resynchronize\n");
     errors++;
   } else {
     printf("Frozen sender: Jitter %d\n", out.timestamp - 100*10);
   }
   jitter_buffer_destroy(jb);

   errors += jitterDestroyCallback();
   return errors != 0;
}
