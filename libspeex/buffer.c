
#include "os_support.h"
#include "misc.h"

struct SpeexBuffer_;
typedef struct SpeexBuffer_ SpeexBuffer;


struct SpeexBuffer_ {
   void *data;
   int   size;
   int   read_ptr;
   int   write_ptr;
   int   available;
};

SpeexBuffer *speex_buffer_init(int size)
{
   SpeexBuffer *st = speex_alloc(sizeof(SpeexBuffer));
   st->data = speex_alloc(size);
   st->size = size;
   st->read_ptr = 0;
   st->write_ptr = 0;
   st->available = 0;
   return st;
}

int speex_buffer_write(SpeexBuffer *st, void *data, int len)
{
   int end;
   int end1;
   if (len > st->size)
   {
      data += len-st->size;
      len = st->size;
   }
   end = st->write_ptr + len;
   end1 = end;
   if (end1 > st->size)
      end1 = st->size;
   speex_move(st->data + st->write_ptr, data, end1 - st->write_ptr);
   if (end > st->size)
   {
      end -= st->size;
      speex_move(st->data, data+end1 - st->write_ptr, end);
   }
   st->available += len;
   if (st->available > st->size)
   {
      st->available = st->size;
      st->read_ptr = st->write_ptr;
   }
   return len;
}

/*int speex_buffer_writezeros(SpeexBuffer *st, int len)
{
   return len;
}*/

int speex_buffer_read(SpeexBuffer *st, void *data, int len)
{
   int end, end1;
   if (len > st->available)
   {
      speex_memset(data+st->available, 0, st->size-st->available);
      len = st->available;
   }
   end = st->read_ptr + len;
   end1 = end;
   if (end1 > st->size)
      end1 = st->size;
   speex_move(data, st->data + st->read_ptr, end1 - st->read_ptr);

   if (end > st->size)
   {
      end -= st->size;
      speex_move(data+end1 - st->read_ptr, st->data, end);
   }
   st->available -= len;
   return len;
}

int speex_buffer_get_available(SpeexBuffer *st)
{
   return st->available;
}

int speex_buffer_resize(SpeexBuffer *st, void *data, int len)
{
   return len;
}
