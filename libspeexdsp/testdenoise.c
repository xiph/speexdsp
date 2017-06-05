#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "speex/speex_preprocess.h"
#include <stdio.h>

#define NN 160*6

int main(int argc, char *argv[])
{
   short in[NN];
   int i;
   SpeexPreprocessState *st;
   int count=0;
   float f;

   freopen(argv[1], "rb", stdin);
   freopen(argv[2], "wb", stdout);

   st = speex_preprocess_state_init(NN, 8000*6);
   i=1;
   speex_preprocess_ctl(st, SPEEX_PREPROCESS_SET_DENOISE, &i);
   i=0;
   speex_preprocess_ctl(st, SPEEX_PREPROCESS_SET_AGC, &i);
   i=8000;
   speex_preprocess_ctl(st, SPEEX_PREPROCESS_SET_AGC_LEVEL, &i);
   i=0;
   speex_preprocess_ctl(st, SPEEX_PREPROCESS_SET_DEREVERB, &i);
   f=.0;
   speex_preprocess_ctl(st, SPEEX_PREPROCESS_SET_DEREVERB_DECAY, &f);
   f=.0;
   speex_preprocess_ctl(st, SPEEX_PREPROCESS_SET_DEREVERB_LEVEL, &f);

   while (1)
   {
      int vad;
      fread(in, sizeof(short), NN, stdin);
      if (feof(stdin))
         break;

      vad = speex_preprocess_run(st, in);
      /*fprintf (stderr, "%d\n", vad);*/
      fwrite(in, sizeof(short), NN, stdout);
      count++;
   }

   speex_preprocess_state_destroy(st);
   return 0;
}
