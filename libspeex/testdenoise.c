#include "speex_preprocess.h"
#include <stdio.h>

#define NN 160

int main()
{
   short in[NN];
   short out[NN];
   float x[NN];
   int i;
   SpeexPreprocessState *st;
   int count=0;
   float f;

   st = speex_preprocess_state_init(NN, 8000);
   i=1;
   speex_preprocess_ctl(st, SPEEX_PREPROCESS_SET_DENOISE, &i);
   i=0;
   speex_preprocess_ctl(st, SPEEX_PREPROCESS_SET_AGC, &i);
   f=8000;
   speex_preprocess_ctl(st, SPEEX_PREPROCESS_SET_AGC_LEVEL, &f);
   while (1)
   {
      int vad;
      fread(in, sizeof(short), NN, stdin);
      if (feof(stdin))
         break;
      for (i=0;i<NN;i++)
         x[i]=in[i];
      vad = speex_preprocess(st, x, NULL);
      for (i=0;i<NN;i++)
         out[i]=x[i];
      /*fprintf (stderr, "%d\n", vad);*/
      fwrite(out, sizeof(short), NN, stdout);
      count++;
   }
   speex_preprocess_state_destroy(st);
   return 0;
}
