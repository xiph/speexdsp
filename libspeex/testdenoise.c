#include "speex_denoise.h"
#include <stdio.h>

#define NN 240

int main()
{
   short in[NN];
   short out[NN];
   float x[NN];
   int i;
   SpeexDenoiseState *st;

   st = speex_denoise_state_init(NN);
   while (1)
   {
      int vad;
      fread(in, sizeof(short), NN, stdin);
      if (feof(stdin))
         break;
      for (i=0;i<NN;i++)
         x[i]=in[i];
      
      vad = speex_denoise(st, x);
      for (i=0;i<NN;i++)
         out[i]=x[i];
      /*fprintf (stderr, "%d\n", vad);*/
      fwrite(out, sizeof(short), NN, stdout);
   }
   speex_denoise_state_destroy(st);
   return 0;
}
