#include "denoise.h"
#include <stdio.h>

#define NN 160

int main()
{
   short in[NN];
   short out[NN];
   float x[NN];
   int i;
   DenoiseState *st;

   st = denoise_state_init(NN);
   while (1)
   {
      fread(in, sizeof(short), NN, stdin);
      if (feof(stdin))
         break;
      for (i=0;i<NN;i++)
         x[i]=in[i];
      
      denoise(st, x);
      for (i=0;i<NN;i++)
         out[i]=x[i];
      
      fwrite(out, sizeof(short), NN, stdout);
   }
   denoise_state_destroy(st);
   return 0;
}
