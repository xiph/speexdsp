#include "speex.h"
#include <stdio.h>
#include <stdlib.h>

#define FRAME_SIZE 160
#include <math.h>
int main(int argc, char **argv)
{
   char *outFile, *bitsFile;
   FILE *fout, *fbits=NULL;
   short out[FRAME_SIZE];
   float output[FRAME_SIZE];
   char cbits[200];
   int i;
   DecState dec;
   FrameBits bits;

   decoder_init(&dec, &mp_nb_mode);
   if (argc != 3)
   {
      fprintf (stderr, "Usage: encode [bits file] [out file]\nargc = %d", argc);
      exit(1);
   }
   bitsFile = argv[1];
   fbits = fopen(bitsFile, "r");
   outFile = argv[2];
   fout = fopen(outFile, "w");
   frame_bits_init(&bits);
   while (!feof(fbits))
   {
      fread(cbits, 1, 37, fbits);
      frame_bits_reset(&bits);
      frame_bits_init_from(&bits, cbits, 37);
      decode(&dec, &bits, output);
      for (i=0;i<FRAME_SIZE;i++)
      {
         if (output[i]>32000)
            output[i]=32000;
         else if (output[i]<-32000)
            output[i]=-32000;
      }
      for (i=0;i<FRAME_SIZE;i++)
         out[i]=output[i];
      fwrite(out, sizeof(short), FRAME_SIZE, fout);
   }
   decoder_destroy(&dec);
   return 1;
}
