#include "speex.h"
#include <stdio.h>
#include <stdlib.h>

#define FRAME_SIZE 160

int main(int argc, char **argv)
{
   char *inFile, *outFile;
   FILE *fin, *fout;
   short in[FRAME_SIZE];
   float input[FRAME_SIZE];
   int i;
   EncState st;
   FrameBits bits;

   encoder_init(&st, &nb_mode);
   if (argc != 3)
   {
      fprintf (stderr, "Usage: encode [in file] [out file]\n");
      exit(1);
   }
   inFile = argv[1];
   fin = fopen(inFile, "r");
   outFile = argv[2];
   fout = fopen(outFile, "w");
   frame_bits_init(&bits);
   while (!feof(fin))
   {
      fread(in, sizeof(short), FRAME_SIZE, fin);
      for (i=0;i<FRAME_SIZE;i++)
         input[i]=in[i];
      encode(&st, input, &bits);
      /* Save the bits here */
      frame_bits_reset(&bits);
      for (i=0;i<FRAME_SIZE;i++)
         in[i]=input[i];
      fwrite(in, sizeof(short), FRAME_SIZE, fout);
   }
   
   encoder_destroy(&st);
   return 1;
}
