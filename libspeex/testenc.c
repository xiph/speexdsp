#include "speex.h"
#include <stdio.h>
#include <stdlib.h>

#define FRAME_SIZE 128

int main(int argc, char **argv)
{
   char *inFile, *outFile;
   FILE *fin, *fout;
   short in[FRAME_SIZE];
   float input[FRAME_SIZE];
   int i;
   EncState st;

   encoder_init(&st);
   if (argc != 3)
   {
      fprintf (stderr, "Usage: encode [in file] [out file]\n");
      exit(1);
   }
   inFile = argv[1];
   fin = fopen(inFile, "r");
   outFile = argv[2];
   fout = fopen(outFile, "w");
   while (!feof(fin))
   {
      fread(in, sizeof(short), FRAME_SIZE, fin);
      for (i=0;i<FRAME_SIZE;i++)
         input[i]=in[i];
      encode(&st, input, NULL, NULL);
   }
   encoder_destroy(&st);
   return 1;
}
