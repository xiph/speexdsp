#include "speex.h"
#include <stdio.h>
#include <stdlib.h>

#define FRAME_SIZE 256
#include <math.h>
int main(int argc, char **argv)
{
   char *inFile, *outFile, *bitsFile;
   FILE *fin, *fout, *fbits=NULL;
   short in[FRAME_SIZE];
   float input[FRAME_SIZE], bak[FRAME_SIZE], bak2[FRAME_SIZE];
   char cbits[200];
   int nbBits;
   int i;
   EncState st;
   DecState dec;
   FrameBits bits;

   for (i=0;i<FRAME_SIZE;i++)
      bak2[i]=0;
   encoder_init(&st, &wb_mode);
   decoder_init(&dec, &wb_mode);
   if (argc != 4 && argc != 3)
   {
      fprintf (stderr, "Usage: encode [in file] [out file] [bits file]\nargc = %d", argc);
      exit(1);
   }
   inFile = argv[1];
   fin = fopen(inFile, "r");
   outFile = argv[2];
   fout = fopen(outFile, "w");
   if (argc==4)
   {
      bitsFile = argv[3];
      fbits = fopen(bitsFile, "w");
   }
   frame_bits_init(&bits);
   while (!feof(fin))
   {
      fread(in, sizeof(short), FRAME_SIZE, fin);
      for (i=0;i<FRAME_SIZE;i++)
         bak[i]=input[i]=in[i];
      frame_bits_reset(&bits);
      encode(&st, input, &bits);
      nbBits = frame_bits_write(&bits, cbits, 200);
      printf ("Encoding frame in %d bits\n", nbBits*8);
      if (argc==4)
         fwrite(cbits, 1, nbBits, fbits);
      {
         float enoise=0, esig=0, snr;
         for (i=0;i<FRAME_SIZE;i++)
         {
            enoise+=(bak2[i]-input[i])*(bak2[i]-input[i]);
            esig += bak2[i]*bak2[i];
         }
         snr = 10*log10((esig+1)/(enoise+1));
         printf ("real SNR = %f\n", snr);
      }

      /* Save the bits here */
      frame_bits_reset(&bits);
      for (i=0;i<FRAME_SIZE;i++)
         in[i]=input[i];
      for (i=0;i<FRAME_SIZE;i++)
         bak2[i]=bak[i];
      fwrite(in, sizeof(short), FRAME_SIZE, fout);
   }
   
   encoder_destroy(&st);
   decoder_destroy(&dec);
   return 1;
}
