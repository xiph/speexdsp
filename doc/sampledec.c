#include <speex.h>
#include <stdio.h>
#include <stdlib.h>

#define FRAME_SIZE 160
int main(int argc, char **argv)
{
   char *outFile;
   FILE *fout;
   short out[FRAME_SIZE];
   float output[FRAME_SIZE];
   char cbits[200];
   int nbBytes;
   void *state;
   SpeexBits bits;
   int i, tmp;


   state = speex_decoder_init(&speex_nb_mode);

   tmp=1;
   speex_decoder_ctl(state, SPEEX_SET_ENH, &tmp);

   outFile = argv[1];
   fout = fopen(outFile, "w");

   speex_bits_init(&bits);
   while (1)
   {
      fread(&nbBytes, sizeof(int), 1, stdin);
      fprintf (stderr, "nbBytes: %d\n", nbBytes);
      if (feof(stdin))
         break;

      fread(cbits, 1, nbBytes, stdin);
      speex_bits_read_from(&bits, cbits, nbBytes);

      speex_decode(state, &bits, output);

      for (i=0;i<FRAME_SIZE;i++)
         out[i]=output[i];

      fwrite(out, sizeof(short), FRAME_SIZE, fout);
   }
   
   speex_encoder_destroy(state);
   speex_bits_destroy(&bits);
   fclose(fout);
   return 0;
}
