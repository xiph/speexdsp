#include <speex.h>
#include <stdio.h>
#include <stdlib.h>

#define FRAME_SIZE 160
int main(int argc, char **argv)
{
   char *inFile;
   FILE *fin;
   short in[FRAME_SIZE];
   float input[FRAME_SIZE];
   char cbits[200];
   int nbBytes;
   void *state;
   SpeexBits bits;
   int i, tmp;


   state = speex_encoder_init(&speex_nb_mode);

   tmp=8;
   speex_encoder_ctl(state, SPEEX_SET_QUALITY, &tmp);

   inFile = argv[1];
   fin = fopen(inFile, "r");

   speex_bits_init(&bits);
   while (1)
   {
      fread(in, sizeof(short), FRAME_SIZE, fin);
      if (feof(fin))
         break;
      for (i=0;i<FRAME_SIZE;i++)
         input[i]=in[i];
      speex_bits_reset(&bits);

      speex_encode(state, input, &bits);
      nbBytes = speex_bits_write(&bits, cbits, 200);

      fwrite(&nbBytes, sizeof(int), 1, stdout);
      fwrite(cbits, 1, nbBytes, stdout);
      speex_bits_rewind(&bits);
      
   }
   
   speex_encoder_destroy(state);
   speex_bits_destroy(&bits);
   fclose(fin);
   return 0;
}
