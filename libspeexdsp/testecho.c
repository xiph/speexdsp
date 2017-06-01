#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "speex/speex_echo.h"
#include "speex/speex_preprocess.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


#define NN 128
#define TAIL 1024

int main(int argc, char **argv)
{
   short echo_buf[NN], ref_buf[NN], e_buf[NN];
   const int sampleRate = 8000;

   if (argc != 4)
   {
      fprintf(stderr, "testecho mic_signal.sw speaker_signal.sw output.sw\n");
      exit(1);
   }

   FILE *echo_fd = fopen(argv[2], "rb");//FAR-end file(to speaker)
   FILE *ref_fd  = fopen(argv[1],  "rb");//mic capture file
   FILE *e_fd  = fopen(argv[3], "wb");//output file

   SpeexEchoState *st = speex_echo_state_init(NN, TAIL);
   SpeexPreprocessState *den = speex_preprocess_state_init(NN, sampleRate);

   speex_echo_ctl(st, SPEEX_ECHO_SET_SAMPLING_RATE, &sampleRate);
   speex_preprocess_ctl(den, SPEEX_PREPROCESS_SET_ECHO_STATE, st);

   while (!feof(ref_fd) && !feof(echo_fd))
   {
      fread(ref_buf, sizeof(short), NN, ref_fd);
      fread(echo_buf, sizeof(short), NN, echo_fd);

      speex_echo_cancellation(st, ref_buf, echo_buf, e_buf);
      speex_preprocess_run(den, e_buf);

      fwrite(e_buf, sizeof(short), NN, e_fd);
   }

   speex_echo_state_destroy(st);
   speex_preprocess_state_destroy(den);
   fclose(e_fd);
   fclose(echo_fd);
   fclose(ref_fd);

   return 0;
}
