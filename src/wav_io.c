/* Copyright (C) 2002 Jean-Marc Valin 
   File: wav_io.c
   Routines to handle wav (RIFF) headers

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.
   
   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <stdio.h>
#include <string.h>
#include "misc.h"

int read_wav_header(FILE *file, int *rate, int *channels, int *format, int *size)
{
   char ch[5];
   int itmp;
   short stmp;
   int bpersec;
   short balign;

   ch[4]=0;
   fread(ch, 1, 4, file);
   if (strcmp(ch, "RIFF")!=0)
   {
      fseek(file, 0, SEEK_SET);
      return 0;
   }

   fread(&itmp, 4, 1, file);
   *size = le_int(itmp-36);

   fread(ch, 1, 4, file);
   if (strcmp(ch, "WAVE")!=0)
   {
      fprintf (stderr, "RIFF file is not a WAVE file\n");
      return -1;
   }

   fread(ch, 1, 4, file);
   if (strcmp(ch, "fmt ")!=0)
   {
      fprintf (stderr, "Corrupted WAVE file: no \"fmt \"\n");
      return -1;
   }
   
   fread(&itmp, 4, 1, file);
   itmp = le_int(itmp);
   if (itmp!=16)
   {
      fprintf (stderr, "Unsupported WAVE file fmt chunk, not PCM?\n");
      return -1;
   }

   fread(&stmp, 2, 1, file);
   stmp = le_short(stmp);
   if (stmp!=1)
   {
      fprintf (stderr, "Only PCM encoding is supported\n");
      return -1;
   }

   fread(&stmp, 2, 1, file);
   stmp = le_short(stmp);
   *channels = stmp;
   
   if (stmp>1)
   {
      fprintf (stderr, "Only mono supported for now\n");
      return -1;
   }

   fread(&itmp, 4, 1, file);
   itmp = le_int(itmp);
   *rate = itmp;
   if (*rate != 8000 && *rate != 16000)
   {
      fprintf (stderr, "Only 8 kHz (narrowband) and 16 kHz (wideband) supported\n");
      return -1;
   }

   fread(&itmp, 4, 1, file);
   bpersec = le_int(itmp);

   fread(&stmp, 2, 1, file);
   balign = le_short(stmp);

   fread(&stmp, 2, 1, file);
   stmp = le_short(stmp);
   if (stmp!=16 && stmp!=8)
   {
      fprintf (stderr, "Only 8/16-bit linear supported\n");
      return -1;
   }
   *format=stmp;

   if (bpersec!=*rate**channels*stmp/8)
   {
      fprintf (stderr, "Corrupted header: ByteRate mismatch\n");
      return -1;
   }

   if (balign!=*channels*stmp/8)
   {
      fprintf (stderr, "Corrupted header: BlockAlign mismatch\n");
      return -1;
   }

   fread(ch, 1, 4, file);
   if (strcmp(ch, "data")!=0)
   {
      fprintf (stderr, "Corrupted WAVE file: no \"data\"\n");
      return -1;
   }

   /*Ignore this for now*/
   fread(&itmp, 4, 1, file);
   itmp = le_int(itmp);


   return 1;
}



void write_wav_header(FILE *file, int rate, int channels, int format, int size)
{
   char ch[5];
   int itmp;
   short stmp;

   ch[4]=0;

   fprintf (file, "RIFF");

   itmp = 0x7fffffff;
   fwrite(&itmp, 4, 1, file);

   fprintf (file, "WAVEfmt ");

   itmp = le_int(16);
   fwrite(&itmp, 4, 1, file);

   stmp = le_short(1);
   fwrite(&stmp, 2, 1, file);

   stmp = le_short(channels);
   fwrite(&stmp, 2, 1, file);

   itmp = le_int(rate);
   fwrite(&itmp, 4, 1, file);

   itmp = le_int(rate*channels*2);
   fwrite(&itmp, 4, 1, file);

   stmp = le_short(2*channels);
   fwrite(&stmp, 2, 1, file);

   stmp = le_short(16);
   fwrite(&stmp, 2, 1, file);

   fprintf (file, "data");

   itmp = le_int(0x7fffffff);
   fwrite(&itmp, 4, 1, file);


}
