/* Copyright (C) 2002 Jean-Marc Valin 
   File: wav_io.c

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

int read_wav_header(FILE *file, int *rate, int *channels, int *format, int *size)
{
   char ch[5];
   int itmp;
   short stmp;

   ch[4]=0;
   fread(ch, 1, 4, file);
   if (strcmp(ch, "RIFF")!=0)
   {
      fseek(file, 0, SEEK_SET);
      return 0;
   }

   fread(&itmp, 4, 1, file);
   /*FIXME: swap bytes*/
   *size = itmp-36;

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
   /*FIXME: swap bytes*/
   if (itmp!=16)
   {
      fprintf (stderr, "Only 16-bit PCM supported\n");
      return -1;
   }

   fread(&stmp, 2, 1, file);
   /*FIXME: swap bytes*/
   if (stmp!=1)
   {
      fprintf (stderr, "Only 16-bit PCM supported\n");
      return -1;
   }

   fread(&stmp, 2, 1, file);
   /*FIXME: swap bytes*/
   *channels = stmp;

   fread(&itmp, 4, 1, file);
   /*FIXME: swap bytes*/
   *rate = itmp;

   fread(&itmp, 4, 1, file);
   /*FIXME: swap bytes*/
   if (itmp!=*rate**channels*2)
   {
      fprintf (stderr, "Corrupted header: ByteRate mismatch\n");
      return -1;
   }

   fread(&stmp, 2, 1, file);
   /*FIXME: swap bytes*/
   if (stmp!=*channels*2)
   {
      fprintf (stderr, "Corrupted header: BlockAlign mismatch\n");
      return -1;
   }

   fread(&stmp, 2, 1, file);
   /*FIXME: swap bytes*/
   if (stmp!=16)
   {
      fprintf (stderr, "Only 16-bit linear supported\n");
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
   /*FIXME: swap bytes*/

   *format=16;

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

   itmp = 16;
   fwrite(&itmp, 4, 1, file);

   stmp = 1;
   fwrite(&stmp, 2, 1, file);

   stmp = channels;
   fwrite(&stmp, 2, 1, file);

   itmp = rate;
   fwrite(&itmp, 4, 1, file);

   itmp = rate*channels*2;
   fwrite(&itmp, 4, 1, file);

   stmp = 2*channels;
   fwrite(&stmp, 2, 1, file);

   stmp = 16;
   fwrite(&stmp, 2, 1, file);

   fprintf (file, "data");

   itmp = 0x7fffffff;
   fwrite(&itmp, 4, 1, file);


}
