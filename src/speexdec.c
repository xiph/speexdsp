/* Copyright (C) 2002 Jean-Marc Valin 
   File: speexdec.c

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
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>

#include "modes.h"
#include "speex.h"

#define MAX_FRAME_SIZE 2000
#define MAX_FRAME_BYTES 1000

void usage()
{
   fprintf (stderr, "speexenc [options] <input file> <output file>\n");
   fprintf (stderr, "options:\n");
   fprintf (stderr, "\t--help       -h    This help\n"); 
   fprintf (stderr, "\t--version    -v    Version information\n"); 
}

void version()
{
   fprintf (stderr, "Speex encoder version " VERSION "\n");
}

int main(int argc, char **argv)
{
   int c;
   int option_index = 0;
   char *inFile, *outFile;
   FILE *fin, *fout;
   short out[MAX_FRAME_SIZE];
   float output[MAX_FRAME_SIZE];
   int frame_size;
   SpeexMode *mode=NULL;
   DecState st;
   FrameBits bits;
   char cbits[MAX_FRAME_BYTES];
   int at_end=0;
   struct option long_options[] =
   {
      {"help", no_argument, NULL, 0},
      {"version", no_argument, NULL, 0},
      {0, 0, 0, 0}
   };
   
   while(1)
   {
      c = getopt_long (argc, argv, "nwhv",
                       long_options, &option_index);
      if (c==-1)
         break;
      
      switch(c)
      {
      case 0:
         if (strcmp(long_options[option_index].name,"help")==0)
         {
            usage();
            exit(0);
         } else if (strcmp(long_options[option_index].name,"version")==0)
         {
            version();
            exit(0);
         }
         break;
      case 'h':
         usage();
         break;
      case 'v':
         version();
         exit(0);
         break;
      case '?':
         usage();
         exit(1);
         break;
      }
   }
   if (argc-optind!=2)
   {
      usage();
      exit(1);
   }
   inFile=argv[optind];
   outFile=argv[optind+1];

   if (strcmp(inFile, "-")==0)
      fin=stdin;
   else 
   {
      fin = fopen(inFile, "r");
      if (!fin)
      {
         perror(inFile);
         exit(1);
      }
   }
   if (strcmp(outFile,"-")==0)
      fout=stdout;
   else 
   {
      fout = fopen(outFile, "w");
      if (!fout)
      {
         perror(outFile);
         exit(1);
      }
   }

   /*This is only the temporary header*/
   {
      char header[6];
      if(fread(header, sizeof(char), 5, fin)!=5)
      {
         perror("cannot read header");
         exit(1);
      }
      header[5]=0;
      if (strcmp(header,"spexn")==0)
         mode=&mp_nb_mode;
      else if (strcmp(header,"spexw")==0)
         mode=&mp_wb_mode;
      else 
      {
         fprintf (stderr, "This does not look like a Speex " VERSION " file\n");
         exit(1);
      }
   }

   decoder_init(&st, mode);
   frame_size=mode->frameSize;
   frame_bits_init(&bits);
   while (1)
   {
      int i, nb_read;
      nb_read=200-((bits.nbBits>>3)-bits.bytePtr);
      if (nb_read>0&&!at_end)
      {
         nb_read=fread(cbits, sizeof(char), nb_read, fin);
         if (feof(fin))
            at_end=1;
         if (nb_read>0 && !at_end)
            frame_bits_read_whole_bytes(&bits, cbits, nb_read);
      }
      
      if (((bits.nbBits>>3)-bits.bytePtr)<2)
         break;

      decode(&st, &bits, output);
      for (i=0;i<frame_size;i++)
      {
         if (output[i]>32000)
            output[i]=32000;
         else if (output[i]<-32000)
            output[i]=-32000;
      }
      for (i=0;i<frame_size;i++)
         out[i]=output[i];
      fwrite(out, sizeof(short), frame_size, fout);
   }
   decoder_destroy(&st);
   exit(0);
   return 1;
}
