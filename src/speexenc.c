/* Copyright (C) 2002 Jean-Marc Valin 
   File: speexenc.c

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
   fprintf (stderr, "\t--narrowband -n    Narrowband (8 kHz) input file\n"); 
   fprintf (stderr, "\t--wideband   -w    Wideband (16 kHz) input file\n"); 
   fprintf (stderr, "\t--help       -h    This help\n"); 
   fprintf (stderr, "\t--version    -v    Version information\n"); 
   fprintf (stderr, "\nInput must be raw audio (no header), 16 bits\n"); 
}

void version()
{
   fprintf (stderr, "Speex encoder version " VERSION "\n");
}

int main(int argc, char **argv)
{
   int c;
   int option_index = 0;
   int narrowband=0, wideband=0;
   char *inFile, *outFile;
   FILE *fin, *fout;
   short in[MAX_FRAME_SIZE];
   float input[MAX_FRAME_SIZE];
   int frame_size;
   int i,nbBytes;
   SpeexMode *mode=NULL;
   EncState st;
   FrameBits bits;
   char cbits[MAX_FRAME_BYTES];
   struct option long_options[] =
   {
      {"wideband", no_argument, NULL, 0},
      {"narrowband", no_argument, NULL, 0},
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
         if (strcmp(long_options[option_index].name,"narrowband")==0)
            narrowband=1;
         else if (strcmp(long_options[option_index].name,"wideband")==0)
               wideband=1;
         else if (strcmp(long_options[option_index].name,"help")==0)
         {
            usage();
            exit(0);
         } else if (strcmp(long_options[option_index].name,"version")==0)
         {
            version();
            exit(0);
         }
         break;
      case 'n':
         narrowband=1;
         break;
      case 'h':
         usage();
         break;
      case 'v':
         version();
         exit(0);
         break;
      case 'w':
         wideband=1;
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

   if (wideband && narrowband)
   {
      fprintf (stderr,"Cannot specify both wideband and narrowband at the same time\n");
      exit(1);
   };
   if (!wideband)
      narrowband=1;
   if (narrowband)
      mode=&mp_nb_mode;
   if (wideband)
      mode=&mp_wb_mode;

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

   /*Temporary header*/
   {
      char *header;
      if (narrowband)
         header="spexn";
      if (wideband)
         header="spexw";
      if (fwrite(header, 1, 5, fout)!=5)
      {
         perror("Cannot write header");
      }
   }

   encoder_init(&st, mode);
   frame_size=mode->frameSize;

   while (1)
   {
      fread(in, sizeof(short), frame_size, fin);
      if (feof(fin))
         break;
      for (i=0;i<frame_size;i++)
         input[i]=in[i];
      encode(&st, input, &bits);
      nbBytes = frame_bits_write_whole_bytes(&bits, cbits, 200);
      fwrite(cbits, 1, nbBytes, fout);
   }
   nbBytes = frame_bits_write(&bits, cbits, 200);
   fwrite(cbits, 1, nbBytes, fout);
   encoder_destroy(&st);
   exit(0);
   return 1;
}
