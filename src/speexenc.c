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

#include "speex.h"
#include <ogg/ogg.h>

/*Write an Ogg page to a file pointer*/
int oe_write_page(ogg_page *page, FILE *fp)
{
   int written;
   written = fwrite(page->header,1,page->header_len, fp);
   written += fwrite(page->body,1,page->body_len, fp);
   
   return written;
}

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
   void *st;
   SpeexBits bits;
   char cbits[MAX_FRAME_BYTES];
   struct option long_options[] =
   {
      {"wideband", no_argument, NULL, 0},
      {"narrowband", no_argument, NULL, 0},
      {"help", no_argument, NULL, 0},
      {"version", no_argument, NULL, 0},
      {0, 0, 0, 0}
   };

   ogg_stream_state os;
   ogg_page 		 og;
   ogg_packet 		 op;
   int bytes_written, ret, result;
   int id=0;

   /*Process command-line options*/
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

   /*Initialize Ogg stream struct*/
   if (ogg_stream_init(&os, 0)==-1)
   {
      fprintf(stderr,"Stream init failed\n");
      exit(1);
   }

   if (wideband && narrowband)
   {
      fprintf (stderr,"Cannot specify both wideband and narrowband at the same time\n");
      exit(1);
   };
   if (!wideband)
      narrowband=1;
   if (narrowband)
      mode=&speex_nb_mode;
   if (wideband)
      mode=&speex_wb_mode;

   /*Initialize Speex encoder*/
   st = speex_encoder_init(mode);

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

   /*Write header (format will change)*/
   {

      if (narrowband)
         op.packet = (unsigned char *)"speex narrowband";
      if (wideband)
         op.packet = (unsigned char *)"speex wideband**";
      op.bytes = 16;
      op.b_o_s = 1;
      op.e_o_s = 0;
      op.granulepos = 0;
      op.packetno = 0;
      ogg_stream_packetin(&os, &op);

      while((result = ogg_stream_flush(&os, &og)))
      {
         if(!result) break;
         ret = oe_write_page(&og, fout);
         if(ret != og.header_len + og.body_len)
         {
            fprintf (stderr,"Failed writing header to output stream\n");
            exit(1);
         }
         else
            bytes_written += ret;
      }
   }

   speex_encoder_ctl(st, SPEEX_GET_FRAME_SIZE, &frame_size);

   /*Main encoding loop (one frame per iteration)*/
   while (1)
   {
      id++;
      /*Read input audio*/
      fread(in, sizeof(short), frame_size, fin);
      if (feof(fin))
         break;
      for (i=0;i<frame_size;i++)
         input[i]=in[i];
      /*Encode current frame*/
      speex_encode(st, input, &bits);

      /*if (id%5!=0)
        continue;*/
      nbBytes = speex_bits_write(&bits, cbits, 500);
      speex_bits_reset(&bits);
      op.packet = (unsigned char *)cbits;
      op.bytes = nbBytes;
      op.b_o_s = 0;
      op.e_o_s = 0;
      op.granulepos = id;
      op.packetno = id;
      ogg_stream_packetin(&os, &op);

      /*Write all new pages (not likely 0 or 1)*/
      while (ogg_stream_pageout(&os,&og))
      {
         ret = oe_write_page(&og, fout);
         if(ret != og.header_len + og.body_len)
         {
            fprintf (stderr,"Failed writing header to output stream\n");
            exit(1);
         }
         else
            bytes_written += ret;
      }
   }
   
   op.packet = (unsigned char *)"END OF STREAM";
   op.bytes = 13;
   op.b_o_s = 0;
   op.e_o_s = 1;
   op.granulepos = id+1;
   op.packetno = id+1;
   ogg_stream_packetin(&os, &op);
   /*Flush all pages left to be written*/
   while (ogg_stream_flush(&os, &og))
   {
      ret = oe_write_page(&og, fout);
      if(ret != og.header_len + og.body_len)
      {
         fprintf (stderr,"Failed writing header to output stream\n");
         exit(1);
      }
      else
         bytes_written += ret;
   }
   

   speex_encoder_destroy(st);
   speex_bits_destroy(&bits);
   ogg_stream_clear(&os);

   exit(0);
   return 1;
}

