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
#include "wav_io.h"
#include "speex_header.h"
#include "misc.h"

/*Write an Ogg page to a file pointer*/
int oe_write_page(ogg_page *page, FILE *fp)
{
   int written;
   written = fwrite(page->header,1,page->header_len, fp);
   written += fwrite(page->body,1,page->body_len, fp);
   
   return written;
}

#define MAX_FRAME_SIZE 2000
#define MAX_FRAME_BYTES 2000

void usage()
{
   fprintf (stderr, "Speex encoder version " VERSION "\n");
   fprintf (stderr, "\n");
   fprintf (stderr, "usage: speexenc [options] input_file output_file\n");
   fprintf (stderr, "\n");
   fprintf (stderr, "input_file can be:\n");
   fprintf (stderr, "  filename.wav          wav file\n");
   fprintf (stderr, "  filename.*            raw PCM file (any extension other than .wav)\n");
   fprintf (stderr, "  -                     stdin\n");
   fprintf (stderr, "\n");  
   fprintf (stderr, "output_file can be:\n");
   fprintf (stderr, "  filename.spx          Speex file\n");
   fprintf (stderr, "  -                     stdout\n");
   fprintf (stderr, "\n");  
   fprintf (stderr, "options:\n");
   fprintf (stderr, "  --narrowband -n    Narrowband (8 kHz) input file\n"); 
   fprintf (stderr, "  --wideband   -w    Wideband (16 kHz) input file\n"); 
   fprintf (stderr, "  --quality n        Encoding quality setting from 0 to 10\n"); 
   fprintf (stderr, "  --lbr              Low bit-rate mode (equivalent to --quality 3)\n"); 
   fprintf (stderr, "  --vbr              Enable variable bit-rate (VBR)\n"); 
   fprintf (stderr, "  --comp n           Set encoding complexity (0-10)\n"); 
   fprintf (stderr, "  --nframes n        Number of frames per Ogg packet\n"); 
   fprintf (stderr, "  --help       -h    This help\n"); 
   fprintf (stderr, "  --version    -v    Version information\n"); 
   fprintf (stderr, "  -V                 Verbose mode (show bit-rate)\n"); 
   fprintf (stderr, "\n");  
   fprintf (stderr, "Input must be mono\n"); 
   fprintf (stderr, "Raw PCM needs to be 16-bit little-endian\n"); 
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
   int vbr_enabled=0;
   int i,nbBytes;
   SpeexMode *mode=NULL;
   void *st;
   SpeexBits bits;
   char cbits[MAX_FRAME_BYTES];
   struct option long_options[] =
   {
      {"wideband", no_argument, NULL, 0},
      {"narrowband", no_argument, NULL, 0},
      {"lbr", no_argument, NULL, 0},
      {"vbr", no_argument, NULL, 0},
      {"quality", required_argument, NULL, 0},
      {"nframes", required_argument, NULL, 0},
      {"comp", required_argument, NULL, 0},
      {"help", no_argument, NULL, 0},
      {"version", no_argument, NULL, 0},
      {0, 0, 0, 0}
   };
   int print_bitrate=0;
   int rate, chan, fmt, size;
   int quality=-1;
   int lbr=0;
   ogg_stream_state os;
   ogg_page 		 og;
   ogg_packet 		 op;
   int bytes_written, ret, result;
   int id=-1;
   SpeexHeader header;
   int nframes=1;
   int complexity=3;
   char *comments = "Encoded with Speex " VERSION;

   /*Process command-line options*/
   while(1)
   {
      c = getopt_long (argc, argv, "nwhvV",
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
         else if (strcmp(long_options[option_index].name,"lbr")==0)
               lbr=1;
         else if (strcmp(long_options[option_index].name,"vbr")==0)
               vbr_enabled=1;
         else if (strcmp(long_options[option_index].name,"quality")==0)
         {
            quality = atoi (optarg);
         } else if (strcmp(long_options[option_index].name,"nframes")==0)
         {
            nframes = atoi (optarg);
            if (nframes<1)
               nframes=1;
            if (nframes>10)
               nframes=10;
         } else if (strcmp(long_options[option_index].name,"comp")==0)
         {
            complexity = atoi (optarg);
         } else if (strcmp(long_options[option_index].name,"help")==0)
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
         exit(0);
         break;
      case 'v':
         version();
         exit(0);
         break;
      case 'V':
         print_bitrate=1;
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

   /*Initialize Ogg stream struct*/
   srand(time(NULL));
   if (ogg_stream_init(&os, rand())==-1)
   {
      fprintf(stderr,"Stream init failed\n");
      exit(1);
   }

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

   rate=0;
   if (strcmp(inFile+strlen(inFile)-4,".wav")==0)
      if (read_wav_header(fin, &rate, &chan, &fmt, &size)==-1)
         exit(1);
   /*fprintf (stderr, "wave info: %d %d %d %d\n", rate, chan, fmt, size);*/

   if (rate==16000)
   {
      wideband=1;
      if (narrowband)
         fprintf (stderr,"Warning: encoding a wideband file in narrowband\n");
   } else if (rate==8000)
   {
      narrowband=1;
      if (wideband)
         fprintf (stderr,"Warning: encoding a narrowband file in wideband\n");
   }

   if (!wideband)
      narrowband=1;
   if (narrowband)
   {
      if (!rate)
         rate = 8000;
      mode=&speex_nb_mode;
   }
   if (wideband)
   {
      if (!rate)
         rate = 16000;
      mode=&speex_wb_mode;
   }

   speex_init_header(&header, rate, 1, mode);
   header.frames_per_packet=nframes;
   header.vbr=vbr_enabled;

   fprintf (stderr, "Encoding %d Hz audio using %s mode\n", 
            header.rate, mode->modeName);
   /*fprintf (stderr, "Encoding %d Hz audio at %d bps using %s mode\n", 
     header.rate, mode->bitrate, mode->modeName);*/

   /*Initialize Speex encoder*/
   st = speex_encoder_init(mode);

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

      op.packet = (unsigned char *)speex_header_to_packet(&header, (int*)&(op.bytes));
      op.b_o_s = 1;
      op.e_o_s = 0;
      op.granulepos = 0;
      op.packetno = 0;
      ogg_stream_packetin(&os, &op);
      free(op.packet);

      op.packet = (unsigned char *)comments;
      op.bytes = strlen((char*)op.packet);
      op.b_o_s = 0;
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
   speex_encoder_ctl(st, SPEEX_SET_COMPLEXITY, &complexity);
   if (vbr_enabled)
   {
      int tmp;
      tmp=1;
      speex_encoder_ctl(st, SPEEX_SET_VBR, &tmp);
   }
   if (lbr || quality != -1)
   {
      int tmp=quality;
      if (quality==-1)
         tmp = 3;
      speex_encoder_ctl(st, SPEEX_SET_QUALITY, &tmp);
      if (vbr_enabled)
         speex_encoder_ctl(st, SPEEX_SET_VBR_QUALITY, &tmp);
   }

   /*Main encoding loop (one frame per iteration)*/
   while (1)
   {
      id++;
      /*Read input audio*/
      fread(in, sizeof(short), frame_size, fin);
      if (feof(fin))
         break;
      for (i=0;i<frame_size;i++)
         input[i]=(short)le_short(in[i]);
      /*Encode current frame*/
      speex_encode(st, input, &bits);
      
      if (print_bitrate) {
         int tmp;
         char ch=13;
         speex_encoder_ctl(st, SPEEX_GET_BITRATE, &tmp);
         fputc (ch, stderr);
         fprintf (stderr, "Bitrate is use: %d bps     ", tmp);
      }
      if ((id+1)%nframes!=0)
         continue;
      nbBytes = speex_bits_write(&bits, cbits, MAX_FRAME_BYTES);
      speex_bits_reset(&bits);
      op.packet = (unsigned char *)cbits;
      op.bytes = nbBytes;
      op.b_o_s = 0;
      op.e_o_s = 0;
      op.granulepos = id*frame_size;
      op.packetno = id/nframes;
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

