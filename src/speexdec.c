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

#include "speex.h"
#include "ogg/ogg.h"

#ifdef HAVE_SYS_SOUNDCARD_H
#include <sys/soundcard.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#endif

#include <string.h>
#include "wav_io.h"

#define MAX_FRAME_SIZE 2000

FILE *out_file_open(char *outFile, int rate)
{
   FILE *fout;
   /*Open output file*/
   if (strlen(outFile)==0)
   {
#ifdef HAVE_SYS_SOUNDCARD_H
      int audio_fd, format, stereo;
      audio_fd=open("/dev/dsp", O_WRONLY);
      
      format=AFMT_S16_LE;
      if (ioctl(audio_fd, SNDCTL_DSP_SETFMT, &format)==-1)
      {
         perror("SNDCTL_DSP_SETFMT");
         close(audio_fd);
         exit(1);
      }
      
      stereo=0;
      if (ioctl(audio_fd, SNDCTL_DSP_STEREO, &stereo)==-1)
      {
         perror("SNDCTL_DSP_STEREO");
         close(audio_fd);
         exit(1);
      }
      if (stereo!=0)
      {
         fprintf (stderr, "Cannot set mono mode\n");
         exit(1);
      }

      if (ioctl(audio_fd, SNDCTL_DSP_SPEED, &rate)==-1)
      {
         perror("SNDCTL_DSP_SPEED");
         close(audio_fd);
         exit(1);
      }
      fout = fdopen(audio_fd, "w");
#else
      fprintf (stderr, "No soundcard support\n");
      exit(1);
#endif
   } else {
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
         if (strcmp(outFile+strlen(outFile)-4,".wav")==0)
            write_wav_header(fout, rate, 1, 0, 0);
      }
   }
   return fout;
}

void usage()
{
   fprintf (stderr, "speexenc [options] <input file> <output file>\n");
   fprintf (stderr, "options:\n");
   fprintf (stderr, "\t--help       -h      This help\n"); 
   fprintf (stderr, "\t--version    -v      Version information\n"); 
   fprintf (stderr, "\t--pf         --pf    Enable post-filter\n"); 
   fprintf (stderr, "\t--no-pf      --no-pf Disable post-filter\n"); 
}

void version()
{
   fprintf (stderr, "Speex decoder version " VERSION "\n");
}

int main(int argc, char **argv)
{
   int c;
   int option_index = 0;
   char *inFile, *outFile;
   FILE *fin, *fout=NULL;
   short out[MAX_FRAME_SIZE];
   float output[MAX_FRAME_SIZE];
   int frame_size=0;
   SpeexMode *mode=NULL;
   void *st=NULL;
   SpeexBits bits;
   int first = 1;
   struct option long_options[] =
   {
      {"help", no_argument, NULL, 0},
      {"version", no_argument, NULL, 0},
      {"pf", no_argument, NULL, 0},
      {"no-pf", no_argument, NULL, 0},
      {0, 0, 0, 0}
   };
   ogg_sync_state oy;
   ogg_page       og;
   ogg_packet     op;
   ogg_stream_state os;
   int pf_enabled;

   pf_enabled = 0;

   /*Process options*/
   while(1)
   {
      c = getopt_long (argc, argv, "hv",
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
         } else if (strcmp(long_options[option_index].name,"pf")==0)
         {
            pf_enabled=1;
         } else if (strcmp(long_options[option_index].name,"no-pf")==0)
         {
            pf_enabled=0;
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
   if (argc-optind!=2 && argc-optind!=1)
   {
      usage();
      exit(1);
   }
   inFile=argv[optind];

   if (argc-optind==2)
      outFile=argv[optind+1];
   else
      outFile = "";
   /*Open input file*/
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


   /*Init Ogg data struct*/
   ogg_sync_init(&oy);
   ogg_stream_init(&os, 0);
   
   speex_bits_init(&bits);
   /*Main decoding loop*/
   while (1)
   {
      char *data;
      int i, nb_read;
      /*Get the ogg buffer for writing*/
      data = ogg_sync_buffer(&oy, 200);
      /*Read bitstream from input file*/
      nb_read = fread(data, sizeof(char), 200, fin);      
      ogg_sync_wrote(&oy, nb_read);

      /*Loop for all complete pages we got (most likely only one)*/
      while (ogg_sync_pageout(&oy, &og)==1)
      {
         /*Add page to the bitstream*/
         ogg_stream_pagein(&os, &og);
         /*Extract all available packets*/
         while (ogg_stream_packetout(&os, &op)==1)
         {
            /*If first packet, process as Speex header*/
            if (first)
            {
               int rate;
               if (strncmp((char *)op.packet, "speex wideband**", 12)==0)
               {
                  rate=16000;
                  mode = &speex_wb_mode;
               } else if (strncmp((char *)op.packet, "speex narrowband", 12)==0)
               {
                  rate=8000;
                  mode = &speex_nb_mode;
               } else if (strncmp((char *)op.packet, "speex narrow-lbr", 12)==0)
               {
                  rate=8000;
                  mode = &speex_nb_lbr_mode;
               } else {
                  fprintf (stderr, "This Ogg file is not a Speex bitstream\n");
                  exit(1);
               }
               /*Initialize Speex decoder*/
               st = speex_decoder_init(mode);
               speex_decoder_ctl(st, SPEEX_SET_PF, &pf_enabled);
               speex_decoder_ctl(st, SPEEX_GET_FRAME_SIZE, &frame_size);

               fout = out_file_open(outFile, rate);

               first=0;
            } else {
               /*End of stream condition*/
               if (strncmp((char *)op.packet, "END OF STREAM", 13)==0)
                  break;
               /*Copy Ogg packet to Speex bitstream*/
               speex_bits_read_from(&bits, (char*)op.packet, op.bytes);
               /*Decode a frame*/
               speex_decode(st, &bits, output, 0);
               
               /*PCM saturation (just in case)*/
               for (i=0;i<frame_size;i++)
               {
                  if (output[i]>32000)
                     output[i]=32000;
                  else if (output[i]<-32000)
                     output[i]=-32000;
               }
               /*Convert to short and save to output file*/
               for (i=0;i<frame_size;i++)
                  out[i]=output[i];
               fwrite(out, sizeof(short), frame_size, fout);
            }
         }
      }
      if (feof(fin))
         break;

   }

   speex_decoder_destroy(st);
   speex_bits_destroy(&bits);
   ogg_stream_clear(&os);
 
   exit(0);
   return 1;
}
