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
#include "speex_header.h"
#include "misc.h"

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

   fprintf (stderr, "Speex decoder version " VERSION "\n");
   fprintf (stderr, "\n");
   fprintf (stderr, "usage: speexdec [options] input_file.spx\n");
   fprintf (stderr, "       speexdec [options] input_file.spx output_file.wav\n");
   fprintf (stderr, "\n");
   fprintf (stderr, "input_file can be:\n");
   fprintf (stderr, "  filename.spx          regular Speex file\n");
   fprintf (stderr, "  -                     stdin\n");
   fprintf (stderr, "\n");  
   fprintf (stderr, "output_file can be:\n");
   fprintf (stderr, "  filename.wav          wav file\n");
   fprintf (stderr, "  filename.*            raw PCM file (any extension other that .wav)\n");
   fprintf (stderr, "  -                     stdout\n");
   fprintf (stderr, "  (nothing)             will be played to soundcard\n");
   fprintf (stderr, "\n");  
   fprintf (stderr, "options:\n");
   fprintf (stderr, "  --enh                 Enable perceptual enhancement\n");
   fprintf (stderr, "  --no-enh              Disable perceptual enhancement (default FOR NOW)\n");
   fprintf (stderr, "  -V                    Verbose mode (show bit-rate)\n"); 
   fprintf (stderr, "  --help       -h       This help\n");
   fprintf (stderr, "  --version    -v       Version information\n");
   fprintf (stderr, "  --pf                  Deprecated, use --pf instead\n");
   fprintf (stderr, "  --no-pf               Deprecated, use --no-pf instead\n");
}

void version()
{
   fprintf (stderr, "Speex decoder version " VERSION "\n");
}

static void *process_header(ogg_packet *op, int enh_enabled, int *frame_size, int *rate, int *nframes)
{
   void *st;
   SpeexMode *mode;
   SpeexHeader *header;
   
   header = speex_packet_to_header((char*)op->packet, op->bytes);
   if (!header)
   {
      fprintf (stderr, "Cannot read header\n");
      return NULL;
   }
   if (header->mode >= SPEEX_NB_MODES)
   {
      fprintf (stderr, "Mode number %d does not (any longer) exist in this version\n", 
               header->mode);
      return NULL;
   }
      
   mode = speex_mode_list[header->mode];
   
   if (mode->bitstream_version < header->mode_bitstream_version)
   {
      fprintf (stderr, "The file was encoded with a newer version of Speex. You need to upgrade in order to play it.\n");
      return NULL;
   }
   if (mode->bitstream_version > header->mode_bitstream_version) 
   {
      fprintf (stderr, "The file was encoded with an older version of Speex. You would need to downgrade the version in order to play it.\n");
      return NULL;
   }
   
   st = speex_decoder_init(mode);
   speex_decoder_ctl(st, SPEEX_SET_ENH, &enh_enabled);
   speex_decoder_ctl(st, SPEEX_GET_FRAME_SIZE, frame_size);
   
   *rate = header->rate;
   *nframes = header->frames_per_packet;
   
   fprintf (stderr, "Decoding %d Hz audio using %s mode", 
            *rate, mode->modeName);

   if (header->vbr)
      fprintf (stderr, " (VBR)\n");
   else
      fprintf(stderr, "\n");
   /*fprintf (stderr, "Decoding %d Hz audio at %d bps using %s mode\n", 
    *rate, mode->bitrate, mode->modeName);*/

   free(header);
   return st;
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
   void *st=NULL;
   SpeexBits bits;
   int packet_count=0;
   int stream_init = 0;
   struct option long_options[] =
   {
      {"help", no_argument, NULL, 0},
      {"version", no_argument, NULL, 0},
      {"enh", no_argument, NULL, 0},
      {"no-enh", no_argument, NULL, 0},
      {"pf", no_argument, NULL, 0},
      {"no-pf", no_argument, NULL, 0},
      {0, 0, 0, 0}
   };
   ogg_sync_state oy;
   ogg_page       og;
   ogg_packet     op;
   ogg_stream_state os;
   int enh_enabled;
   int nframes=2;
   int print_bitrate=0;
   int close_in=0;
   int eos=0;

   enh_enabled = 0;

   /*Process options*/
   while(1)
   {
      c = getopt_long (argc, argv, "hvV",
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
         } else if (strcmp(long_options[option_index].name,"enh")==0)
         {
            enh_enabled=1;
         } else if (strcmp(long_options[option_index].name,"no-enh")==0)
         {
            enh_enabled=0;
         } else if (strcmp(long_options[option_index].name,"pf")==0)
         {
            fprintf (stderr, "--pf is deprecated, use --enh instead\n");
            enh_enabled=1;
         } else if (strcmp(long_options[option_index].name,"no-pf")==0)
         {
            fprintf (stderr, "--no-pf is deprecated, use --no-enh instead\n");
            enh_enabled=0;
         }
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
      close_in=1;
   }


   /*Init Ogg data struct*/
   ogg_sync_init(&oy);
   
   speex_bits_init(&bits);
   /*Main decoding loop*/
   while (1)
   {
      char *data;
      int i, j, nb_read;
      /*Get the ogg buffer for writing*/
      data = ogg_sync_buffer(&oy, 200);
      /*Read bitstream from input file*/
      nb_read = fread(data, sizeof(char), 200, fin);      
      ogg_sync_wrote(&oy, nb_read);

      /*Loop for all complete pages we got (most likely only one)*/
      while (ogg_sync_pageout(&oy, &og)==1)
      {
         if (stream_init == 0) {
            ogg_stream_init(&os, ogg_page_serialno(&og));
            stream_init = 1;
         }
         /*Add page to the bitstream*/
         ogg_stream_pagein(&os, &og);
         /*Extract all available packets*/
         while (!eos && ogg_stream_packetout(&os, &op)==1)
         {
            /*If first packet, process as Speex header*/
            if (packet_count==0)
            {
               int rate;
               st = process_header(&op, enh_enabled, &frame_size, &rate, &nframes);
               if (!nframes)
                  nframes=1;
               if (!st)
                  exit(1);
               fout = out_file_open(outFile, rate);

            } else if (packet_count==1){
               fprintf (stderr, "File comments: ");
               fwrite(op.packet, 1, op.bytes, stderr);
               fprintf (stderr, "\n");
            } else {

               /*End of stream condition*/
               if (op.e_o_s)
                  eos=1;

               /*Copy Ogg packet to Speex bitstream*/
               speex_bits_read_from(&bits, (char*)op.packet, op.bytes);
               for (j=0;j<nframes;j++)
               {
                  /*Decode frame*/
                  speex_decode(st, &bits, output, 0);
               
                  if (print_bitrate) {
                     int tmp;
                     char ch=13;
                     speex_decoder_ctl(st, SPEEX_GET_BITRATE, &tmp);
                     fputc (ch, stderr);
                     fprintf (stderr, "Bitrate is use: %d bps     ", tmp);
                  }
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
                     out[i]=(short)le_short(output[i]);
                  fwrite(out, sizeof(short), frame_size, fout);
               }
            }
            packet_count++;
         }
      }
      if (feof(fin))
         break;

   }

   if (st)
      speex_decoder_destroy(st);
   speex_bits_destroy(&bits);
   ogg_sync_clear(&oy);
   ogg_stream_clear(&os);

   if (close_in)
      fclose(fin);
   fclose(fout);
   return 1;
}
