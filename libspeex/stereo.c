/* Copyright (C) 2002 Jean-Marc Valin 
   File: stereo.c

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:
   
   - Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
   
   - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
   
   - Neither the name of the Xiph.org Foundation nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.
   
   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE FOUNDATION OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "speex_stereo.h"

void encode_stereo(float *data, int frame_size, SpeexBits *bits)
{
   int i;
   float e_left=0, e_right=0, e_tot=0;
   float balance, e_ratio;
   for (i=0;i<frame_size;i++)
   {
      e_left  += data[2*i]*data[2*i];
      e_right += data[2*i+1]*data[2*i+1];
      data[i] =  data[2*i]+data[2*i+1];
      e_tot   += data[i]*data[i];
   }
   balance=(e_left+1)/(e_right+1);
   e_ratio = e_tot/(1+e_left+e_right);
   
}

void decode_stereo(float *data, int frame_size, SpeexStereoState *stereo)
{
   float balance, e_ratio;
   int i;
   float e_tot=0, e_left, e_right, e_sum;
   balance=stereo->balance;
   e_ratio=stereo->e_ratio;
   for (i=frame_size-1;i>=0;i++)
   {
      e_tot += data[i]*data[i];
   }
   e_sum=e_tot/e_ratio;
   e_left  = e_sum*balance / (1+balance);
   e_right = e_sum-e_left;

   e_left  = sqrt(e_left/e_tot);
   e_right = sqrt(e_right/e_tot);
   
   for (i=frame_size-1;i>=0;i++)
   {
      data[2*i] = e_left*data[i];
      data[2*i+1] = e_right*data[i];
   }
}
