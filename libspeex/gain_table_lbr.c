/* Copyright (C) 2002 Jean-Marc Valin 
   File: gain_table_lbr.c
   Codebook for 3-tap pitch prediction gain (32 entries)
  
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:

   1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.  

   2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

   3. The name of the author may not be used to endorse or promote products
   derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
   IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
   OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
   INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
   HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
   STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

float gain_cdbk_lbr[] = {
0,0,0,
0.019578,  -0.411369,  0.250244,
-0.141413,  0.127455,  -0.177373,
-0.379174,  0.154715,  -0.359933,
0.295340,  1.014952,  -0.144606,
0.431555,  -0.107415,  0.360701,
-0.141305,  0.735394,  0.312635,
0.382416,  0.267769,  0.318738,
0.511146,  0.524061,  -0.190435,
0.153482,  -0.531485,  -0.149959,
-0.094091,  0.930054,  0.139366,
0.164167,  0.711936,  -0.077780,
0.503705,  0.823130,  -0.273699,
-0.330264,  -0.613346,  0.085310,
-0.083597,  0.481953,  0.201470,
0.195682,  0.429066,  0.059682,
0.598746,  1.523378,  -0.189717,
-0.010502,  -0.257728,  -0.018047,
-0.132438,  1.383543,  0.280042,
0.234771,  0.555249,  -0.210053,
0.010973,  1.090455,  -0.009557,
0.141315,  0.930896,  -0.128939,
-0.168645,  0.950529,  0.314244,
-0.028768,  0.695554,  0.133637,
0.246305,  0.740436,  0.073124,
0.280190,  -0.787092,  0.268726,
0.010162,  0.894487,  0.006648,
0.177218,  0.572144,  0.427882,
-0.237882,  -0.484537,  -0.303846,
-0.211570,  0.684685,  0.539195,
0.064373,  0.236576,  0.042304,
0.347794,  0.726175,  -0.126887,
};
