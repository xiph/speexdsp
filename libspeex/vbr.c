/* Copyright (C) 2002 Jean-Marc Valin 
   File: vbr.c

   VBR-related routines

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

#include "vbr.h"

/*
  This function should analyse the signal and decide how critical the
  coding error will be perceptually. The following factors should be
  taken into account:

  -Attacks (positive energy derivative) should be coded with more bits

  -Stationary voiced segments should receive more bits

  -Segments with (very) low absolute energy should receive less bits (maybe
  only shaped noise?)

  -DTX for near-zero energy?

  -Stationary fricative segments should have less bits

  -Temporal masking: when energy slope is decreasing, decrease the bit-rate

  -Decrease bit-rate for males (low pitch)?

  -(wideband only) less bits in the high-band when signal is very 
  non-stationary (harder to notice high-frequency noise)???

*/
void vbr_analysis(float *sig, int len)
{
   
}
