/* Copyright (C) 2002 Jean-Marc Valin 
   File: mics.h
   Various utility routines for Speex

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

#ifndef MISC_H
#define MISC_H

#ifndef VERSION
#define VERSION "unknown version"
#endif

unsigned int be_int(unsigned int i);
unsigned int le_int(unsigned int i);


unsigned short be_short(unsigned short s);
unsigned short le_short(unsigned short s);

void *speex_alloc (int size);
void speex_free (void *ptr);
void *speex_move (void *dest, void *src, int n);

#endif
