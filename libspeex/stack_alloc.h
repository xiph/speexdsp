/* Copyright (C) 2002 Jean-Marc Valin 
   File: stack_alloc.h
   
   Temporary memory allocation on stack

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

#ifndef STACK_ALLOC_H
#define STACK_ALLOC_H


#ifdef __GNUC__
#define VAR_ARRAY
#endif


#if defined (VAR_ARRAY)  /* Prefered method is variable-size arrays is supported */

#define DYN_VEC(type, num, var) type var[num];

#elif defined (HAVE_ALLOCA_H)  /* Second best: alloca */

#include <alloca.h>

#define DYN_VEC(type, num, var) type *var=(type*)alloca((num)*sizeof(type));

#elif defined WIN32  /* On Win32, it's _alloca */

#include <malloc.h>
#define DYN_VEC(type, num, var) type *var=(type*)_alloca((num)*sizeof(type));

#else  /* When all else fails, allocate on the heap (but it's going to be slow) */

#error Cannot allocate memory on stack

#endif



#endif
