/* Copyright (C) 2002 Jean-Marc Valin 
   File: matrix.h

   Matrix stuff

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


void solve(float *A, float *b, float *x, int N)
{
   int i,j,k;
   for (i=0;i<N;i++)
      x[i]=b[i];
   
   for (i=0;i<N;i++)
   {
      float d,d_1;
      d=A[i*N+i];
      d_1=1/d;
      for (j=i+1;j<N;j++)
      {
         float fact=A[j*N+i]*d_1;
         for (k=0;k<N;k++)
         {
            A[j*N+k]-=fact*A[i*N+k];
         }
         x[j]-=fact*x[i];
      }
   }
   
   
   for (i=N-1;i>=0;i--)
   {
      float d=A[i*N+i];
      for (j=0;j<i;j++)
      {
         x[j]-=A[j*N+i]/d*x[i];
      }
   }
   
   for (i=0;i<N;i++)
   {
      x[i]/=A[i*N+i];
   }
   
  
}
