/*---------------------------------------------------------------------------*\

	FILE........: AK2LSPD.H
	TYPE........: Turbo C header file
	COMPANY.....: Voicetronix
	AUTHOR......: James Whitehall
	DATE CREATED: 21/11/95

\*---------------------------------------------------------------------------*/
/* Speex License:

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

#ifndef __AK2LSPD__
#define __AK2LSPD__

int lpc_to_lsp (float *a, int lpcrdr, float *freq, int nb, float delta, float *stack);
void lsp_to_lpc(float *freq, float *ak, int lpcrdr, float *stack);

/*Added by JMV*/
void lsp_enforce_margin(float *lsp, int len, float margin);


#endif	/* __AK2LSPD__ */
