/* Copyright (C) 2002 Jean-Marc Valin 
   File: modes.c

   Describes the different modes of the codec

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

#include <stdlib.h>
#include "modes.h"
#include "ltp.h"
#include "quant_lsp.h"
#include "cb_search.h"
#include "mpulse.h"
#include "sb_celp.h"
#include "nb_celp.h"

extern float gain_cdbk_nb[];
extern float exc_gains_table[];
extern float exc_table[];
extern float exc_wb_table[];
extern float exc_gains_wb_table[];
extern float exc_sb_table[];
extern float hexc_table[];

static ltp_params ltp_params_nb = {
   gain_cdbk_nb,
   7,
   7
};

static ltp_params ltp_params_wb = {
   gain_cdbk_nb,
   7,
   8
};

static split_cb_params split_cb_nb = {
   8,               /*subvect_size*/
   5,               /*nb_subvect*/
   exc_table,       /*shape_cb*/
   7,               /*shape_bits*/
};

static split_cb_params split_cb_sb = {
   5,               /*subvect_size*/
   8,              /*nb_subvect*/
   exc_sb_table,    /*shape_cb*/
   7,               /*shape_bits*/
};

static split_cb_params split_cb_wb = {
   8,               /*subvect_size*/
   10,              /*nb_subvect*/
   exc_wb_table,    /*shape_cb*/
   7,               /*shape_bits*/
};

static split_cb_params split_cb_high = {
   8,               /*subvect_size*/
   5,               /*nb_subvect*/
   hexc_table,       /*shape_cb*/
   8,               /*shape_bits*/
};


static mpulse_params mpulse_nb = {
   15,     /*nb_pulse*/
   5,      /*nb_tracks*/
   4,    /*gain_coef*/
   10
};


static mpulse_params mpulse_sb = {
   15,     /*nb_pulse*/
   5,      /*nb_tracks*/
   3.8,    /*gain_coef*/
   10
};


static mpulse_params mpulse_wb = {
   24,     /*nb_pulse*/
   4,      /*nb_tracks*/
   2.2,    /*gain_coef*/
   26,
};


static SpeexNBMode mp_nb_mode = {
   160,    /*frameSize*/
   40,     /*subframeSize*/
   320,    /*windowSize*/
   10,     /*lpcSize*/
   640,    /*bufSize*/
   17,     /*pitchStart*/
   144,    /*pitchEnd*/
   0.9,    /*gamma1*/
   0.6,    /*gamma2*/
   .005,   /*lag_factor*/
   1.0001, /*lpc_floor*/
   0.0,    /*preemph*/
   /*LSP quantization*/
   lsp_quant_nb,
   lsp_unquant_nb,
   /*Pitch quantization*/
   pitch_search_3tap,
   pitch_unquant_3tap,
   &ltp_params_nb,
   /*Innovation quantization*/
   mpulse_search,
   mpulse_unquant,
   &mpulse_nb
};

static SpeexNBMode low_sb_mode = {
   160,    /*frameSize*/
   40,     /*subframeSize*/
   320,    /*windowSize*/
   10,     /*lpcSize*/
   640,    /*bufSize*/
   17,     /*pitchStart*/
   144,    /*pitchEnd*/
   .9,    /*gamma1*/
   0.6,    /*gamma2*/
   .002,   /*lag_factor*/
   1.00005, /*lpc_floor*/
   0.0,    /*preemph*/
   /*LSP quantization*/
   lsp_quant_nb,
   lsp_unquant_nb,
   /*Pitch quantization*/
   pitch_search_3tap,
   pitch_unquant_3tap,
   &ltp_params_nb,
   /*Innovation quantization*/
#if 1
   split_cb_search2,
   split_cb_unquant,
   &split_cb_sb
#else
   mpulse_search,
   mpulse_unquant,
   &mpulse_sb
#endif
};

SpeexMode low_wb_mode = {
   &low_sb_mode,
   &nb_encoder_init,
   &nb_encoder_destroy,
   &nb_encode,
   &nb_decoder_init,
   &nb_decoder_destroy,
   &nb_decode,
   160
};

static SpeexSBMode sb_wb_mode = {
   &low_wb_mode,
   160,    /*frameSize*/
   40,     /*subframeSize*/
   320,    /*windowSize*/
   8,     /*lpcSize*/
   640,    /*bufSize*/
   .9,    /*gamma1*/
   0.6,    /*gamma2*/
   .002,   /*lag_factor*/
   1.0001, /*lpc_floor*/
   0.0,    /*preemph*/
   /*LSP quantization*/
   lsp_quant_high,
   lsp_unquant_high,
   /*Innovation quantization*/
   split_cb_search_nogain,
   split_cb_nogain_unquant,
   &split_cb_high
};


SpeexMode nb_mode = {
   &mp_nb_mode,
   &nb_encoder_init,
   &nb_encoder_destroy,
   &nb_encode,
   &nb_decoder_init,
   &nb_decoder_destroy,
   &nb_decode,
   160
};

SpeexMode wb_mode = {
   &sb_wb_mode,
   &sb_encoder_init,
   &sb_encoder_destroy,
   &sb_encode,
   &sb_decoder_init,
   &sb_decoder_destroy,
   &sb_decode,
   320
};




void *encoder_init(SpeexMode *mode)
{
   return mode->enc_init(mode);
}

void *decoder_init(SpeexMode *mode)
{
   return mode->dec_init(mode);
}

void encoder_destroy(void *state)
{
   (*((SpeexMode**)state))->enc_destroy(state);
}

void encode(void *state, float *in, FrameBits *bits)
{
   (*((SpeexMode**)state))->enc(state, in, bits);
}

void decoder_destroy(void *state)
{
   (*((SpeexMode**)state))->dec_destroy(state);
}

void decode(void *state, FrameBits *bits, float *out)
{
   (*((SpeexMode**)state))->dec(state, bits, out);
}
