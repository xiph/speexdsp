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
#include "sb_celp.h"
#include "nb_celp.h"
#include "post_filter.h"


SpeexMode *speex_mode_list[] = {&speex_nb_mode, &speex_nb_lbr_mode, &speex_wb_mode};

/* Extern declarations for all codebooks we use here */
extern float gain_cdbk_nb[];
extern float gain_cdbk_lbr[];
extern float hexc_table[];
extern float exc_5_256_table[];
extern float exc_5_64_table[];
extern float exc_10_32_table[];

/* Post-filter parameters for narrowband */
pf_params pf_params_nb = {
   0.65,      /* formant enhancement numerator */
   0.7,      /* formant enhancement denominator */
   0.4       /* pitch enhancement factor */
};

/* Post-filter parameters for low bit-rate narrowband */
pf_params pf_params_lbr = {
   0.65,      /* formant enhancement numerator */
   0.7,      /* formant enhancement denominator */
   0.4       /* pitch enhancement factor */
};

/* Post-filter parameters for wideband */
pf_params pf_params_sb = {
   0.65,      /* formant enhancement numerator */
   0.68,      /* formant enhancement denominator */
   0.3       /* pitch enhancement factor */
};

/* Parameters for Long-Term Prediction (LTP)*/
ltp_params ltp_params_nb = {
   gain_cdbk_nb,
   7,
   7
};

/* Parameters for Long-Term Prediction (LTP)*/
ltp_params ltp_params_lbr = {
   gain_cdbk_lbr,
   5,
   4
};

/* Split-VQ innovation parameters for low bit-rate narrowband */
split_cb_params split_cb_nb_lbr = {
   10,               /*subvect_size*/
   4,               /*nb_subvect*/
   exc_10_32_table, /*shape_cb*/
   5,               /*shape_bits*/
};


/* Split-VQ innovation parameters narrowband */
split_cb_params split_cb_nb = {
   5,               /*subvect_size*/
   8,               /*nb_subvect*/
   exc_5_64_table, /*shape_cb*/
   6,               /*shape_bits*/
};

/* Split-VQ innovation for low-band wideband */
split_cb_params split_cb_sb = {
   5,               /*subvect_size*/
   8,              /*nb_subvect*/
   exc_5_256_table,    /*shape_cb*/
   8,               /*shape_bits*/
};

/* Split-VQ innovation for high-band wideband */
static split_cb_params split_cb_high = {
   8,               /*subvect_size*/
   5,               /*nb_subvect*/
   hexc_table,       /*shape_cb*/
   7,               /*shape_bits*/
};

/* Default mode for narrowband */
SpeexNBMode nb_mode = {
   160,    /*frameSize*/
   40,     /*subframeSize*/
   320,    /*windowSize*/
   10,     /*lpcSize*/
   640,    /*bufSize*/
   17,     /*pitchStart*/
   144,    /*pitchEnd*/
   0,      /*lbr_pitch*/
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
   split_cb_search_nogain2,
   split_cb_nogain_unquant,
   &split_cb_nb,
   nb_post_filter,
   &pf_params_nb
};


/* Default mode for narrowband */
SpeexNBMode nb_lbr_mode = {
   160,    /*frameSize*/
   40,     /*subframeSize*/
   320,    /*windowSize*/
   10,     /*lpcSize*/
   640,    /*bufSize*/
   17,     /*pitchStart*/
   144,    /*pitchEnd*/
   1,      /*lbr_pitch*/
   0.9,    /*gamma1*/
   0.6,    /*gamma2*/
   .005,   /*lag_factor*/
   1.0001, /*lpc_floor*/
   0.0,    /*preemph*/
   /*LSP quantization*/
   lsp_quant_lbr,
   lsp_unquant_lbr,
   /*Pitch quantization*/
   pitch_search_3tap,
   pitch_unquant_3tap,
   &ltp_params_lbr,
   /*Innovation quantization*/
   split_cb_search_nogain2,
   split_cb_nogain_unquant,
   &split_cb_nb_lbr,
   nb_post_filter,
   &pf_params_lbr
};



/* Narrowband mode used for split-band wideband CELP*/
static SpeexNBMode low_sb_mode = {
   160,    /*frameSize*/
   40,     /*subframeSize*/
   320,    /*windowSize*/
   10,     /*lpcSize*/
   640,    /*bufSize*/
   17,     /*pitchStart*/
   144,    /*pitchEnd*/
   0,      /*lbr_pitch*/
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
   split_cb_search_nogain2,
   split_cb_nogain_unquant,
   &split_cb_sb,
   nb_post_filter,
   &pf_params_sb
};

SpeexMode low_wb_mode = {
   &low_sb_mode,
   -1,
   -1,
   &nb_encoder_init,
   &nb_encoder_destroy,
   &nb_encode,
   &nb_decoder_init,
   &nb_decoder_destroy,
   &nb_decode,
   &nb_encoder_ctl,
   &nb_decoder_ctl,
   160,
   17950,
   0
};

SpeexMode speex_nb_mode = {
   &nb_mode,
   0,
   0,
   &nb_encoder_init,
   &nb_encoder_destroy,
   &nb_encode,
   &nb_decoder_init,
   &nb_decoder_destroy,
   &nb_decode,
   &nb_encoder_ctl,
   &nb_decoder_ctl,
   160,
   14750,
   0
};

SpeexMode speex_nb_lbr_mode = {
   &nb_lbr_mode,
   1,
   0,
   &nb_encoder_init,
   &nb_encoder_destroy,
   &nb_encode,
   &nb_decoder_init,
   &nb_decoder_destroy,
   &nb_decode,
   &nb_encoder_ctl,
   &nb_decoder_ctl,
   160,
   7900,
   0
};

/* Split-band wideband CELP mode*/
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
   split_cb_search_shape_sign,
   split_cb_shape_sign_unquant,
   &split_cb_high
};


SpeexMode speex_wb_mode = {
   &sb_wb_mode,
   2,
   0,
   &sb_encoder_init,
   &sb_encoder_destroy,
   &sb_encode,
   &sb_decoder_init,
   &sb_decoder_destroy,
   &sb_decode,
   &sb_encoder_ctl,
   &sb_decoder_ctl,
   320,
   27350,
   0
};








void *speex_encoder_init(SpeexMode *mode)
{
   return mode->enc_init(mode);
}

void *speex_decoder_init(SpeexMode *mode)
{
   return mode->dec_init(mode);
}

void speex_encoder_destroy(void *state)
{
   (*((SpeexMode**)state))->enc_destroy(state);
}

void speex_encode(void *state, float *in, SpeexBits *bits)
{
   (*((SpeexMode**)state))->enc(state, in, bits);
}

void speex_decoder_destroy(void *state)
{
   (*((SpeexMode**)state))->dec_destroy(state);
}

void speex_decode(void *state, SpeexBits *bits, float *out, int lost)
{
   (*((SpeexMode**)state))->dec(state, bits, out, lost);
}


void speex_encoder_ctl(void *state, int request, void *ptr)
{
   (*((SpeexMode**)state))->enc_ctl(state, request, ptr);
}

void speex_decoder_ctl(void *state, int request, void *ptr)
{
   (*((SpeexMode**)state))->dec_ctl(state, request, ptr);
}
