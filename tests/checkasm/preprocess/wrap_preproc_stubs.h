/* Link stubs for the fftwrap/filterbank/mdf symbols preprocess.c
 * references (renamed per TU by wrap_preproc_rename.h). The kernel tests
 * never run the preprocessor itself, so they are never called. Include
 * right after ../wrap_impl.h. */

void *spx_fft_init(int size) { (void) size; return 0; }
void spx_fft_destroy(void *table) { (void) table; }
void spx_fft(void *table, spx_word16_t *in, spx_word16_t *out)
{ (void) table; (void) in; (void) out; }
void spx_ifft(void *table, spx_word16_t *in, spx_word16_t *out)
{ (void) table; (void) in; (void) out; }
FilterBank *filterbank_new(int banks, spx_word32_t sampling, int len, int type)
{ (void) banks; (void) sampling; (void) len; (void) type; return 0; }
void filterbank_destroy(FilterBank *bank) { (void) bank; }
void filterbank_compute_bank32(FilterBank *bank, spx_word32_t *ps, spx_word32_t *mel)
{ (void) bank; (void) ps; (void) mel; }
void filterbank_compute_psd16(FilterBank *bank, spx_word16_t *mel, spx_word16_t *psd)
{ (void) bank; (void) mel; (void) psd; }
void speex_echo_get_residual(SpeexEchoState *st, spx_word32_t *Yout, int len)
{ (void) st; (void) Yout; (void) len; }
