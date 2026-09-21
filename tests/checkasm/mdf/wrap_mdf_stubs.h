/* Link stubs for the fftwrap symbols mdf.c references (renamed per TU by
 * wrap_mdf_rename.h). The kernel tests never run the echo canceller, so
 * they are never called. Include right after ../wrap_impl.h. */

void *spx_fft_init(int size) { (void) size; return 0; }
void spx_fft_destroy(void *table) { (void) table; }
void spx_fft(void *table, spx_word16_t *in, spx_word16_t *out)
{ (void) table; (void) in; (void) out; }
void spx_ifft(void *table, spx_word16_t *in, spx_word16_t *out)
{ (void) table; (void) in; (void) out; }
