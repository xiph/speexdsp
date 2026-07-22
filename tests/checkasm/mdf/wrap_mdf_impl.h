/* #includes libspeexdsp/mdf.c privately into a wrap_mdf_*.c TU.
 * Include AFTER wrap_mdf_rename.h + wrap.h and any per-variant
 * #undef/#define. #undef HAVE_CONFIG_H stops config.h from re-defining the
 * USE_* macros the caller just cleared. */

#undef HAVE_CONFIG_H

#if defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunused-function"
#endif

#include "../../../libspeexdsp/mdf.c"

#if defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif

/* The kernels under test never touch the FFT; these (renamed, per-TU)
 * stubs only satisfy the link for mdf.c's unused public API. */
void *spx_fft_init(int size) { (void) size; return 0; }
void spx_fft_destroy(void *table) { (void) table; }
void spx_fft(void *table, spx_word16_t *in, spx_word16_t *out)
{ (void) table; (void) in; (void) out; }
void spx_ifft(void *table, spx_word16_t *in, spx_word16_t *out)
{ (void) table; (void) in; (void) out; }
