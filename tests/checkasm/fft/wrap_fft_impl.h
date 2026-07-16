/* #includes libspeexdsp/kiss_fft.c privately into a wrap_fft_*.c TU.
 * Include AFTER wrap_fft_rename.h + wrap.h and any per-variant
 * #undef/#define. #undef HAVE_CONFIG_H stops config.h from re-defining the
 * USE_* macros the caller just cleared. */

#undef HAVE_CONFIG_H

#if defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunused-function"
#endif

#include "../../../libspeexdsp/kiss_fft.c"

#if defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif
