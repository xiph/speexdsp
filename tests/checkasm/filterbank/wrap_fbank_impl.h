/* #includes libspeexdsp/filterbank.c privately into a wrap_fbank_*.c
 * TU. Include AFTER wrap_fbank_rename.h + wrap.h and any per-variant
 * #undef/#define. #undef HAVE_CONFIG_H stops config.h from re-defining
 * the USE_* macros the caller just cleared. */

#undef HAVE_CONFIG_H

#if defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunused-function"
#endif

#include "../../../libspeexdsp/filterbank.c"

#if defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif
