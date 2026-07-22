/* #includes libspeexdsp/smallft.c privately into a wrap_smallft_*.c TU.
 * Include AFTER wrap_smallft_rename.h + wrap.h and any per-variant
 * #undef/#define. #undef HAVE_CONFIG_H stops config.h from re-defining
 * the USE_* macros the caller just cleared. smallft.c is self-contained
 * (allocators are static inline), so no stubs are needed. */

#undef HAVE_CONFIG_H

#if defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunused-function"
#endif

#include "../../../libspeexdsp/smallft.c"

#if defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif
