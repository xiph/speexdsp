/* #includes libspeexdsp/resample.c privately into a wrap_resample_*.c TU.
 *
 * Include this AFTER wrap_resample_rename.h + wrap.h (so the public API is
 * renamed and the prototypes are in scope) and AFTER any per-variant
 * USE_SSE/USE_SSE2/USE_NEON #undefs that select the kernel build.
 *
 * #undef HAVE_CONFIG_H stops resample.c / os_support.h from re-#including
 * config.h, which would otherwise re-define the USE_* macros the caller just
 * cleared (HAVE_CONFIG_H is set on the command line, but the #undef still takes
 * effect for the rest of this TU). The full file pulls in the whole public API
 * plus internal helpers; in a given TU many go uncalled, so the surrounding
 * pragma quiets the resulting -Wunused-function noise. */

#undef HAVE_CONFIG_H

#if defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunused-function"
#endif

#include "../../../libspeexdsp/resample.c"

#if defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif
