/* #includes a libspeexdsp module privately into a wrap_<mod>_<variant>.c TU.
 * The module's wrap_<mod>_rename.h sets CKA_MODULE_SRC to its file name
 * (found through the libspeexdsp include path). Include this AFTER that
 * rename header + wrap.h and any per-variant #undef USE_* / *_FORCE_ON.
 *
 * #undef HAVE_CONFIG_H stops the module (or os_support.h) from re-including
 * config.h, which would re-define the USE_* macros the caller just cleared;
 * HAVE_CONFIG_H is set on the command line, but the #undef still takes
 * effect for the rest of this TU. The full module pulls in its whole
 * public API plus internal helpers, many of which a given TU never calls,
 * so the pragma quiets the resulting -Wunused-function noise. */

#ifndef CKA_MODULE_SRC
#  error "CKA_MODULE_SRC not set; include the module's wrap_*_rename.h first"
#endif

#undef HAVE_CONFIG_H

#if defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunused-function"
#endif

#include CKA_MODULE_SRC

#if defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif
