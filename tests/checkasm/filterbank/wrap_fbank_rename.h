/* Rename filterbank.c's public API to a per-TU-unique prefix so each
 * wrap_fbank_*.c can #include the whole file without link-time
 * collisions -- same scheme as ../preprocess/wrap_preproc_rename.h.
 * #define CKA_PREFIX and include this BEFORE wrap.h; it also names the
 * module source for ../wrap_impl.h. */
#include "../wrap_common.h"
#define CKA_MODULE_SRC "filterbank.c"

#define filterbank_new            CKA_CAT(CKA_PREFIX, filterbank_new)
#define filterbank_destroy        CKA_CAT(CKA_PREFIX, filterbank_destroy)
#define filterbank_compute_bank32 CKA_CAT(CKA_PREFIX, filterbank_compute_bank32)
#define filterbank_compute_psd16  CKA_CAT(CKA_PREFIX, filterbank_compute_psd16)
#define filterbank_compute_bank   CKA_CAT(CKA_PREFIX, filterbank_compute_bank)
#define filterbank_compute_psd    CKA_CAT(CKA_PREFIX, filterbank_compute_psd)
#define filterbank_psy_smooth     CKA_CAT(CKA_PREFIX, filterbank_psy_smooth)
