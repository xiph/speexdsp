/* RVV build of the smallft radix stages for checkasm: base-ISA C, the V
 * instructions live in smallft_rvv_asm.S (linked alongside).
 * SMALLFT_RVV_FORCE_ON pins the dispatch on so the asm is tested
 * deterministically, not per the build host's getauxval. */
#define CKA_PREFIX ckasrvv_
#include "wrap_smallft_rename.h"
#include "wrap.h"

#ifdef USE_RVV
#  define SMALLFT_RVV_FORCE_ON
#endif

#include "wrap_smallft_impl.h"

#ifdef HAVE_RVV_SMALLFT

void smallft_stage_rvv(struct drft_lookup *l, const struct drft_stage *st,
                       int backward, float *cc, float *ch)
{
    float *wa = l->trigcache + l->n;
    float *wa1 = wa + st->iw - 1;
    float *wa2 = wa + st->iw + st->ido - 1;
    float *wa3 = wa + st->iw + 2 * st->ido - 1;

    if (st->ip == 4) {
        if (backward)
            SPX_DRADB4(st->ido, st->l1, cc, ch, wa1, wa2, wa3);
        else
            SPX_DRADF4(st->ido, st->l1, cc, ch, wa1, wa2, wa3);
    } else {
        if (backward)
            SPX_DRADB2(st->ido, st->l1, cc, ch, wa1);
        else
            SPX_DRADF2(st->ido, st->l1, cc, ch, wa1);
    }
}

void smallft_forward_rvv(struct drft_lookup *l, const struct drft_stage *st,
                         int backward, float *data, float *unused)
{
    (void) st;
    (void) unused;
    if (backward)
        spx_drft_backward(l, data);
    else
        spx_drft_forward(l, data);
}

#endif /* HAVE_RVV_SMALLFT */
