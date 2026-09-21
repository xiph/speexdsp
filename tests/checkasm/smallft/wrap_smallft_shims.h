/* External-linkage shims around smallft.c's static radix stages and the
 * whole transform. Included by wrap_smallft_c.c and wrap_smallft_rvv.c
 * after the module, with CKA_VARIANT naming the TU; wrap.h carries the
 * matching prototypes. The stage shim derives the wa pointers from the
 * lookup exactly as drftf1/drftb1 do and calls the (SPX_DRAD*-dispatched)
 * stage. */

void CKA_SHIM(smallft_stage)(struct drft_lookup *l, const struct drft_stage *st,
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

/* Whole transform, in place; st/unused only pad the signature to the
 * stage shim's. */
void CKA_SHIM(smallft_forward)(struct drft_lookup *l, const struct drft_stage *st,
                               int backward, float *data, float *unused)
{
    (void) st;
    (void) unused;
    if (backward)
        spx_drft_backward(l, data);
    else
        spx_drft_forward(l, data);
}
