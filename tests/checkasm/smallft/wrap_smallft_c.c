/* Pure-C reference build of the smallft radix stages plus the shared
 * lookup and stage helpers. #undef USE_RVV keeps the stages scalar; this
 * TU is the benchmark baseline, so it is built with the no-autovec flags
 * (checkasm_c_ref_args in tests/meson.build). */
#define CKA_PREFIX ckasc_
#include "wrap_smallft_rename.h"
#include "wrap.h"

#undef USE_RVV

#include "wrap_smallft_impl.h"

/* ------------- lookup + stage helpers ------------- */

struct drft_lookup *smallft_make_lookup(int n)
{
    struct drft_lookup *l = (struct drft_lookup *) speex_alloc(sizeof(*l));
    if (l)
        spx_drft_init(l, n);
    return l;
}

void smallft_free_lookup(struct drft_lookup *l)
{
    spx_drft_clear(l);
    speex_free(l);
}

/* Mirrors drftf1's (backward==0) / drftb1's (backward==1) stage walk,
 * recording the radix-2/4 stages with the arguments the driver passes.
 * dradfg/dradbg stages are walked (their iw advance matters) but not
 * recorded. */
int smallft_stages(struct drft_lookup *l, int backward,
                   struct drft_stage *out, int max)
{
    const int *ifac = l->splitcache;
    const int n = l->n, nf = ifac[1];
    int cnt = 0, k1;

    if (!backward) {
        int l2 = n, iw = n;
        for (k1 = 0; k1 < nf; k1++) {
            const int ip = ifac[nf - k1 + 1];
            const int l1 = l2 / ip;
            const int ido = n / l2;
            iw -= (ip - 1) * ido;
            if ((ip == 2 || ip == 4) && cnt < max) {
                out[cnt].ip = ip; out[cnt].ido = ido;
                out[cnt].l1 = l1; out[cnt].iw = iw;
                cnt++;
            }
            l2 = l1;
        }
    } else {
        int l1 = 1, iw = 1;
        for (k1 = 0; k1 < nf; k1++) {
            const int ip = ifac[k1 + 2];
            const int l2 = ip * l1;
            const int ido = n / l2;
            if ((ip == 2 || ip == 4) && cnt < max) {
                out[cnt].ip = ip; out[cnt].ido = ido;
                out[cnt].l1 = l1; out[cnt].iw = iw;
                cnt++;
            }
            l1 = l2;
            iw += (ip - 1) * ido;
        }
    }
    return cnt;
}

/* ------------- Stage driver -------------
 * Derives the wa pointers from the lookup exactly as drftf1/drftb1 do
 * and calls the (SPX_DRAD*-dispatched) stage. */
void smallft_stage_c(struct drft_lookup *l, const struct drft_stage *st,
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

/* ------------- Whole transform (in place; st/unused ignored) ------------- */
void smallft_forward_c(struct drft_lookup *l, const struct drft_stage *st,
                       int backward, float *data, float *unused)
{
    (void) st;
    (void) unused;
    if (backward)
        spx_drft_backward(l, data);
    else
        spx_drft_forward(l, data);
}
