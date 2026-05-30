#include "../internal.h"
#include "wrap.h"

void checkasm_check_resample(void)
{
    test_inner_product_single();
    test_inner_product_double();
    test_interpolate_product_single();
    test_interpolate_product_double();
}
