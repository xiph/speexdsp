#include "../internal.h"
#include "wrap.h"

void checkasm_check_resample(void)
{
    test_resampler_basic_direct_single();
    test_resampler_basic_direct_double();
    test_resampler_basic_interpolate_single();
    test_resampler_basic_interpolate_double();
}
