if not get_option('tests').enabled()
    subdir_done()
endif

checkasm_sources = files(
    'checkasm/checkasm.c',
    'checkasm/cpu.c',
    'checkasm/resample/resample.c',
    'checkasm/resample/resampler_basic_direct_single.c',
    'checkasm/resample/resampler_basic_direct_double.c',
    'checkasm/resample/resampler_basic_interpolate_single.c',
    'checkasm/resample/resampler_basic_interpolate_double.c',
    'checkasm/resample/wrap_resample_c.c',
)

if has_neon
    checkasm_sources += files('checkasm/resample/wrap_resample_neon.c')
endif
if has_sse and not get_option('fixed-point').enabled()
    checkasm_sources += files('checkasm/resample/wrap_resample_sse.c')
endif
if has_sse2 and not get_option('fixed-point').enabled()
    checkasm_sources += files('checkasm/resample/wrap_resample_sse2.c')
endif

checkasm = executable('checkasm',
    checkasm_sources,
    include_directories: [
        include_directories('../libspeexdsp'),
        include_directories('../include'),
        include_directories('..'),
    ],
    dependencies: [
        checkasm_dep,
    ],
)
