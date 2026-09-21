if not get_option('tests').enabled()
    subdir_done()
endif

checkasm_sources = files(
    'checkasm/checkasm.c',
    'checkasm/cpu.c',
    'checkasm/resample/resample_kernels.c',
    'checkasm/resample/resampler_basic_direct_single.c',
    'checkasm/resample/resampler_basic_direct_double.c',
    'checkasm/resample/resampler_basic_interpolate_single.c',
    'checkasm/resample/resampler_basic_interpolate_double.c',
    'checkasm/resample/resample_process.c',
    'checkasm/fft/fft_bfly.c',
    'checkasm/fft/fft_transform.c',
    'checkasm/mdf/mdf_kernels.c',
    'checkasm/preprocess/preproc_kernels.c',
    'checkasm/smallft/smallft_drft.c',
    'checkasm/fftwrap/fftwrap_scale.c',
    'checkasm/filterbank/fbank_psd16.c',
)

checkasm_c_ref_args = []
if get_option('checkasm-no-autovec').enabled()
    foreach arg : ['-fno-tree-vectorize', '-fno-tree-slp-vectorize',
                   '-fno-vectorize', '-fno-slp-vectorize',
                   '-mno-implicit-float',
                   '--param=riscv-autovec-preference=none']
        if cc.has_argument(arg)
            checkasm_c_ref_args += arg
        endif
    endforeach
endif

if has_neon
    checkasm_sources += files('checkasm/resample/wrap_resample_neon.c')
endif
if has_sse and not get_option('fixed-point').enabled()
    checkasm_sources += files('checkasm/resample/wrap_resample_sse.c')
endif
if has_sse2 and not get_option('fixed-point').enabled()
    checkasm_sources += files('checkasm/resample/wrap_resample_sse2.c')
endif

checkasm_incs = [
    include_directories('../libspeexdsp'),
    config_inc,
    speexdsp_incs,
]

# Each wrapper #includes a renamed full copy of its module's .c file, so each
# needs its own static_library: the C-ref libs carry checkasm_c_ref_args
# (scalar-forcing flags only under -Dcheckasm-no-autovec=enabled); the RVV libs
# are base-ISA C and link the module's *_rvv_asm.S to supply the kernels.
# [wrapper name, checkasm/ subdir, libspeexdsp module]
checkasm_wrap_modules = [
    ['resample', 'resample',   'resample'],
    ['fft',      'fft',        'kiss_fft'],
    ['mdf',      'mdf',        'mdf'],
    ['preproc',  'preprocess', 'preprocess'],
    ['smallft',  'smallft',    'smallft'],
    ['fbank',    'filterbank', 'filterbank'],
]
checkasm_extra_libs = []
if has_rvv_asm
    # cpu.c reuses the library's RVV probe
    checkasm_sources += files('../libspeexdsp/rvv_cpu.c')
endif
foreach m : checkasm_wrap_modules
    checkasm_extra_libs += static_library('checkasm_' + m[0] + '_c',
        files('checkasm/' + m[1] + '/wrap_' + m[0] + '_c.c'),
        speexdsp_config_types_h,
        include_directories: checkasm_incs,
        c_args: checkasm_c_ref_args,
        dependencies: [checkasm_dep],
    )
    if has_rvv_asm
        checkasm_extra_libs += static_library('checkasm_' + m[0] + '_rvv',
            files('checkasm/' + m[1] + '/wrap_' + m[0] + '_rvv.c',
                  '../libspeexdsp/' + m[2] + '_rvv_asm.S'),
            speexdsp_config_types_h,
            include_directories: checkasm_incs,
            dependencies: [checkasm_dep],
        )
    endif
endforeach

checkasm = executable('checkasm',
    checkasm_sources,
    speexdsp_config_types_h,
    include_directories: checkasm_incs,
    link_with: checkasm_extra_libs,
    dependencies: [
        checkasm_dep,
    ],
)

# SIMD Correctness gate: meson test -C build
test('checkasm', checkasm, timeout: 300)

# SIMD Benchmark: meson test --benchmark -C build
benchmark('checkasm', checkasm, timeout: 0, verbose: true, 
    args: ['--bench', '--duration=30000', '--repeat=5'])
