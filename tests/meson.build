if not get_option('tests').enabled()
    subdir_done()
endif

checkasm_sources = files(
    'checkasm/checkasm.c',
    'checkasm/resample.c',
    'checkasm/resample_c_wrap.c',
    'checkasm/resample_neon_wrap.c',
)

checkasm = executable('checkasm',
    checkasm_sources,
    include_directories: [
        include_directories('../libspeexdsp'),
        include_directories('../include'),
        include_directories('checkasm'),
        include_directories('.'),
        include_directories('../'),
    ],
    dependencies: [
        checkasm_dep,
    ],
)
