# Sources for compilation in the library
speexdsp_sources = files([
  'preprocess.c',
  'jitter.c',
  'mdf.c',
  'fftwrap.c',
  'filterbank.c',
  'resample.c',
  'buffer.c',
  'scal.c',
])

if fft_used == 'kiss'
  # Use KISS Fast Fourier Transform
  speexdsp_sources += files(['kiss_fft.c', 'kiss_fftr.c'])
elif fft_used == 'smallft'
  # Use FFT from OggVorbis
  speexdsp_sources += files(['smallft.c'])
endif

if has_rvv_asm
  speexdsp_sources += files(['resample_rvv_asm.S', 'mdf_rvv_asm.S'])
  if fft_used == 'kiss'
    speexdsp_sources += files(['kiss_fft_rvv_asm.S'])
  elif fft_used == 'smallft'
    speexdsp_sources += files(['smallft_rvv_asm.S'])
  endif
endif

libm = cc.find_library('m', required: false)

libspeexdsp = library('speexdsp', speexdsp_sources, speexdsp_config_types_h,
  include_directories: [config_inc, speexdsp_incs],
  dependencies: [fft_dep, libm],
  gnu_symbol_visibility: 'hidden',
  vs_module_defs: 'speexdsp.def',
  version: libversion,
  install: true)

speexdsp_dep = declare_dependency(link_with: libspeexdsp,
  include_directories: speexdsp_incs)

meson.override_dependency('speexdsp', speexdsp_dep)

# pkg-config file
pkg = import('pkgconfig')
pkg.generate(libspeexdsp,
  name: 'speexdsp',
  description: 'SpeexDSP is a library for processing audio',
  version: meson.project_version())

# Test binaries
if not get_option('examples').disabled()
  foreach bin : ['testdenoise', 'testecho', 'testjitter', 'testresample', 'testresample2', 'bench_echo']
    executable(bin, '@0@.c'.format(bin),
      include_directories: config_inc,
      dependencies: [speexdsp_dep, libm],
      install: false)
  endforeach
endif