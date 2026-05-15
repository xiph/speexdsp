install_headers([
  'speex_echo.h',
  'speex_jitter.h',
  'speex_preprocess.h',
  'speex_resampler.h',
  'speexdsp_types.h',
], subdir: 'speex')

# speexdsp_config_types.h
ctypes_cdata = configuration_data()
ctypes_cdata.set('INCLUDE_STDINT', '#include <stdint.h>')
ctypes_cdata.set('SIZE16', 'int16_t')
ctypes_cdata.set('USIZE16', 'uint16_t')
ctypes_cdata.set('SIZE32', 'int32_t')
ctypes_cdata.set('USIZE32', 'uint32_t')

speexdsp_config_types_h = configure_file(input: 'speexdsp_config_types.h.in',
  output: 'speexdsp_config_types.h',
  configuration: ctypes_cdata,
  install_dir: get_option('includedir') / 'speex',
  install: true)