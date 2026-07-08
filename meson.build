project('speexdsp', 'c', version: '1.2.1',
  meson_version: '>= 0.63')

libversion = '1.5.2'

speexdsp_version = meson.project_version()
version_arr = speexdsp_version.split('.')
speexdsp_version_major = version_arr[0].to_int()
speexdsp_version_minor = version_arr[1].to_int()

# extra version might be "-svn" or ".2" or empty
if version_arr.length() == 4
  speexdsp_version_micro = version_arr[2].to_int()
  speexdsp_version_extra = '.' + version_arr[3]
elif version_arr[2].contains('-')
  speexdsp_version_micro = version_arr[2].split('-')[0].to_int()
  speexdsp_version_extra = '-' + version_arr[2].split('-')[1]
else
  speexdsp_version_micro = version_arr[2].to_int()
  speexdsp_version_extra = ''
endif

host_system = host_machine.system()
host_cpu_family = host_machine.cpu_family()

cc = meson.get_compiler('c')

config_inc = include_directories('.')

speexdsp_incs = include_directories('include', 'include/speex')

cdata = configuration_data()

add_project_arguments(
  '-Wno-unused-but-set-variable',
  '-Wno-unused-variable',
  '-Wno-pointer-sign',
  language : 'c'  # build warning suppression
)

add_project_arguments('-DHAVE_CONFIG_H=1', language: 'c')

if host_machine.endian() == 'big'
  cdata.set('WORDS_BIGENDIAN', true)
endif

cdata.set_quoted('SPEEXDSP_VERSION', speexdsp_version,
                 description: 'Complete version string')

cdata.set('SPEEXDSP_MAJOR_VERSION', speexdsp_version_major,
          description: 'Version major')
cdata.set('SPEEXDSP_MINOR_VERSION', speexdsp_version_minor,
          description: 'Version minor')
cdata.set('SPEEXDSP_MICRO_VERSION', speexdsp_version_micro,
          description: 'Version micro')
cdata.set_quoted('SPEEXDSP_EXTRA_VERSION', speexdsp_version_extra,
          description: 'Version extra (version suffix string)')

# Test whether the compiler supports the C99 varaiable size arrays.
have_var_arrays = cc.compiles('void func (int n) {int arr[n]; }',
                              name: 'C99 variable-size arrays')

cdata.set('VAR_ARRAYS', have_var_arrays,
          description: 'Use C99 variable-size arrays')

has_sse = false
if host_cpu_family in ['x86', 'x86_64'] and get_option('sse').allowed()
  intrin_check = '''#include <xmmintrin.h>
                    __m128 testfunc(float *a, float *b) {
                      return _mm_add_ps(_mm_loadu_ps(a), _mm_loadu_ps(b));
                    }
                    int main (int argc, char ** argv) {
                      float a = 2.0, b = 99.0;
                      testfunc(&a, &b);
                      return 0;
                    }
                 '''
  # Intrinsics arguments are not available with MSVC-like compilers
  intrin_args = cc.get_argument_syntax() == 'msvc' ? [] : ['-msse']
  if cc.links(intrin_check, name: 'Compiler supports SSE intrinsics', args: intrin_args)
    has_sse = true
  endif

  # Mimicing the logic in win32/config.h
  if cc.get_id() == 'msvc' and host_cpu_family == 'x86'
    m_ix86_fp = cc.get_define('_M_IX86_FP')
    if m_ix86_fp == '' or m_ix86_fp.to_int() < 1
      has_sse = false
    endif
  endif
endif

has_sse2 = false
if host_cpu_family in ['x86', 'x86_64'] and get_option('sse2').allowed()
  intrin_check = '''#include <emmintrin.h>
                    __m128d testfunc(double *a, double *b) {
                      return _mm_add_pd(_mm_loadu_pd(a), _mm_loadu_pd(b));
                    }
                    int main (int argc, char ** argv) {
                      double a = 2.0, b = 99.0;
                      testfunc(&a, &b);
                      return 0;
                    }
                 '''
  # Intrinsics arguments are not available with MSVC-like compilers
  intrin_args = cc.get_argument_syntax() == 'msvc' ? [] : ['-msse2']
  if cc.links(intrin_check, name: 'Compiler supports SSE2 intrinsics', args: intrin_args)
    has_sse2 = true
  endif

  # Mimicing the logic in win32/config.h
  if cc.get_id() == 'msvc' and host_cpu_family == 'x86'
    m_ix86_fp = cc.get_define('_M_IX86_FP')
    if m_ix86_fp == '' or m_ix86_fp.to_int() < 2
      has_sse2 = false
    endif
  endif
endif

has_neon = false
if host_cpu_family in ['arm', 'aarch64'] and get_option('neon').allowed()
  intrin_check = '''#include <arm_neon.h>
                    int32x4_t testfunc(int16_t *a, int16_t *b) {
                      return vmull_s16(vld1_s16(a), vld1_s16(b));
                    }
                    int main (int argc, char ** argv) {
                      int16_t a[4] = {2, 0, 0, 0}, b[4] = {99, 0, 0, 0};
                      testfunc(a, b);
                      return 0;
                    }
                 '''
  neon_args = []
  if host_cpu_family == 'arm'
    neon_args = cc.get_argument_syntax() == 'msvc' ? [] : ['-mfpu=neon']
  endif
  if cc.links(intrin_check, name: 'Compiler supports NEON intrinsics', args: neon_args)
    has_neon = true
    if host_cpu_family == 'arm' and neon_args.length() > 0
      add_project_arguments(neon_args, language: 'c')
    endif
  endif
endif

# The .S enables RVV via ".option arch, +v" rather than a "-march" flag on the C
# compiler; probe that the assembler accepts it. RVV 1.0 has no per-instruction
# subsetting, so a single instruction suffices.
has_rvv_asm = false
if host_cpu_family.startswith('riscv') and get_option('rvv').allowed()
  rvv_check = '''__asm__ (
".option arch, +v\n"
"vsetivli zero, 0, e8, m1, ta, ma"
);
'''
  if cc.compiles(rvv_check, name: 'RVV (.option arch, +v)')
    has_rvv_asm = true
  endif
endif


if cc.get_id() != 'msvc' and cc.has_argument('-fvisibility=hidden')
  add_project_arguments('-fvisibility=hidden', language: 'c')
  cdata.set('EXPORT', '__attribute__((visibility("default")))')
else
  cdata.set('EXPORT', '/* */')
endif

# FFT implementation

fft_dep = []

use_fixed_point = get_option('fixed-point').enabled()

if use_fixed_point
  default_fft = 'kiss'
  has_sse = false
  has_sse2 = false
  cdata.set('FIXED_POINT', true)
else
  default_fft = 'smallft'
  cdata.set('FLOATING_POINT', true)
endif

fft_opt = get_option('fft')

if fft_opt == 'auto'
  fft_used = default_fft
  if default_fft == 'kiss'
    cdata.set('USE_KISS_FFT', true)
  elif default_fft == 'smallft'
    cdata.set('USE_SMALLFT', true)
  endif
elif fft_opt == 'kiss'
  fft_used = 'kiss'
  cdata.set('USE_KISS_FFT', true)
elif fft_opt == 'smallft'
  fft_used = 'smallft'
  cdata.set('USE_SMALLFT', true)
elif fft_opt == 'gpl-fftw3'
  fft_used = 'gpl-fftw3'
  fft_dep = dependency('fftw3f', required: true)
  cdata.set('USE_GPL_FFTW3', true)
elif fft_opt == 'proprietary-intel-mkl'
  error('FIXME: add support for proprietary-intel-mkl FFT')
else
  error('Unexpected FFT option', fft_opt)
endif

if not has_sse and get_option('sse').enabled()
  error('No SSE support, but SSE support required via options')
endif

cdata.set('USE_SSE', has_sse)

if not has_sse2 and get_option('sse2').enabled()
  error('No SSE2 support, but SSE2 support required via options')
endif

cdata.set('USE_SSE2', has_sse2)

cdata.set('USE_NEON', has_neon)

if not has_rvv_asm and get_option('rvv').enabled()
  error('No RVV assembler support, but RVV support required via options')
endif

cdata.set('USE_RVV', has_rvv_asm)

opt_arm4_asm = get_option('arm4-asm').enabled()
opt_arm5e_asm = get_option('arm5e-asm').enabled()
opt_blackfin_asm = get_option('blackfin-asm').enabled()
# opt_ti_c55x = get_option('ti-c55x').enabled()

cdata.set('ARM4_ASM', opt_arm4_asm)
cdata.set('ARM5E_ASM', opt_arm5e_asm)
cdata.set('BFIN_ASM', opt_blackfin_asm)
# cdata.set('TI_C55X', opt_ti_c55x) # has_char16=yes

opt_disable_float_api = get_option('float-api').disabled()
opt_fixed_point_debug = get_option('fixed-point-debug').enabled()
opt_resample_full_sinc_table = get_option('resample-full-sinc-table').enabled()

cdata.set('DISABLE_FLOAT_API', opt_disable_float_api)
cdata.set('FIXED_POINT_DEBUG', opt_fixed_point_debug)
cdata.set('RESAMPLE_FULL_SINC_TABLE', opt_resample_full_sinc_table)

# config.h
configure_file(output: 'config.h', configuration: cdata)

subdir('doc')
subdir('include/speex')
subdir('libspeexdsp')

if get_option('tests').enabled()
    checkasm_dep = dependency('checkasm',
        fallback: ['checkasm', 'checkasm_dep'],
        required: false
    )

    if checkasm_dep.found()
        subdir('tests')
    endif
endif

# summary

asms = {
  'has_sse': 'SSE',
  'has_sse2': 'SSE2',
  'has_neon': 'NEON',
  'has_rvv_asm': 'RVV (runtime-dispatched)',
  'opt_arm4_asm': 'ARM4 assembly',
  'opt_arm5e_asm': 'ARMv5e assembly',
  'opt_blackfin_asm': 'Blackfin assembly',
  #'opt_ti_c55x': 'TI C55X DSP',
}
optimizations = []
foreach var, desc : asms
  if get_variable(var, false)
    optimizations += [desc]
  endif
endforeach

summary(
  {
    'Assembly optimizations': optimizations.length() > 0 ? optimizations : false,
    'Floating point support': not use_fixed_point,
    'Fixed point debugging': use_fixed_point ? opt_fixed_point_debug : 'N/A',
    'Float API': not opt_disable_float_api,
    'FFT': fft_used,
  },
  bool_yn: true,
  list_sep: ', ',
)
