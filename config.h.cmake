// Use C99 variable-size arrays
#cmakedefine01 VAR_ARRAYS

// Make use of alloca
#cmakedefine01 USE_ALLOCA

// Symbol visibility prefix
#cmakedefine EXPORT

// Enable valgrind extra checks
#cmakedefine ENABLE_VALGRIND

// Compile as fixed-point
#cmakedefine FIXED_POINT

// Compile as floating-point
#cmakedefine FLOATING_POINT

// Enable SSE support
#cmakedefine USE_SSE

// Enable NEON support
#cmakedefine USE_NEON

// Enable SSE2 support
#cmakedefine USE_SSE2

// Disable all parts of the API that are using floats
#cmakedefine DISABLE_FLOAT_API

// Make use of ARM4 assembly optimizations
#cmakedefine ARM4_ASM

// Make use of ARM5E assembly optimizations
#cmakedefine ARM5E_ASM

// Resample with full SINC table (no interpolation)
#cmakedefine RESAMPLE_FULL_SINC_TABLE

// Use KISS Fast Fourier Transform
#cmakedefine USE_KISS_FFT

// Use FFT from OggVorbis
#cmakedefine USE_SMALLFT

// Use FFTW3 for FFT
#cmakedefine USE_GPL_FFTW3

// Use Intel Math Kernel Library for FFT
#cmakedefine USE_INTEL_MKL

// Used to avoid potential symbol clashes
#cmakedefine RANDOM_PREFIX
