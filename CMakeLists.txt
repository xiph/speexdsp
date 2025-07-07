cmake_minimum_required(VERSION 3.10)

project(speexdsp LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(USE_FIXED_POINT "Use fixed-point arithmetic" ON)

# Define compile-time options
if(USE_FIXED_POINT)
    add_compile_definitions(FIXED_POINT=1)
else()
    add_compile_definitions(FLOATING_POINT=1)
endif()

# 'M_PI': undeclared identifier
if(WIN32)
    add_compile_definitions(_USE_MATH_DEFINES=1)
endif()


add_compile_definitions(EXPORT= USE_KISS_FFT=1)

# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/libspeexdsp
)

# Source files
file(GLOB SPEEXDSP_SOURCES libspeexdsp/*.c)
add_library(speexdsp ${SPEEXDSP_SOURCES})

# Installation rules
install(TARGETS speexdsp
    EXPORT speexdspConfig
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
install(DIRECTORY include/ DESTINATION include)
