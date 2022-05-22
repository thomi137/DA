cmake_minimum_required(VERSION 3.10)

# Set Project
project(BEC)

SET(CMAKE_C_COMPILER gcc)
SET(CMAKE_CXX_COMPILER g++)
file(GLOB SOURCES "*.cpp")
SET(CMAKE_CXX_FLAGS "-lm -lfftw3")
SET(CMAKE_C_FLAGS   "-lm -lfftw3")
INCLUDE_DIRECTORIES(/usr/local/include)
LINK_DIRECTORIES(/usr/lib)
add_library(fftw3 STATIC IMPORTED)
set(CMAKE_C_OUTPUT_EXTENSION_REPLACE 1)
set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE 1)

add_executable(bec_try BEC_try.cpp)