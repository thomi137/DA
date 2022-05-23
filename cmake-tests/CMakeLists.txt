cmake_minimum_required(VERSION 3.10)

# Set Project
project(BEC)

SET(CMAKE_C_COMPILER gcc)
SET(CMAKE_CXX_COMPILER g++)
file(GLOB SOURCES "*.cpp")
SET(CMAKE_CXX_FLAGS "-lm -llapack -lfftw3")
SET(CMAKE_C_FLAGS   "-lm -llapack -lfftw3")

set (BOOST_INCLUDEDIR "/usr/local/include/boost")

FIND_PACKAGE(Boost 1.78.0)

INCLUDE_DIRECTORIES(/usr/local/include)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

LINK_DIRECTORIES(/usr/local/lib)
LINK_DIRECTORIES(/usr/local/opt/lapack/lib)

add_library(fftw3 STATIC IMPORTED)
set(CMAKE_C_OUTPUT_EXTENSION_REPLACE 1)
set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE 1)

add_executable(bec_try BEC_try.cpp)
add_executable(bec_try2 BEC_2try.cpp)
