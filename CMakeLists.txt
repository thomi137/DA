cmake_minimum_required(VERSION 3.10)

# Set Project
project(bec_lattice_kick
        VERSION 0.3.0
        DESCRIPTION "Simulation of a Kicked Bose-Einstein Condensate"
        LANGUAGES CXX
        )

set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_STANDARD 11)
Set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-lm -lfftw3")
set(CMAKE_C_FLAGS   "-lm -lfftw3")

file(GLOB SOURCES "*.cpp")

include_directories(include)

# NOTE: Invoke cmake with -DBLAS_LIBRARIES=<path_to_blas> -DLAPACK_LIBRARIES=<path_to_blas
# This is necessary as the built-in module does not find Apple/Homebrew specifics.
find_package(LAPACK)
if(LAPACK_FOUND)
        message(STATUS "Lapack found in ${LAPACK_LIBRARIES}")
else()
        message(FATAL_ERROR "Lapack not found on system")
endif()

find_package(Boost 1.78.0 REQUIRED COMPONENTS program_options)
if(Boost_FOUND)
        message(STATUS "Boost found in ${Boost_DIR}, includes in ${Boost_INCLUDE_DIR}")
        message(STATUS "Boost libraries: ${Boost_LIBRARIES}")
        include_directories(${Boost_INCLUDE_DIR})
else()
        message(FATAL_ERROR "Boost not found on system")
endif()

link_directories(/usr/local/lib)

add_library(fftw3 STATIC IMPORTED)

add_library(bec_groundstate STATIC src/bec_groundstate/bec_groundstate.cpp)
target_include_directories(bec_groundstate PUBLIC include/bec_groundstate)
target_link_libraries(
        bec_groundstate
        PUBLIC Boost::program_options
        LAPACK
        )


set(CMAKE_C_OUTPUT_EXTENSION_REPLACE 1)
set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE 1)

add_executable(${PROJECT_NAME} src/main.cpp)
set_target_properties(
        ${PROJECT_NAME}
        PROPERTIES
        CXX_STANDARD 11
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
        )

target_link_libraries(
        ${PROJECT_NAME}
        bec_groundstate
        )
#[[
target_include_directories(
        ${PROJECT_NAME}
        include)]]
