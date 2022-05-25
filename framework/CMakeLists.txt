cmake_minimum_required(VERSION 3.10)

# Set Project
project(bec_lattice_kick
        VERSION 0.3.0
        DESCRIPTION "Simulation of a Kicked Bose-Einstein Condensate"
        LANGUAGES CXX
        )

set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)
file(GLOB SOURCES "*.cpp")
set(CMAKE_CXX_FLAGS "-lm -llapack -lfftw3")
set(CMAKE_C_FLAGS   "-lm -llapack -lfftw3")

find_package(Boost 1.78.0 COMPONENTS program_options REQUIRED)
if(Boost_FOUND)
        message(STATUS "Boost found in ${Boost_DIR}, includes in ${Boost_INCLUDE_DIR}")
        include_directories(${Boost_INCLUDE_DIR})
else()
        message(FATAL_ERROR "Boost not found on system")
endif()

include_directories(/usr/local/include)

link_directories(/usr/local/lib)
link_directories(/usr/local/opt/lapack/lib)

add_library(fftw3 STATIC IMPORTED)
add_library(bec_groundstate STATIC src/bec_groundstate/bec_groundstate.cpp)
target_include_directories(bec_groundstate PRIVATE include/bec_groundstate/bec_groundstate.h)
set_target_properties(
        bec_groundstate
        PROPERTIES
        CXX_STANDARD 11
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
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

target_link_libraries(${PROJECT_NAME}
        ${Boost_LIBRARIES}
        bec_groundstate
        )
#[[
target_include_directories(
        ${PROJECT_NAME}
        include)]]
