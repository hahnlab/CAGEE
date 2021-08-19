cmake_minimum_required(VERSION 3.13)

include(CheckIncludeFile)

# set the project name
project(CAGEE VERSION 1.0)

set(OPTIMIZER_HIGH_PRECISION 1e-6)
set(OPTIMIZER_LOW_PRECISION 1e-3)
set(PHASED_OPTIMIZER_PHASE1_ATTEMPTS 4)
set(NUM_OPTIMIZER_INITIALIZATION_ATTEMPTS 100)
set(LAMBDA_PERTURBATION_STEP_SIZE 1)

CHECK_INCLUDE_FILE(getopt.h HAVE_GETOPT_H)
find_package(BLAS)
find_package(Boost COMPONENTS program_options)

configure_file(config.h.in config.h)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# add the executable
add_executable(cagee main.cpp)
add_executable(tests EXCLUDE_FROM_ALL test.cpp)

target_compile_definitions(cagee PRIVATE DOCTEST_CONFIG_DISABLE)

add_subdirectory(src)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(cagee PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(tests PRIVATE OpenMP::OpenMP_CXX)
endif()

if(Boost_FOUND)
    target_link_libraries(cagee PRIVATE Boost::program_options)
    target_link_libraries(tests PRIVATE Boost::program_options)
endif()

set(CMAKE_CXX_FLAGS "-include config.h")

install(TARGETS cagee DESTINATION bin)


