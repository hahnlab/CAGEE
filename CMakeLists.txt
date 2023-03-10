cmake_minimum_required(VERSION 3.13)

include(CheckIncludeFile)

# set the project name
project(CAGEE VERSION 1.1)

set (CMAKE_PREFIX_PATH      "${PROJECT_SOURCE_DIR}/cmake")

set(OPTIMIZER_HIGH_PRECISION 1e-6)
set(OPTIMIZER_LOW_PRECISION 1e-3)
set(PHASED_OPTIMIZER_PHASE1_ATTEMPTS 4)
set(NUM_OPTIMIZER_INITIALIZATION_ATTEMPTS 1)
set(LAMBDA_PERTURBATION_STEP_SIZE 1)
set(LOG_OFFSET 1.0)

if (CMAKE_HOST_WIN32)
	set(NOMINMAX 1)
endif()

find_package(MKL QUIET)
find_package(GSL QUIET)
find_package(OpenBLAS QUIET)

CHECK_INCLUDE_FILE(getopt.h HAVE_GETOPT_H)

find_package (Eigen3 3.4 REQUIRED NO_MODULE)


find_package(Boost REQUIRED COMPONENTS program_options)

configure_file(config.h.in config.h)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(OPTIMIZER_STRATEGY "NelderMead" CACHE STRING "Optimizer strategy (LBFGS, Nelder-Mead, some others)")
set(DISCRETIZATION_RANGE "200" CACHE STRING "How finely to divide the continuous range")
option(MODEL_GENE_EXPRESSION_LOGS "Perform calculations using log values rather than linear values" ON)
option(USE_MAX_PROBABILITY "Evaluate probability as maximum value rather than cumulative" OFF)
option(INCLUDE_GIT_IN_VERSION_INFO "Write information about the last GIT commit into version info" OFF)

if (MODEL_GENE_EXPRESSION_LOGS)
  add_definitions(-DMODEL_GENE_EXPRESSION_LOGS)
endif()

if (USE_MAX_PROBABILITY)
  add_definitions(-DUSE_MAX_PROBABILITY)
endif()

add_definitions(-DOPTIMIZER_STRATEGY=${OPTIMIZER_STRATEGY})
add_definitions(-DDISCRETIZATION_RANGE=${DISCRETIZATION_RANGE})


if (INCLUDE_GIT_IN_VERSION_INFO)
  # Define the two required variables before including
  # the source code for watching a git repository.
  set(PRE_CONFIGURE_FILE "git_version.cpp.in")
  set(POST_CONFIGURE_FILE "${CMAKE_CURRENT_BINARY_DIR}/git_version.cpp")
  include(git_watcher.cmake)

  # Create a library out of the compiled post-configure file.
  add_definitions(-DINCLUDE_GIT_IN_VERSION_INFO)
  add_library(git STATIC ${POST_CONFIGURE_FILE})
  target_include_directories(git PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
  add_dependencies(git check_git)
endif()

# add the executable
add_executable(cagee main.cpp)
add_executable(tests EXCLUDE_FROM_ALL test.cpp)
add_executable(diffmat_precalc diffmat_precalc.cpp)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(LBFGSpp)

target_compile_definitions(cagee PRIVATE DOCTEST_CONFIG_DISABLE)
target_compile_definitions(diffmat_precalc PRIVATE DOCTEST_CONFIG_DISABLE)
target_compile_definitions(tests PRIVATE SILENT ELPP_NO_CHECK_MACROS)
target_compile_options(tests PRIVATE -Wno-range-loop-construct)
add_subdirectory(src)

link_libraries(Eigen3::Eigen)

if (INCLUDE_GIT_IN_VERSION_INFO)
  target_link_libraries(cagee PRIVATE git)
  target_link_libraries(tests PRIVATE git)
endif()

find_package(ZLIB)
target_link_libraries(cagee PRIVATE ZLIB::ZLIB)
target_link_libraries(tests PRIVATE ZLIB::ZLIB)
target_link_libraries(diffmat_precalc PRIVATE ZLIB::ZLIB)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(cagee PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(tests PRIVATE OpenMP::OpenMP_CXX)
endif()

if(Boost_FOUND)
    target_link_libraries(cagee PRIVATE Boost::program_options)
    target_link_libraries(tests PRIVATE Boost::program_options)
endif()

if(MKL_FOUND)
    target_include_directories(cagee PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
    target_include_directories(tests PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
    target_include_directories(diffmat_precalc PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)

    target_compile_options(cagee PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
    target_compile_options(tests PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
    target_compile_options(diffmat_precalc PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)

    target_link_libraries(cagee PUBLIC $<LINK_ONLY:MKL::MKL>)
    target_link_libraries(tests PUBLIC $<LINK_ONLY:MKL::MKL>)
    target_link_libraries(diffmat_precalc PUBLIC $<LINK_ONLY:MKL::MKL>)
    add_compile_definitions(BLAS_AVAILABLE)
elseif(OpenBLAS_FOUND)
    target_include_directories(cagee PUBLIC ${OpenBLAS_INCLUDE_DIRS})
    target_include_directories(tests PUBLIC ${OpenBLAS_INCLUDE_DIRS})
    target_include_directories(diffmat_precalc PUBLIC ${OpenBLAS_INCLUDE_DIRS})

    target_link_libraries(cagee PUBLIC ${OpenBLAS_LIBRARIES})
    target_link_libraries(tests PUBLIC ${OpenBLAS_LIBRARIES})
    target_link_libraries(diffmat_precalc PUBLIC ${OpenBLAS_LIBRARIES})
    add_compile_definitions(BLAS_AVAILABLE)
elseif (GSL_FOUND)
    target_include_directories(cagee PUBLIC ${GSL_INCLUDE_DIRS})
    target_include_directories(tests PUBLIC ${GSL_INCLUDE_DIRS})
    target_include_directories(diffmat_precalc PUBLIC ${GSL_INCLUDE_DIRS})

    target_link_libraries(cagee PUBLIC ${GSL_LIBRARIES})
    target_link_libraries(tests PUBLIC ${GSL_LIBRARIES})
    target_link_libraries(diffmat_precalc PUBLIC ${GSL_LIBRARIES})
    add_compile_definitions(BLAS_AVAILABLE)
endif()

# GCC or Clang
add_definitions(-include config.h)
target_compile_options(tests PRIVATE -g -rdynamic -Wall -Werror)


install(TARGETS cagee DESTINATION bin)
install(TARGETS diffmat_precalc DESTINATION bin)
install(DIRECTORY examples DESTINATION .)
install(DIRECTORY docs/manual DESTINATION .)
install(FILES README.md LICENSE CHANGELOG.md DESTINATION .)

include(CPack)

