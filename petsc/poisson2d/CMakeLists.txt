#
#  This is a sample CMakeLists.txt suitable for direct use with a new PETSc application or
#  you can add fragments of material below to an existing applications CMakeLists.txt
#
#  Usage:  To build ex1 from ex1.c
#    rm -fr build
#    mkdir build
#    cd build
#    cmake ..
#    make
#    ./ex1
#
#  By default it gets both the compiler information and the library information from PETSc
#
cmake_minimum_required(VERSION 3.1.0)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
PROJECT(poisson2d CXX)
include(CTest)
include(CMakePrintHelpers)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 11)
set(MAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_VERBOSE_MAKEFILE TRUE)
set(CMAKE_VERBOSE_MAKEFILE FALSE)

# set root of location to find PETSc's pkg-config
set(PETSC $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})
set(ENV{PKG_CONFIG_PATH} ${PETSC}/lib/pkgconfig)

# Remove the lines below if you do not wish to have PETSc determine the compilers
#execute_process ( COMMAND pkg-config PETSc --variable=ccompiler COMMAND tr -d '\n' OUTPUT_VARIABLE C_COMPILER)
#SET(CMAKE_C_COMPILER ${C_COMPILER})
#execute_process ( COMMAND pkg-config PETSc --variable=cxxcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE CXX_COMPILER)
#if (CXX_COMPILER)
#  SET(CMAKE_CXX_COMPILER ${CXX_COMPILER})
#endif (CXX_COMPILER)
#execute_process ( COMMAND pkg-config PETSc --variable=fcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE FORTRAN_COMPILER)
#if (FORTRAN_COMPILER)
#  SET(CMAKE_Fortran_COMPILER ${FORTRAN_COMPILER})
#  enable_language(Fortran)
#endif (FORTRAN_COMPILER)

#  tells CMake to build the application ex1 from the source file ex1.c
#  this must appear AFTER the compilers are set
project(poisson2d)
add_executable(poisson2d main.cc)
target_compile_options(poisson2d PRIVATE
    $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
    -Wall -Wextra -pedantic>
    $<$<CXX_COMPILER_ID:MSVC>:
    /W4>)

find_package(PkgConfig REQUIRED)
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
find_package(MPI REQUIRED)

set(EXTRA_DEBUG_FLAGS STRING "-errors -fsanitize=address -fno-omit-frame-pointer -fsanitize-recover=address")

list(APPEND COMPILE_OPTIONS ${MPI_CFLAGS})
# TODO it's possible that these flags are never used. Will have to check once something starts segfaulting again.
list(APPEND CMAKE_CXX_COMPILE_OPTIONS_DEBUG ${EXTRA_DEBUG_FLAGS})
include_directories(${MPI_INCLUDE_PATH})
set(LINK_FLAGS "${LINK_FLAGS} ${MPI_LDFLAGS} ${EXTRA_DEBUG_FLAGS}")
list(APPEND LIBRARIES ${MPI_C_LIBRARIES})
set(CMAKE_REQUIRED_FLAGS ${MPI_CFLAGS})   
set(CMAKE_REQUIRED_INCLUDES "${MPI_INCLUDE_DIRS}")

target_link_libraries(poisson2d PkgConfig::PETSC MPI::MPI_CXX ${LIBRARIES})

# Tests 
add_test(NAME Compile_and_Run  COMMAND ./poisson2d)
