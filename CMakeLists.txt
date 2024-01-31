cmake_minimum_required(VERSION 3.18)

project(PBNB)

add_executable(mkl_test mkl_test.cpp mkl_wrap.hpp mkl_wrap.cpp)

set(CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math -DNDEBUG")


set(BLA_VENDOR Intel10_64lp)
find_package(BLAS REQUIRED)

target_link_libraries(mkl_test BLAS::BLAS)