cmake_minimum_required(VERSION 3.18)

project(PBAB)

add_executable(mkl_test mkl_test.c mkl_wrap.c)

set(C_STANDARD 11)
set(CMAKE_C_FLAGS_DEBUG "-g -Wall -DDEBUG")
set(CMAKE_C_FLAGS_RELEASE "-O3 -ffast-math -DNDEBUG")


set(BLA_VENDOR Intel10_64lp)
find_package(BLAS REQUIRED)

target_link_libraries(mkl_test BLAS::BLAS)