cmake_minimum_required(VERSION 3.20)
project(Tridiagonal)

find_package(GTest REQUIRED)

enable_testing()

add_executable(Tridiagonal test.cpp solve.hpp)
target_link_libraries(Tridiagonal GTest::gtest_main)

add_test(Alltest Tridiagonal)
