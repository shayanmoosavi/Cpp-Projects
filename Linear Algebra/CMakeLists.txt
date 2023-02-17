cmake_minimum_required(VERSION 3.24)
project(Linear_Algebra)

set(CMAKE_CXX_STANDARD 17)
set(SOURCE_FILES dep/include/vector.h source/vector.cpp
        dep/include/matrix.h source/matrix.cpp dep/include/linear_algebra.h
        source/linear_algebra.cpp)
set(CMAKE_STATIC_LIBRARY_PREFIX "")

include_directories(dep/include)

add_library(Lin_Alg STATIC ${SOURCE_FILES})
