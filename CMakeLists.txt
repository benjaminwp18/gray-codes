cmake_minimum_required(VERSION 3.10)
project(gray-codes LANGUAGES CXX)

# set(CMAKE_BUILD_TYPE Debug)
# set(DCMAKE_C_FLAGS_DEBUG, "-g -O0")
# set(DCMAKE_CXX_FLAGS_DEBUG, "-g -O0")

include_directories("include")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(gray-codes src/main.cpp include/main.hpp)
