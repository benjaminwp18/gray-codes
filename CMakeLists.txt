cmake_minimum_required(VERSION 3.10)
project(gray-codes CXX)

find_package(OpenMP REQUIRED)

# For debugging
# set(CMAKE_BUILD_TYPE Debug)
# set(DCMAKE_C_FLAGS_DEBUG, "-g -O0")
# set(DCMAKE_CXX_FLAGS_DEBUG, "-g -O0")

include_directories("include")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(gray_codes src/gray_codes.cpp include/gray_codes.hpp)
# target_link_libraries(gray_codes PUBLIC OpenMP::OpenMP_CXX)

add_executable(main src/main.cpp)
target_link_libraries(main PRIVATE gray_codes OpenMP::OpenMP_CXX)
