# 
# root/CMakeLists.txt
# 
cmake_minimum_required (VERSION 3.0)

project (twMD)
# The version number.
set (twMD_VERSION_MAJOR 0)
set (twMD_VERSION_MINOR 1)

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIB_OUTPUT_PATH ${CMAKE_BINARY_DIR})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

set(CMAKE_CXX_STANDARD 11)  # enable C++11 standard

add_subdirectory(bin)
# add_subdirectory(src)