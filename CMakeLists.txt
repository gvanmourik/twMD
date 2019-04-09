# 
# root/CMakeLists.txt
# 
cmake_minimum_required (VERSION 3.12.0)

project (twMD)

# The version number.
set (twMD_VERSION_MAJOR 0)
set (twMD_VERSION_MINOR 1)

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIB_OUTPUT_PATH ${CMAKE_BINARY_DIR})


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# Enable C++11 standard
set(CMAKE_CXX_STANDARD 11)

# Boost setup
find_package(Boost REQUIRED)
if(Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS})
endif()

# MPI check
find_package(MPI REQUIRED)
if(MPI_FOUND)
  include_directories(SYSTEM ${MPI_INCLUDE_PATH})
endif()

# OpenMP check
if(APPLE)
    # Set the clang compiler
    message(STATUS "Setting compiler...")
    set(CMAKE_C_COMPILER "/usr/local/Cellar/llvm/8.0.0/bin/clang")
    set(CMAKE_CXX_COMPILER "/usr/local/Cellar/llvm/8.0.0/bin/clang++")
    set(OPENMP_LIBRARIES "/usr/local/Cellar/llvm/8.0.0/lib")
    set(OPENMP_INCLUDES "/usr/local/Cellar/llvm/8.0.0/include")

    # message(${CMAKE_C_COMPILER_ID})
    if(CMAKE_C_COMPILER_ID MATCHES "Clang")
        message(STATUS "C compiler: " ${CMAKE_C_COMPILER_ID})
       	set(OpenMP_C "${CMAKE_C_COMPILER}")
        set(OpenMP_C_FLAGS
           "-Xpreprocessor -I/usr/local/opt/libomp/include -fopenmp=libomp -Wno-unused-command-line-argument")
        set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        message(STATUS "CXX compiler: " ${CMAKE_C_COMPILER_ID})
      	set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
      	set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
      	set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
      	set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
      	set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
      	set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    endif()
endif()

find_package(OpenMP REQUIRED)
if (APPLE AND OPENMP_FOUND)
    include_directories("${OPENMP_INCLUDES}")
    link_directories("${OPENMP_LIBRARIES}")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


# Verify include directories
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
	message(STATUS "\tinclude_dir='${dir}'")
endforeach()


add_subdirectory(bin)