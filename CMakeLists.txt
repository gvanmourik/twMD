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

# Set the clang compiler
set(CMAKE_C_COMPILER "/usr/local/Cellar/llvm/8.0.0/bin/clang")
set(CMAKE_CXX_COMPILER "/usr/local/Cellar/llvm/8.0.0/bin/clang++")
set(OPENMP_LIBRARIES "/usr/local/Cellar/llvm/8.0.0/lib")
set(OPENMP_INCLUDES "/usr/local/Cellar/llvm/8.0.0/include")


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

set(CMAKE_CXX_STANDARD 11)  # enable C++11 standard

# Boost setup
	# IF boost components need to be add use:
	# find_package( Boost REQUIRED COMPONENTS <component1> <component2> )
	# AND add the following line in the directory where the main executable is being created:
	# target_link_libraries( run ${Boost_<COMPONENT1>_LIBRARY} ${Boost_<COMPONENT2>_LIBRARY} )
find_package(Boost REQUIRED)
if(Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS})
endif()

# OpenMP check
# if(APPLE)
    if(CMAKE_C_COMPILER_ID MATCHES "Clang")
       	set(OpenMP_C "${CMAKE_C_COMPILER}")
        set(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
        set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      	set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
      	set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
      	set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
      	set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
      	set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
      	set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    endif()
# endif()

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    include_directories("${OPENMP_INCLUDES}")
    link_directories("${OPENMP_LIBRARIES}")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif(OPENMP_FOUND)



# Verify include directories
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
	message(STATUS "\tinclude_dir='${dir}'")
endforeach()


add_subdirectory(bin)