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

# Boost setup
	# IF boost components need to be add use:
	# find_package( Boost REQUIRED COMPONENTS <component1> <component2> )
	# AND add the following line in the directory where the main executable is being created:
	# target_link_libraries( run ${Boost_<COMPONENT1>_LIBRARY} ${Boost_<COMPONENT2>_LIBRARY} )
find_package(Boost REQUIRED)
if(Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS})
else()
	message(STATUS "Boost NOT found...install Boost and try again!")
endif()

# Verify include directories
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
	message(STATUS "\tinclude_dir='${dir}'")
endforeach()


add_subdirectory(bin)