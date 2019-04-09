# 
# bin/CMakeLists.txt
# 
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

add_executable(twMD main.cpp)
target_link_libraries(twMD ${MPI_LIBRARIES})
TARGET_LINK_LIBRARIES( twMD ${Boost_LIBRARIES} )

install(TARGETS twMD 
		DESTINATION ${PROJECT_BINARY_DIR}/bin
		PERMISSIONS WORLD_EXECUTE)