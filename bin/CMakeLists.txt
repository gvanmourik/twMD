# 
# bin/CMakeLists.txt
# 
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

add_executable(twMD main.cpp)
# use target_link_libraries to ink libraries

install(TARGETS twMD 
		DESTINATION ${PROJECT_BINARY_DIR}/bin
		PERMISSIONS WORLD_EXECUTE)