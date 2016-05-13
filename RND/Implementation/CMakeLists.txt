cmake_minimum_required (VERSION 2.8.11)
project (RND)

set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER cmake)

#Changing target directories of cmake outputs
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Eigen)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Eigen)


include_directories(${CMAKE_SOURCE_DIR})


#Compile project located in src
include_directories (src)
add_subdirectory (src)
