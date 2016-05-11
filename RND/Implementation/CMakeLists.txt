cmake_minimum_required (VERSION 2.8.11)
project (RND)

file(GLOB SOURCES
    src/*.h
	src/*.cpp
)

add_executable(RND ${SOURCES})
