cmake_minimum_required(VERSION 3.7)
project(untitled)

set(CMAKE_CXX_STANDARD 14)

find_package(Eigen3)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}-Wall -O0 -ansi -pedantic -Wno-unused-but-set-variable -I$(SRC) -I.. -DTRACE -DENAGLO -DENAQZS -DENACMP -DENAGAL -DENAIRN -DNFREQ=3 -DSVR_REUSEADDR -g")

include_directories(${PROJECT_SOURCE_DIR})
aux_source_directory(. DIR_SRCS)

add_executable(untitled ${DIR_SRCS})

target_link_libraries(untitled m rt pthread)