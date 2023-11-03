cmake_minimum_required(VERSION 3.10)
set(CMAKE_C_STANDARD 11)
project(test)
aux_source_directory(${PROJECT_SOURCE_DIR} SRC)
set(CMAKE_C_FLAGS "-O3 -march=native -Wall")

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/reference/)
include_directories(${PROJECT_SOURCE_DIR}/reference/libXKCP.a.headers)
link_libraries(XKCP)

add_executable(rsign ${SRC}
        verify.c
        verify.h)




