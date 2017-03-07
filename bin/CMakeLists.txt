cmake_minimum_required(VERSION 3.6)

project(bin)

set(CMAKE_CXX_STANDARD 14)

set(TEST_FILES test1.cpp)
add_executable(test1 ${TEST_FILES})

target_link_libraries(test1 lib)