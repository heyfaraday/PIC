cmake_minimum_required(VERSION 3.6)

project(PIC)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES main.cpp)
add_executable(main ${SOURCE_FILES})

include_directories(lib)
add_subdirectory(lib)
target_link_libraries(main lib)

include_directories(bin)
add_subdirectory(bin)