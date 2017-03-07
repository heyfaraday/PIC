cmake_minimum_required(VERSION 3.6)

project(lib)

set(SOURCE_LIB boris_push.cpp constants.cpp parameters.cpp field.cpp)
set(HEADER_LIB boris_push.hpp constants.hpp parameters.hpp field.hpp)

add_library(lib STATIC ${SOURCE_LIB} ${HEADER_LIB})