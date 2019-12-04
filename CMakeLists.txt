cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/ext/glm/")
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
project(Stanford3dRepoConverter)

add_executable(${PROJECT_NAME} "main.cpp" "tinyply.cpp" "tinyply.h" "range_scan_reader.h" "range_scan_reader.cpp")
