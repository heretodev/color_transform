cmake_minimum_required (VERSION 3.16)

project(color_transform VERSION 0.1.0) #Version vars are now autoset by project command.

#interface target for CImg
add_library(CImg INTERFACE)
target_include_directories(CImg INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/deps>)
target_link_options(CImg INTERFACE -pthread)
target_link_libraries(CImg INTERFACE X11)

add_executable(color_transform color_transform.cpp)
target_link_libraries(color_transform PRIVATE CImg)
