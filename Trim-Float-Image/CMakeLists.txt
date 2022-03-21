cmake_minimum_required(VERSION 2.8)

project(trimFloatImage)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(trimFloatImage trimFloatImage.cxx)

target_link_libraries(trimFloatImage ${ITK_LIBRARIES})
