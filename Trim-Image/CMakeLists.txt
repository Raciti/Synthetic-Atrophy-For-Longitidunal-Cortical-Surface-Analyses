cmake_minimum_required(VERSION 2.8)

project(trimImage)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(trimImage trimImage.cxx)

target_link_libraries(trimImage ${ITK_LIBRARIES})
