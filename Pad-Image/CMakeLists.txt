cmake_minimum_required(VERSION 2.8)

project(padImage)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(padImage padImage.cxx)

target_link_libraries(padImage ${ITK_LIBRARIES})
