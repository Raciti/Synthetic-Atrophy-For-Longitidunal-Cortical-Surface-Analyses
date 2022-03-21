cmake_minimum_required(VERSION 2.8)

project(padVectorImage)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(padVectorImage padVectorImage.cxx)

target_link_libraries(padVectorImage ${ITK_LIBRARIES})
