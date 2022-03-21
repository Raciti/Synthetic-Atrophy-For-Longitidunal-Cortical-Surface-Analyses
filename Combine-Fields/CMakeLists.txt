cmake_minimum_required(VERSION 2.8)

project(combineFields)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(combineFields combineFields.cxx)

target_link_libraries(combineFields ${ITK_LIBRARIES})
