cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

set(CMAKE_SYSTEM_VERSION 10.0 CACHE STRING "" FORCE)
set(ROOT_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
set(KIDS_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Trilateral/Mesh/off")
set(TOSCA_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Trilateral/Mesh/SCB/Data/TOSCA/Meshes")
set(SCAPE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Trilateral/Mesh/SCB/Data/SCAPE/Meshes")
set(WATERTIGHT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Trilateral/Mesh/SCB/Data/watertight_shrec07/Meshes")
add_compile_definitions(ROOT_PATH="${ROOT_PATH}")
add_compile_definitions(KIDS_PATH="${KIDS_PATH}")
add_compile_definitions(TOSCA_PATH="${TOSCA_PATH}")
add_compile_definitions(SCAPE_PATH="${SCAPE_PATH}")
add_compile_definitions(WATERTIGHT_PATH="${WATERTIGHT_PATH}")

project(Trilateral CXX)

add_subdirectory(Trilateral)

