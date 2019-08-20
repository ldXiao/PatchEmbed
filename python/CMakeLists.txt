cmake_minimum_required(VERSION 3.1)
project(pybcclean)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_CGAL            "Use CGAL"  ON)

find_package(LIBIGL REQUIRED QUIET)


add_subdirectory(${CMAKE_SOURCE_DIR}/external/pybind11)

pybind11_add_module(pybcclean pybcclean/pybcclean.cpp)
