cmake_minimum_required(VERSION 3.1)
project(pybcclean)

add_subdirectory(${CMAKE_SOURCE_DIR}/pybind11)
message(${CMAKE_SOURCE_DIR}/pybind11)
pybind11_add_module(pybcclean pybcclean/pybcclean.cpp)
