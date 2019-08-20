cmake_minimum_required(VERSION 3.1)
project(pybcclean)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_CGAL              "Use CGAL"           ON)
option(LIBIGL_WITH_EMBREE            "Use Embree"         ON)
find_package(LIBIGL REQUIRED QUIET)

#nanoflann
add_library(nanoflann INTERFACE)
target_include_directories(nanoflann SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/external/nanoflann/include)

#pybind11
add_subdirectory(${CMAKE_SOURCE_DIR}/external/pybind11)

#bcclean_lib
set(BCCLEAN_SOURCES
    pybcclean/graphcut_cgal.h
    pybcclean/bcclean.cpp
    pybcclean/bcclean.h
    pybcclean/kdtree_NN_Eigen.hpp
)
add_library(bcclean_lib ${BCCLEAN_SOURCES})
target_link_libraries(bcclean_lib PUBLIC igl::core igl::cgal igl::embree nanoflann)


#binding modules
pybind11_add_module(pybcclean MODULE pybcclean/pybcclean.cpp)
target_link_libraries(pybcclean PRIVATE pybind11::module bcclean_lib)
set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                        SUFFIX "${PYTHON_MODULE_EXTENSION}")
