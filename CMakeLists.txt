cmake_minimum_required(VERSION 3.1)
project(otmapping)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(DownloadProject)
download_project(
    PROJ         libigl
    SOURCE_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/external/libigl
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/.cache/libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG       v2.1.0)
# libigl
option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR       "Use AntTweakBar"    OFF)
option(LIBIGL_WITH_CGAL              "Use CGAL"           ON)
option(LIBIGL_WITH_COMISO            "Use CoMiso"         OFF)
option(LIBIGL_WITH_CORK              "Use Cork"           OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"         OFF)
option(LIBIGL_WITH_LIM               "Use LIM"            OFF)
option(LIBIGL_WITH_MATLAB            "Use Matlab"         OFF)
option(LIBIGL_WITH_MOSEK             "Use MOSEK"          OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          ON)
option(LIBIGL_WITH_PNG               "Use PNG"            OFF)
option(LIBIGL_WITH_PYTHON            "Use Python"         OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"       OFF)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_XML               "Use XML"            OFF)

find_package(LIBIGL REQUIRED QUIET)
download_project(
        PROJ         json
        SOURCE_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/external/json
        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/.cache/json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG        v3.6.1
)
#json
add_library(json INTERFACE)
target_include_directories(json SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/external/json/include)

download_project(
        PROJ         cxxopts
        SOURCE_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/external/cxxopts
        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/.cache/cxxopts
        GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
        GIT_TAG        v2.1.2)



#cxxopts
#add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/cxxopts)
add_library(cxxopts INTERFACE)
target_include_directories(cxxopts SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/external/cxxopts/include)

download_project(
        PROJ         tinyexpr
        SOURCE_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/external/tinyexpr
        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/.cache/tinyexpr
        GIT_REPOSITORY https://github.com/codeplea/tinyexpr.git
        GIT_TAG        master)
add_library(tinyexpr ${CMAKE_CURRENT_SOURCE_DIR}/external/tinyexpr/tinyexpr.c)
target_include_directories(tinyexpr SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/external/tinyexpr)

#
set(SOURCES
    src/CutMesh.cpp
    src/CutMesh.h
    src/Sinkhorn.hpp
    src/Skeleton.cpp
    src/Skeleton.h
    src/exportCutMesh.cpp
    src/Kruskal.hpp
    src/ExpressionValue.cpp
    src/ExpressionValue.hpp)
add_library(otmapping ${SOURCES})
target_link_libraries(otmapping PUBLIC
        json
        igl::core
        igl::cgal
        igl::opengl_glfw
        igl::opengl_glfw_imgui
        cxxopts
        tinyexpr)
# Add your project files

add_executable(${PROJECT_NAME}_bin src/main.cpp)
target_link_libraries(${PROJECT_NAME}_bin PUBLIC otmapping)
add_executable(cuttest_bin src/main2.cpp src/Sinkhorn.hpp)
target_link_libraries(cuttest_bin igl::core igl::cgal igl::opengl igl::opengl_glfw)
add_executable(cutmeshtest_bin
        src/CutMeshTest.cpp
        )
target_link_libraries(cutmeshtest_bin otmapping tinyexpr)
add_executable(sinkhorntest_bin
        src/Sinkhorn.hpp
        src/sinkhorntest.cpp
        src/ExpressionValue.cpp)
target_link_libraries(sinkhorntest_bin igl::core json tinyexpr)
add_executable(quadratic_OT_bin
        src/quadratic_OT.hpp
        src/quadratic_OT_test.cpp)
target_link_libraries(quadratic_OT_bin igl::core)
add_executable(exportCutMesh_bin
        src/exportCutMesh.cpp)
target_link_libraries(exportCutMesh_bin otmapping)

