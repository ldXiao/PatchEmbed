cmake_minimum_required(VERSION 3.1)
project(bcclean)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/)
set_property(GLOBAL PROPERTY ALLOW_DUPLICATE_CUSTOM_TARGETS 1)
# include_directories(${PYTHON_INCLUDE_DIRS})
# find_package(PythonLibs 3.6)
# message(${PYTHON_LIBRARIES})
# find_package(PythonInterp 3.6)
# message(${PYTHON_INCLUDE_DIRS})


include(DownloadProject)
# download_project(
    # PROJ         libigl
    # SOURCE_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/external/libigl
    # DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/.cache/libigl
    # GIT_REPOSITORY https://github.com/libigl/libigl.git
    # GIT_TAG       master)

option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR       "Use AntTweakBar"    OFF)
option(LIBIGL_WITH_CGAL              "Use CGAL"           ON)
option(LIBIGL_WITH_COMISO            "Use CoMiso"         OFF)
option(LIBIGL_WITH_CORK              "Use Cork"           OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"         ON)
option(LIBIGL_WITH_LIM               "Use LIM"            OFF)
option(LIBIGL_WITH_MATLAB            "Use Matlab"         OFF)
option(LIBIGL_WITH_MOSEK             "Use MOSEK"          OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         OFF)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           OFF)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          OFF)
option(LIBIGL_WITH_PNG               "Use PNG"            OFF)
option(LIBIGL_WITH_PYTHON            "Use Python"         OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"       OFF)
option(LIBIGL_WITH_PREDICATES        "Use exact predicates" ON)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  OFF)
option(LIBIGL_WITH_XML               "Use XML"            OFF)
# find_package(LIBIGL REQUIRED QUIET)

# # pybind11
# download_project(
#         PROJ         pybind11
#         SOURCE_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/external/pybind11
#         DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/.cache/pybind11
#         GIT_REPOSITORY https://github.com/pybind/pybind11.git
#         GIT_TAG        v2.3.0
# )
# set(PYBIND11_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/pybind11)
# add_subdirectory(${PYBIND11_DIR}/ pybind11)

download_project(
PROJ         nanoflann
SOURCE_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/external/spdlog
DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/.cache/spdlog
GIT_REPOSITORY https://github.com/gabime/spdlog.git
GIT_TAG        v1.5.0
)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/spdlog)

# nanoflann

download_project(
PROJ         nanoflann
SOURCE_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/external/nanoflann
DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/.cache/nanoflann
GIT_REPOSITORY https://github.com/jlblancoc/nanoflann.git
GIT_TAG        v1.3.0
)
add_library(nanoflann INTERFACE)
target_include_directories(nanoflann SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/external/nanoflann/include)
## tetwild
#download_project(
#        PROJ         tetwild
#        SOURCE_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/external/tetwild
#        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/.cache/tetwild
#        GIT_REPOSITORY https://github.com/Yixin-Hu/TetWild.git
#)
#add_library(tetwild ${CMAKE_CURRENT_SOURCE_DIR}/external/tetwild/src/tetwild)
#target_include_directories(tetwild SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/external/tetwild/include)

download_project(
        PROJ         FloatTetwild
        SOURCE_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/external/fTetwild
        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/.cache/fTetwild
        GIT_REPOSITORY https://github.com/wildmeshing/fTetWild.git
        GIT_TAG        2308d2f
)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/fTetwild)
# option(FLOAT_TETWILD_ENABLE_TBB         "Enable TBB"                    OFF)
#json 


# download_project(
#         PROJ         json
#         SOURCE_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/external/json
#         DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/.cache/json
#         GIT_REPOSITORY https://github.com/nlohmann/json.git
#         GIT_TAG        v3.6.1
# )

# #json
# # add_library(json INTERFACE)
# target_include_directories(json SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/external/json/include)

download_project(
        PROJ         cxxopts
        SOURCE_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/external/cxxopts
        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/.cache/cxxopts
        GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
        GIT_TAG        v2.2.0)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/cxxopts)
#cxxopts

# # yaml
download_project(
        PROJ         yaml-cpp
        SOURCE_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/external/yaml-cpp
        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/.cache/yaml-cpp
        GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
        GIT_TAG        master)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/yaml-cpp)


set(SOURCES
    src/graphcut_cgal.h
    src/bcclean.cpp
    src/bcclean.h
    src/kdtree_NN_Eigen.hpp
    src/node.h
    src/node.cpp
    src/edge.h
    src/edge.cpp    
    src/patch.h 
    src/patch.cpp
    src/planar_cut_simply_connect.h 
    src/planar_cut_simply_connect.cpp 
    src/patch_cut_relabel.h 
    src/patch_cut_relabel.cpp 
    src/proj_node.h
    src/proj_node.cpp
    src/loop_colorize.h
    src/loop_colorize.cpp
    src/Match_Maker_Tree.h 
    src/Match_Maker_Tree.cpp 
    src/Kruskal.h
    src/Kruskal.cpp
    src/fTetwild.h 
    src/fTetwild.cpp
    src/polyline_distance.h 
    src/polyline_distance.cpp 
    src/params.h
    src/orientation_check.h 
    src/orientation_check.cpp 
    src/Edge_Dijkstra.h 
    src/Edge_Dijkstra.cpp
    src/non_vertex_manifold_relabel.h 
    src/non_vertex_manifold_relabel.cpp 
    src/backtrack_diff.h 
    src/backtrack_diff.cpp
    src/recursive_cut.h 
    src/recursive_cut.cpp
    src/CellularGraph.h
    src/CellularGraph.cpp
    src/TraceComplex.h
    src/TraceComplex.cpp
    src/MatchMakerDynamic.h 
    src/MatchMakerDynamic.cpp
)
add_library(bcclean ${SOURCES})
target_link_libraries(bcclean PUBLIC
        igl::core
        igl::cgal
        igl::embree
        cxxopts
        yaml-cpp
        nanoflann
        FloatTetwild
        )

add_executable(btcmm src/main.cpp)
target_link_libraries(btcmm bcclean)
set_target_properties(btcmm PROPERTIES CXX_STANDARD 17)

add_executable(gen_feat src/gen_feat.cpp)
target_link_libraries(gen_feat bcclean)
set_target_properties(gen_feat PROPERTIES CXX_STANDARD 17)