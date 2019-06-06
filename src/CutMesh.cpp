//
// Created by Lind Xiao on 6/6/19.
//

#include "CutMesh.h"
#include <igl/random_points_on_mesh.h>
#include <igl/opengl/glfw/Viewer.h>

void CutMesh::PlotCutMesh(igl::opengl::glfw::Viewer &viewer, std::string options){
    viewer.data().clear();
    viewer.data().set_mesh(this->Vertices, this->Faces);
}