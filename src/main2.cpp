//
// Created by Lind Xiao on 6/6/19.
//

#include <igl/readOBJ.h>
//#undef IGL_STATIC_LIBRARY
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <igl/copyleft/cgal/intersect_with_half_space.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/slice_mask.h>
#include <Eigen/Core>
#include <iostream>

std::string DATA_PATH = "../data/";

Eigen::MatrixXd VA,VB,VC;
Eigen::VectorXi J,I;
Eigen::MatrixXi FA,FB,FC;
igl::MeshBooleanType boolean_type(
        igl::MESH_BOOLEAN_TYPE_UNION);


const char * MESH_BOOLEAN_TYPE_NAMES[] =
        {
                "Union",
                "Intersect",
                "Minus",
                "XOR",
                "Resolve",
        };

void update(igl::opengl::glfw::Viewer &viewer)
{
    igl::copyleft::cgal::mesh_boolean(VA,FA,VB,FB,boolean_type,VC,FC,J);
    Eigen::MatrixXd VE;
    Eigen::MatrixXi  FE;
    {
        Eigen::Vector3d p;
        p << 0,0,0.5;
        Eigen::Vector3d n;
        n << 0,0,1;
        std::cout << n<<std::endl;
        Eigen::VectorXi J;
        igl::copyleft::cgal::intersect_with_half_space(VC, FC, p, n, VE, FE,J);
    }

    viewer.data().clear();
    viewer.data().set_mesh(VE,FE);
    Eigen::MatrixXd C(FC.rows(),3);
    Eigen::MatrixXd VD;
    Eigen::MatrixXi  FD(FC.rows(),3);
    int count = 0;
    for(size_t f = 0;f<C.rows();f++)
    {
        if(J(f)<FA.rows())
        {
            FD.row(count) = FC.row(f);
            C.row(f) = Eigen::RowVector3d(1,0,0);
            count+=1;
        }else
        {
            C.row(f) = Eigen::RowVector3d(0,1,0);
        }
    }
    FD.conservativeResize(count, 3);
    {
        Eigen::VectorXi I;
        igl::remove_unreferenced(VC, Eigen::MatrixXi(FD), VD, FD, I);
    }
//    viewer.data().set_colors(C);
    std::cout<<"A "<<MESH_BOOLEAN_TYPE_NAMES[boolean_type]<<" B."<<std::endl;
    }

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods)
{
    switch(key)
    {
        default:
            return false;
        case '.':
            boolean_type =
                    static_cast<igl::MeshBooleanType>(
                            (boolean_type+1)% igl::NUM_MESH_BOOLEAN_TYPES);
            break;
        case ',':
            boolean_type =
                    static_cast<igl::MeshBooleanType>(
                            (boolean_type+igl::NUM_MESH_BOOLEAN_TYPES-1)%
                            igl::NUM_MESH_BOOLEAN_TYPES);
            break;
        case '[':
            viewer.core.camera_dnear -= 0.1;
            return true;
        case ']':
            viewer.core.camera_dnear += 0.1;
            return true;

    }
    update(viewer);
    return true;
}

int main(int argc, char *argv[])
{
    using namespace Eigen;
    using namespace std;
    igl::readOBJ(DATA_PATH+"/cube.obj",VA,FA);
    igl::readOBJ(DATA_PATH+"/polysphere.obj",VB,FB);
    // Plot the mesh with pseudocolors
    igl::opengl::glfw::Viewer viewer;

    // Initialize
    update(viewer);

    viewer.data().show_lines = true;
    viewer.callback_key_down = &key_down;
    viewer.core.camera_dnear = 3.9;
    cout<<
        "Press '.' to switch to next boolean operation type."<<endl<<
        "Press ',' to switch to previous boolean operation type."<<endl<<
        "Press ']' to push near cutting plane away from camera."<<endl<<
        "Press '[' to pull near cutting plane closer to camera."<<endl<<
        "Hint: investigate _inside_ the model to see orientation changes."<<endl;
    viewer.launch();
}

