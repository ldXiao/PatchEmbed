//
// Copyright (C) 2019 Lind Xiao on 4/4/19 <lx471@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/local_basis.h>
//#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/jet.h>
#include <igl/grad.h>
#include <fstream>
#include "hard_soft_nrosy.h"
#include "scalar_from_vector.h"

#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/lscm.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>

// Mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;

// 2d mesh for uv plane
Eigen::MatrixXd V_uv;

// 3D mesh for the edited uv plane, pad with zero
Eigen::MatrixXd V_uv_edit;

//indices to store boundary indices
Eigen::VectorXi bnd, bnd_head(2,1);
Eigen::MatrixXd bnd_head_uv(2,2);

// Triangle-triangle adjacency
Eigen::MatrixXi TT;
Eigen::MatrixXi TTi;

// Constrained faces id
Eigen::VectorXi b;

// Cosntrained faces representative vector
Eigen::MatrixXd bc;

// file_name for many ASCII dump
std::string file_body;

int uv_scale=10;

//control options
bool show_uv=false;
bool harmonic=false;
bool poisson_err = false;

// Currently selected face
int selected;

// Degree of the N-RoSy field
int N;

bool save_flips(std::string, std::vector<int>);

bool find_boundary(
        const Eigen::MatrixXi & F,
        Eigen::VectorXi & bnd){
    try {
        igl::boundary_loop(F, bnd);
    }
    catch(...) {
        std::cerr << "Please check that the model imported indeed has a boundary" << std::endl;
        return false;
    }
    return true;
}


// Local basis
Eigen::MatrixXd B1, B2, B3;

// Texture image (grayscale)
Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> texture_I;

void line_texture() {
    int size = 128;              // Texture size
    int w = 7;                // Line width
    int pos = size / 2 - w / 2; // Center the line
    texture_I.setConstant(size, size, 255);
    texture_I.block(0, pos, size, w).setZero();
    texture_I.block(pos, 0, w, size).setZero();

}

// Converts a representative vector per face in the full set of vectors that describe
// an N-RoSy field
void representative_to_nrosy(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        const Eigen::MatrixXd &R,
        const int N,
        Eigen::MatrixXd &Y) {
    using namespace Eigen;
    using namespace std;

    Y.resize(F.rows() * N, 3);
    for (unsigned i = 0; i < F.rows(); ++i) {
        double x = R.row(i) * B1.row(i).transpose();
        double y = R.row(i) * B2.row(i).transpose();
        double angle = atan2(y, x);

        for (unsigned j = 0; j < N; ++j) {
            double anglej = angle + 2 * M_PI * double(j) / double(N);
            double xj = cos(anglej);
            double yj = sin(anglej);
            Y.row(i * N + j) = xj * B1.row(i) + yj * B2.row(i);
            Y.row(i * N + j) = Y.row(i * N + j) * R.row(i).norm();
        }

    }
}

// Plots the mesh with an N-RoSy field
// The constrained faces (b) are colored in red.
void plot_mesh_nrosy(
        igl::opengl::glfw::Viewer &viewer,
        Eigen::MatrixXd &V,
        Eigen::MatrixXi &F,
        int N,
        Eigen::MatrixXd &PD1,
        Eigen::VectorXi &b) {
    using namespace Eigen;
    using namespace std;
    // Clear the mesh
    viewer.data().clear();
    viewer.data().set_mesh(V, F);
    viewer.data().set_texture(texture_I, texture_I, texture_I);

    // Expand the representative vectors in the full vector set and plot them as lines
    double avg = igl::avg_edge_length(V, F);
    MatrixXd Y;
    representative_to_nrosy(V, F, PD1, N, Y);

    MatrixXd B;
    igl::barycenter(V, F, B);

    MatrixXd Be(B.rows() * N, 3);
    for (unsigned i = 0; i < B.rows(); ++i)
        for (unsigned j = 0; j < N; ++j)
            Be.row(i * N + j) = B.row(i);

    viewer.data().add_edges(Be, Be + Y * (avg / 2), RowVector3d(0, 0, 1));

    // Highlight in red the constrained faces
    MatrixXd C = MatrixXd::Constant(F.rows(), 3, 1);
    for (unsigned i = 0; i < b.size(); ++i)
        C.row(b(i)) << 1, 0, 0;
    viewer.data().set_colors(C);
}

std::vector<int> detect_flip(
        const Eigen::MatrixXd & V_uv,
        const Eigen::MatrixXd & V_uv_edit,
        const Eigen::MatrixXi& F){
    std::vector<int> flip_idx;
    try{
        //calculate local basis
        Eigen::MatrixXd V_uv3d(V_uv.rows(), 3);
        Eigen::MatrixXd V_uv_edit3d(V_uv_edit.rows(), 3);
        for(int i =0; i< V_uv.rows(); ++i){
            V_uv3d(i,0) = V_uv(i,0);
            V_uv3d(i,1) = V_uv(i,1);
            V_uv3d(i,2) = 0;
            V_uv_edit3d(i,0) = V_uv_edit(i,0);
            V_uv_edit3d(i,1) = V_uv_edit(i,1);
            V_uv_edit3d(i,2) = 0;

        }

        Eigen::MatrixXd c1, c2, c3;
//        igl::local_basis(V_uv3d,F, b1,b2, b3);
        igl::local_basis(V_uv_edit3d, F, c1, c2, c3);
        double sum = c3.rowwise().sum()(2);
        for(int i =0; i < F.rows(); ++i){
//            double prod = b3.row(i).dot( c3.row(i));
            if( sum * c3(i,2)<0){
                flip_idx.push_back(i);
            }
        }
    }
    catch(...){std::cerr<<"flip detection failed"<<std::endl;}
    return flip_idx;

}

// It allows to change the degree of the field when a number is pressed
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier) {
    using namespace Eigen;
    using namespace std;
    if (key >= '1' && key <= '9') {

        if (key == '1') {
            // plot the plain mesh
            // Simple constraints
            b.resize(2);
            b(0) = 0;
            b(1) = F.rows() - 1;
            bc.resize(2, 3);
            bc << 1, 1, 1, 0, 1, 1;

            selected = 0;



            viewer.data().clear();
            viewer.data().set_mesh(V, F);
            viewer.data().set_texture(texture_I, texture_I, texture_I);
            viewer.data().show_texture = false;
            // clear the global variables
        }
        if (key == '2') {
            // plot the nrosy field with degree N
            MatrixXd R = Hard_nrosy(V, F, TT, b, bc, N);
            plot_mesh_nrosy(viewer, V, F, N, R, b);
            viewer.data().show_texture = false;
        }

        if(key == '3') {
            if(!poisson_err) {
                viewer.data().clear();
                viewer.data().set_mesh(V, F);
                viewer.data().set_texture(texture_I, texture_I, texture_I);
                // plot scalar function and gradient field
                MatrixXd R = Hard_nrosy(V, F, TT, b, bc, N);
                plot_mesh_nrosy(viewer, V, F, N, R, b);
                Eigen::MatrixXd s = scalar_from_vector(V, F, TT, R);
                igl::writeDMAT(file_body + ".dmat", s);
                // Compute pseudocolor for original function
                Eigen::MatrixXd C;
                igl::jet(s, true, C);
                viewer.data().set_colors(C);
                Eigen::MatrixXd BC;
                igl::barycenter(V, F, BC);
                Eigen::MatrixXd Gs = gradient_scalar(V, F, s);
                double avg = igl::avg_edge_length(V, F);
                const RowVector3d red(1, 0, 0);
                viewer.data().add_edges(BC, BC + Gs * (avg / 2), red);
                viewer.data().show_texture = false;
            }
            else{
                viewer.data().clear();
                viewer.data().set_mesh(V, F);
                viewer.data().set_texture(texture_I, texture_I, texture_I);
                // plot scalar function and gradient field
                MatrixXd R = Hard_nrosy(V, F, TT, b, bc, N);
                plot_mesh_nrosy(viewer, V, F, N, R, b);
                Eigen::MatrixXd s = scalar_from_vector(V, F, TT, R);
                // Compute pseudocolor for original function
                Eigen::MatrixXd BC;
                igl::barycenter(V, F, BC);
                Eigen::MatrixXd Gs = gradient_scalar(V, F, s);
                Eigen::MatrixXd poisson(F.rows(), 1);
                for(int i =0; i< F.rows(); ++i){
                    Eigen::VectorXd d = Gs.row(i) - R.row(i);
                    poisson(i,0)= d.norm();
                }
                Eigen::MatrixXd C;
                igl::jet(poisson, true, C);
                viewer.data().set_colors(C);
                double avg = igl::avg_edge_length(V, F);
                const RowVector3d red(1, 0, 0);
                viewer.data().add_edges(BC, BC + Gs * (avg / 2), red);
                viewer.data().show_texture = false;

            }

        }
    }

    if (key == 'Q' || key == 'W') {
        if (selected >= b.size() || selected < 0)
            return false;

        bc.row(selected) = bc.row(selected) * (key == 'Q' ? 3. / 2. : 2. / 3.);

        MatrixXd R = Hard_nrosy(V, F, TT, b, bc, N);
        plot_mesh_nrosy(viewer, V, F, N, R, b);
    }
    return false;
}

bool mouse_down(igl::opengl::glfw::Viewer &viewer, int, int) {
    int fid_ray;
    Eigen::Vector3f bary;
    // Cast a ray in the view direction starting from the mouse position
    double x = viewer.current_mouse_x;
    double y = viewer.core.viewport(3) - viewer.current_mouse_y;
    if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view,
                                 viewer.core.proj, viewer.core.viewport, V, F, fid_ray, bary)) {
        bool found = false;
        for (int i = 0; i < b.size(); ++i) {
            if (b(i) == fid_ray) {
                found = true;
                selected = i;
            }
        }

        if (!found) {
            b.conservativeResize(b.size() + 1);
            b(b.size() - 1) = fid_ray;
            bc.conservativeResize(bc.rows() + 1, bc.cols());
            bc.row(bc.rows() - 1) << 1, 1, 1;
            selected = bc.rows() - 1;

            Eigen::MatrixXd R = Hard_nrosy(V, F, TT, b, bc, N);
            plot_mesh_nrosy(viewer, V, F, N, R, b);
        }

        return true;
    }
    return false;
};


bool read_constraints(std::string file_name, Eigen::VectorXi & b, Eigen::MatrixXd & bc){
    using namespace std;
    fstream file;
    file.open(file_name);
    if (file.is_open()){
        int max_size = F.rows();
        b.resize(max_size);
        bc.resize(max_size, 3);
        string line;
        int count =0;
        while (getline(file, line)) {

            stringstream line_stream(line);
            // read the first integer into b
            int temp_b;
            line_stream >> temp_b;
            b[count] = temp_b;
            // read the remaining 3 doubles into bc
            for(int i =0; i< 3; ++ i){
                double temp_bc;
                line_stream >> temp_bc;
                bc.row(count)[i] = temp_bc;
            }
            count += 1;
        }
        b.conservativeResize(count);
        bc.conservativeResize(count, 3);
        return true;
    }
    else {
        cerr << "Can not open the file path " + file_name << "please check and retry"<<endl;
        return false;
    }
}

bool save_interpolations(std::string file_name, Eigen::VectorXi b, Eigen::MatrixXd bc){
    using namespace std;
    fstream file;
    file.open(file_name, fstream::out);
    int n = bc.rows();
    file << " Input constraints" << '\n';
    for(int i = 0; i < n; ++i){
        file << b[i] << ' ';
        file << bc.row(i)[0] << ' ';
        file << bc.row(i)[1] << ' ';
        file << bc.row(i)[2] << '\n';
    }
    file << '\n';

    Eigen::MatrixXd R = Hard_nrosy(V, F, TT, b, bc, N);
    // Compute pseudocolor for original function
    file << " Hard nrosy results" << '\n';
    for(int i = 0; i < R.rows(); ++i){
        file << i << ' ';
        file << R.row(i)[0] << ' ';
        file << R.row(i)[1] << ' ';
        file << R.row(i)[2] << '\n';
    }
    file.close();
    return true;
}

bool save_scalar(std::string file_name, Eigen::MatrixXd s, Eigen::MatrixXd s1){
    using namespace std;
    fstream file;
    file.open(file_name, fstream::out);
    file << " index input output" << '\n';
    for(int i =0 ; i < V.rows() ; ++i){
        file << i << "     ";
        file << s(i,0) << "     ";
        file << s1(i,0) << '\n';
    }
    file.close();
    return true;
}

bool save_flips(std::string file_name, std::vector<int> flip_idx){
    using namespace std;
    fstream file;
    file.open(file_name, fstream::out);
    file << " flips" << '\n';
    for(int i =0 ; i < flip_idx.size() ; ++i){
        file << flip_idx[i] << '\n';
    }
    file.close();
    return true;
}

bool save_scalar(std::string file_name, Eigen::MatrixXd s){
    using namespace std;
    fstream file;
    file.open(file_name, fstream::out);
    file << " index output" << '\n';
    for(int i =0 ; i < V.rows() ; ++i){
        file << i << "     ";
        file << s(i,0) << '\n';
    }
    file.close();
    return true;
}

bool save_test_scalar(std::string file_name, const Eigen::MatrixXd& V){
    Eigen::MatrixXd s(V.rows(),1);
    for(int i = 0; i < V.rows(); ++i){
        s(i,0) = V(i,1);
    }
    igl::writeDMAT(file_name, s);
    return true;
}

int main(int argc, char *argv[]) {
    using namespace std;
    using namespace Eigen;
    if (argc != 2) {
        cout << "Usage ex2_bin mesh.obj" << endl;
        exit(0);
    }

    string input = argv[1];
    int len = input.length();
    file_body = input.substr(0,len - 4);




    N = 1;
    // Load a mesh in OBJ format
    igl::readOBJ(argv[1], V, F);
    line_texture();
    // Triangle-triangle adjacency
    igl::triangle_triangle_adjacency(F, TT, TTi);
    // Compute the local_basis
    igl::local_basis(V, F, B1, B2, B3);

    // Simple constraints
    b.resize(2);
    b(0) = 0;
    b(1) = F.rows() - 1;
    bc.resize(2, 3);
    bc << 1, 1, 1, 0, 1, 1;

    selected = 0;

    igl::opengl::glfw::Viewer viewer;

    // Interpolate the field and plot
    key_down(viewer, '1', 0);

    // Plot the mesh
    viewer.data().set_mesh(V, F);
    viewer.data().set_texture(texture_I, texture_I, texture_I);

    // Register the callbacks
    viewer.callback_key_down = &key_down;
    viewer.callback_mouse_down = &mouse_down;

    // Disable wireframe
    viewer.data().show_lines = false;


//    viewer.data().show_texture = true;

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);

    menu.callback_draw_viewer_menu = [&]() {
        // Add widgets to the sidebar.
        if (ImGui::CollapsingHeader("Reconstruction Options", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::InputScalar("uv_scale", ImGuiDataType_U32, &uv_scale);
//            if (ImGui::Button("Reset Fields", ImVec2(-1, 0))) {
//                // Recreate the grid
//                key_down(viewer, '2', 0);
//            }
            if (ImGui::Button("Load Constraints", ImVec2(-2, 0))) {
                // Load the constraint from file
                viewer.data().clear();

                string file_name = file_body + ".constraints";

                read_constraints(file_name, b, bc);

                key_down(viewer, '2', 0);

            }
            if (ImGui::Button("Save Interpolation", ImVec2(-3, 0))) {
                // Save the constraints
                string file_name = file_body + ".results";
                save_interpolations(file_name, b, bc);
            }
            if (ImGui::Button("Load Test Scalar", ImVec2(-4, 0))) {
                // Load the constraint from file
//                viewer.data().clear();
                key_down(viewer, '1', 0);
                string file_name = file_body+ "test_scalar" + ".dmat";
                Eigen::MatrixXd s;
                save_test_scalar(file_name, V);

//                read_constraints(file_name, b, bc);
                igl::readDMAT(file_name, s);
//                key_down(viewer, '2', 0);

                Eigen::MatrixXd BC;
                igl::barycenter(V,F,BC);
                Eigen::MatrixXd Gs = gradient_scalar(V, F, s);
                double avg = igl::avg_edge_length(V, F);
                const RowVector3d blue(0,0,1);
                viewer.data().add_edges(BC,BC+ Gs * (avg / 2), blue);

                Eigen::MatrixXd s1 = scalar_from_vector(V, F, TT,  Gs);

                igl::writeDMAT(file_body+".dmat", s1);
//                cout << s1 - s << endl;
                // Compute pseudocolor for original function
                Eigen::MatrixXd Gs1 = gradient_scalar(V, F, s1);
                Eigen::MatrixXd C;
                igl::jet(s1, true,C);
                viewer.data().set_colors(C);
                const RowVector3d red(1,0,0);
                viewer.data().add_edges(BC,BC+ Gs1 * (avg / 2), red);
                save_scalar(file_body+".scalars", s, s1);
            }
            if (ImGui::Button("Toggle show_uv", ImVec2(-5, 0))) {
                show_uv = !show_uv;
            }
            if (ImGui::Button("Toggle poisson_err", ImVec2(-6, 0))) {
                poisson_err = !poisson_err;
            }
        }

    };


    // Launch the viewer
    viewer.launch();
}
