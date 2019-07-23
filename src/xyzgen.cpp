#include <igl/read_triangle_mesh.h>
#include <igl/copyleft/cgal/CSGTree.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/jet.h>
#include <igl/writeOBJ.h>
#include <igl/writeDMAT.h>
#include <Eigen/Core>

int main(int argc, char * argv[])
{
    using namespace Eigen;
    using namespace igl::copyleft::cgal;
    using namespace std;
    using namespace igl;
    cout<<R"(
[,]  Toggle between boolean sub-tree operations
)";

    MatrixXi FA,FB,FC,FD,FE;
    MatrixXd VA,VB,VC,VD,VE;
    // Read in inputs as double precision floating point meshes
    read_triangle_mesh( "../data/cube.obj"     ,VA,FA);
    read_triangle_mesh("../data/sphere.obj"   ,VB,FB);
    read_triangle_mesh( "../data/xcylinder.obj",VC,FC);
    read_triangle_mesh( "../data/ycylinder.obj",VD,FD);
    read_triangle_mesh( "../data/zcylinder.obj",VE,FE);
    igl::opengl::glfw::Viewer viewer;

    int num_views = 5+4;
    int view_id = num_views-1;
    view_id = 7;
    const auto & update = [&]() {
        viewer.data().clear();
        // CSGTree templated on type of F
        VectorXd I;
        const auto & set_mesh =
                [&](const MatrixXd & V, const MatrixXi & F, const int i)
                {
                    viewer.data().set_mesh(V,F);
                    I = VectorXd::Constant(F.rows(),1,i);
                };
        CSGTree M;
        Matrix<MatrixXi::Index,Dynamic,1> J;

        // Compute result of (C ∪ D) ∪ E
        M = {{{VC,FC},{VD,FD},"u"},{VE,FE},"u"};
        J = M.J().array()+FA.rows()+FB.rows();
        igl::writeOBJ("../data/joint.obj", M.cast_V<MatrixXd>(), M.F() );

        viewer.data().set_mesh(M.cast_V<MatrixXd>(),M.F());
        I.resize(M.F().rows(),1);
        // Compute colors based on original facets
        for(int f = 0;f<M.F().rows();f++)
        {
            const int j = J(f);
            I(f) =
                    (int)(j<FA.rows())+
                    (int)(j<FA.rows()+FB.rows())+
                    (int)(j<FA.rows()+FB.rows()+FC.rows())+
                    (int)(j<FA.rows()+FB.rows()+FC.rows()+FD.rows())+
                    (int)(j<FA.rows()+FB.rows()+FC.rows()+FD.rows()+FE.rows());
        }


        MatrixXd C;
        jet(I,1,5,C);
        viewer.data().set_colors(C);
        igl::writeDMAT("../data/jointC.dmat", C);
    };
    update();
    viewer.launch();
}
