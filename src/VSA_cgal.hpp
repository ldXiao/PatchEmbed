//
// Created by Lind Xiao on 7/23/19.
//

#ifndef OTMAPPING_VSA_CGAL_H
#define OTMAPPING_VSA_CGAL_H
#include <iostream>
#include <fstream>
#include <map>
#include <tuple>
#include <nlohmann/json.hpp>
#include <Eigen/Core>
#include <igl/copyleft/cgal/mesh_to_polyhedron.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh_approximation/approximate_triangle_mesh.h>
#include <CGAL/boost/graph/copy_face_graph.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef Mesh::Property_map<face_descriptor, std::size_t> Face_proxy_pmap;
namespace VSA = CGAL::Surface_mesh_approximation;
namespace OTMapping {
    typedef std::tuple<int, int, int> face_tuple;

    int vsa_compute(
            const Eigen::MatrixXd &V,
            const Eigen::MatrixXi &F,
            int num_proxies,
            Eigen::MatrixXi & face_proxies,
            int & count) {
        Polyhedron plhd;
        Mesh mesh;
        // convert to polyhedron
        igl::copyleft::cgal::mesh_to_polyhedron(V, F, plhd);
        // The output will be an indexed triangle mesh

        // transfer the polyhedron into suface_mesh
        CGAL::copy_face_graph(plhd, mesh);

        // output face proxy index property map
        Face_proxy_pmap fpxmap =
                mesh.add_property_map<face_descriptor, std::size_t>("f:proxy_id", 0).first;
        // output planar proxies
        std::vector<Kernel::Vector_3> proxies;

        // output indexed triangle mesh
        std::vector<Kernel::Point_3> anchors;
        std::vector<CGAL::cpp11::array<std::size_t, 3> > triangles;
        // free function interface with named parameters
        VSA::approximate_triangle_mesh(mesh,
                                       CGAL::parameters::verbose_level(VSA::MAIN_STEPS). // seeding with minimum error drop
                                               max_number_of_proxies(num_proxies).
                                               face_proxy_map(fpxmap). // get face partition map
                                               proxies(std::back_inserter(proxies)). // output proxies
                                               anchors(std::back_inserter(anchors)). // output anchor points
                                               triangles(std::back_inserter(triangles))); // output indexed triangles
        std::cout << "#anchor points: " << anchors.size() << std::endl;
        std::cout << "#triangles: " << triangles.size() << std::endl;

        // iterates over faces and outputs segment id to console
        count  = 0;
        face_proxies.resize(F.rows(),1);
        BOOST_FOREACH(face_descriptor f, faces(mesh)) {
            face_proxies(count,0)= fpxmap[f];
            count +=1;
        }
        return EXIT_SUCCESS;
    }
}
#endif //OTMAPPING_VSA_CGAL_H
