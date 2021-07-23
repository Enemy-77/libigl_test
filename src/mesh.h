#pragma once

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Utils/getopt.h>

#include <Eigen/core>

using namespace OpenMesh;
// ----------------------------------------------------------------------------
typedef TriMesh_ArrayKernelT<>  MyMesh;

class Mesh {
public:
    explicit Mesh(const std::string& file_name);
    void mean_curvature_flow(float lambda = 0.1);
    void update_vertices();

    MyMesh mesh_;
    IO::Options ropt_, wopt_;

    Eigen::MatrixXd V_;
    Eigen::MatrixXi F_;

private:
    void init_vertices();
};
