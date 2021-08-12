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
    Mesh() = default;
    explicit Mesh(const std::string& file_name);
    void local_laplacian(float lambda = 0.1);
    void global_laplacian();

    void update_vertices();

    MyMesh mesh_;
    IO::Options ropt_, wopt_;

    Eigen::MatrixXd V_;
    Eigen::MatrixXi F_;

private:
    void init_vertices();
};
