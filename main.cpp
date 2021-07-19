#include <igl/opengl/glfw/Viewer.h>

// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Utils/getopt.h>

#include <iostream>

#include <Eigen/core>

using namespace OpenMesh;
// ----------------------------------------------------------------------------
typedef TriMesh_ArrayKernelT<>  MyMesh;

int main(int argc, char *argv[])
{
    MyMesh mesh;
    IO::Options ropt, wopt;
    if (!read_mesh(mesh, "assets/models/Nefertiti_face.obj", ropt)) {
        std::cerr << "Error loading file " << std::endl;
    }
    int num_vertices = mesh.n_vertices();
    int num_faces = mesh.n_faces();
    Eigen::MatrixXd V(num_vertices, 3);
    Eigen::MatrixXi F(num_faces, 3);
    int row = 0;
    for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
        MyMesh::Point p = mesh.point(v_it);
        V(row, 0) = p[0];
        V(row, 1) = p[1];
        V(row, 2) = p[2];
        row++;
    }
    row = 0;
    for (MyMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
        int col = 0;
        for (MyMesh::FaceVertexIter fvi = mesh.fv_begin(*f_it); fvi.is_valid(); ++fvi) {
            MyMesh::Point p = mesh.point(*fvi);
            F(row, col) = fvi->idx();
            col++;
        }
        row++;
    }

  // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    viewer.data().set_face_based(true);
    viewer.launch();
}
