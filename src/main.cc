#include <igl/opengl/glfw/Viewer.h>
#include <igl/gaussian_curvature.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>

// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Utils/getopt.h>

#include <iostream>

#include <Eigen/core>

#include "mesh.h"

using namespace OpenMesh;
// ----------------------------------------------------------------------------
typedef TriMesh_ArrayKernelT<>  MyMesh;

Mesh mesh;

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) {
    std::cout << "Key: " << key << " " << (unsigned int)key << std::endl;
    if (key == '1') {
        for (int i = 0; i < 100; ++i) {
            mesh.local_laplacian(0.01);
        }
        mesh.update_vertices();
        viewer.data().clear();
        viewer.data().set_mesh(mesh.V_, mesh.F_);
        viewer.data().set_face_based(true);
    }
    return true;
}

int main() {
    mesh = Mesh("../assets/models/Balls.obj");

    igl::opengl::glfw::Viewer viewer;
//    viewer.callback_key_down = &key_down;

    mesh.global_laplacian();
    mesh.update_vertices();


    viewer.data().set_mesh(mesh.V_, mesh.F_);
    viewer.data().set_face_based(true);
    viewer.launch();
}
