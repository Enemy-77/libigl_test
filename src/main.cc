#include <igl/opengl/glfw/Viewer.h>

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

int main() {
    Mesh mesh("../assets/models/Nefertiti_face.obj");

    for (int i = 0; i < 50; ++i) {
        mesh.mean_curvature_flow();
    }

    

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(mesh.V_, mesh.F_);
    viewer.data().set_face_based(true);
    viewer.launch();
}
