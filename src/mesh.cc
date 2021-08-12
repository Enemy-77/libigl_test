#include "mesh.h"
#include <Eigen/sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>

using SpMatf = Eigen::SparseMatrix<double>;
using Tplf = Eigen::Triplet<double>;

Mesh::Mesh(const std::string& file_name) {
    if (!read_mesh(mesh_, file_name, ropt_)) {
        std::cerr << "Error loading file " << std::endl;
    }
    init_vertices();
}

void Mesh::update_vertices() {
    init_vertices();
}

void Mesh::init_vertices() {
    int num_vertices = mesh_.n_vertices();
    int num_faces = mesh_.n_faces();
    V_ = Eigen::MatrixXd(num_vertices, 3);
    F_ = Eigen::MatrixXi(num_faces, 3);

    int row = 0;
    for (MyMesh::VertexIter v_it = mesh_.vertices_begin(); v_it != mesh_.vertices_end(); ++v_it) {
        MyMesh::Point p = mesh_.point(v_it);
        V_(row, 0) = p[0];
        V_(row, 1) = p[1];
        V_(row, 2) = p[2];
        row++;
    }
    row = 0;
    for (MyMesh::FaceIter f_it = mesh_.faces_begin(); f_it != mesh_.faces_end(); ++f_it) {
        int col = 0;
        for (MyMesh::FaceVertexIter fvi = mesh_.fv_begin(*f_it); fvi.is_valid(); ++fvi) {
            MyMesh::Point p = mesh_.point(*fvi);
            F_(row, col) = fvi->idx();
            col++;
        }
        row++;
    }
}

static double vf_area(const MyMesh& mesh, const MyMesh::VertexIter& vertex) {
    double area = 0;
    typename MyMesh::VertexFaceIter vf_it;
    typename MyMesh::FaceVertexIter fv_it;
    for (vf_it = mesh.cvf_iter(vertex); vf_it; ++vf_it) {
        fv_it = mesh.cfv_iter(*vf_it);
        MyMesh::Point p = mesh.point(*fv_it); fv_it++;
        MyMesh::Point q = mesh.point(*fv_it); fv_it++;
        MyMesh::Point r = mesh.point(*fv_it);
        area += ((q - p) % (r - p)).norm() * 0.5;
    }
    return area;
}

static double angle(const MyMesh& mesh, const MyMesh::VertexOHalfedgeIter& voh_it, const MyMesh::Point& point) {
    MyMesh::Point to_point = mesh.point(mesh.to_vertex_handle(voh_it));
    HalfedgeHandle next = mesh.next_halfedge_handle(voh_it);
    MyMesh::Point next_point = mesh.point(mesh.to_vertex_handle(next));
    double d_1 = (point - to_point).sqrnorm();
    double d_2 = (point - next_point).sqrnorm();
    double d_3 = (to_point - next_point).sqrnorm();
    double cos_result = (d_2 + d_3 - d_1) / (2 * std::sqrt(d_2) * std::sqrt(d_3));
    return std::acos(cos_result);
}

static double angle(const MyMesh& mesh, const MyMesh::HalfedgeHandle& voh_it, const MyMesh::Point& point) {
    MyMesh::Point to_point = mesh.point(mesh.from_vertex_handle(voh_it));
    HalfedgeHandle next = mesh.next_halfedge_handle(voh_it);
    MyMesh::Point next_point = mesh.point(mesh.to_vertex_handle(next));
    double d_1 = (point - to_point).sqrnorm();
    double d_2 = (point - next_point).sqrnorm();
    double d_3 = (to_point - next_point).sqrnorm();
    double cos_result = (d_2 + d_3 - d_1) / (2 * std::sqrt(d_2) * std::sqrt(d_3));
    return std::acos(cos_result);
}

void Mesh::local_laplacian(float lambda) {
    MyMesh::VertexIter v_it, v_begin, v_end;
    v_begin = mesh_.vertices_begin();
    v_end = mesh_.vertices_end();
    for (v_it = v_begin; v_it != v_end; ++v_it) {
        if (v_it->is_boundary()) {
            continue;
        }
        double area = vf_area(mesh_, v_it);
        Vec3f sum(0, 0, 0);
        MyMesh::Point this_point = mesh_.point(v_it); // 当前的点
        MyMesh::VertexOHalfedgeIter voh_it; // outgoing 的 half_edge，该顶点指向的半边
        for (voh_it = mesh_.voh_iter(v_it); voh_it.is_valid(); ++voh_it) {
            MyMesh::Point to_point = mesh_.point(mesh_.to_vertex_handle(voh_it));  // Q点
            Vec3f q_p = to_point - this_point;
            MyMesh::HalfedgeHandle opposite = mesh_.opposite_halfedge_handle(voh_it);
            double alpha = angle(mesh_, voh_it, this_point);
            double beta = angle(mesh_, opposite, this_point);
            sum += (1 / std::tan(alpha) + 1 / std::tan(beta)) * q_p;
        }
        sum /= (4 * area);
        MyMesh::Point change = this_point + lambda * sum;
        mesh_.set_point(v_it, change);
    }
}

void Mesh::global_laplacian() {
    using pair_vec = std::vector<std::pair<int, double>>;
    int n = V_.rows();
    SpMatf sm(n, n);
    Eigen::VectorXd b_x(n);
    Eigen::VectorXd b_y(n);
    Eigen::VectorXd b_z(n);
   
    MyMesh::VertexIter v_it, v_begin, v_end;
    v_begin = mesh_.vertices_begin();
    v_end = mesh_.vertices_end();
    std::vector<Tplf> triplets;
    for (v_it = v_begin; v_it != v_end; ++v_it) {
        triplets.emplace_back(v_it->idx(), v_it->idx(), 1.f);
        if (v_it->is_boundary()) {
            Vec3f p = mesh_.point(v_it);
            b_x(v_it->idx()) = p[0];
            b_y(v_it->idx()) = p[1];
            b_z(v_it->idx()) = p[2];
        } else {
            MyMesh::Point this_point = mesh_.point(v_it);
            double sum = 0;
            pair_vec idx_sum;
            MyMesh::VertexOHalfedgeIter voh_it;
            for (voh_it = mesh_.voh_iter(v_it); voh_it.is_valid(); ++voh_it) {
                int idx = voh_it->to().idx();
                MyMesh::HalfedgeHandle opposite = voh_it->opp();
                double alpha = angle(mesh_, voh_it, this_point);
                double beta = angle(mesh_, opposite, this_point);
                double sum_up = 1.0f / std::tan(alpha) + 1.0f / std::tan(beta);
                idx_sum.emplace_back(idx, sum_up);
                sum += sum_up;
            }
            for (auto& pv : idx_sum) {
                triplets.emplace_back(v_it->idx(), pv.first, -pv.second / sum);
            }
            b_x(v_it->idx()) = 0.f;
            b_y(v_it->idx()) = 0.f;
            b_z(v_it->idx()) = 0.f;
        }
    }
    sm.setFromTriplets(triplets.begin(), triplets.end());
    Eigen::VectorXd result_x;
    Eigen::VectorXd result_y;
    Eigen::VectorXd result_z;
    //Eigen::SparseLU<SpMatf> solver;
    Eigen::SparseQR<SpMatf, Eigen::COLAMDOrdering<int>> solver;
    solver.compute(sm);
    result_x = solver.solve(b_x);
    result_y = solver.solve(b_y);
    result_z = solver.solve(b_z);
    for (MyMesh::VertexIter v_it = mesh_.vertices_begin(); v_it != mesh_.vertices_end(); ++v_it) {
        if (!v_it->is_boundary()) {
            int idx = (*v_it).idx();
            MyMesh::Point tp(result_x(idx), result_y(idx), result_z(idx));
            mesh_.set_point(*v_it, tp);
        }
    }
}

