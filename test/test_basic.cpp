#include <cassert>
#include <iostream>
#include <Eigen/Dense>
#include "pde2d/Mesh2D.hpp"
#include "pde2d/problem_syb.hpp"
#include "pde2d/fem_2D.hpp"
#include "pde2d/solver_newton.hpp"

using namespace pde2d;

//
// Unit Test 1: 网格节点数量
// 测试生成 2x2 网格时节点数和单元数是否正确
//
void test_mesh_basic() {
    mesh_2d::Mesh2D mesh(1.0, 1.0, 2, 2, mesh_2d::ElementType::TRIANGULAR);
    assert(mesh.getNumNodes() == 9);    // (nx+1)*(ny+1) = 9
    assert(mesh.getNumElements() == 8); // 每个小方格两个三角形，2*2*2=8
    std::cout << "[PASS] test_mesh_basic\n";
}

//
// Unit Test 2: 简单的常数边界条件
// 测试 Dirichlet 边界条件是否正确施加
//
void test_dirichlet_bc() {
    mesh_2d::Mesh2D mesh(1.0, 1.0, 1, 1, mesh_2d::ElementType::RECTANGULAR);
    problem::ProblemExpr2D prob("0", "0"); // f(u)=0, df/du=0
    prob.set_domain_bounds(0,1,0,1);
    prob.set_dirichlet_left("1");  // 左边界值=1
    prob.set_dirichlet_right("1");
    prob.set_dirichlet_bottom("1");
    prob.set_dirichlet_top("1");

    Eigen::VectorXd u = Eigen::VectorXd::Zero(mesh.getNumNodes());
    Eigen::SparseMatrix<double> K;
    Eigen::VectorXd R;

    fem_2d::assembleForIteration(mesh, u, K, R, true, prob);

    // 所有边界点的残差 R(i) 应该等于 (u(i) - 1) = -1
    auto nodes = mesh.getNodes();
    for (int i=0; i<mesh.getNumNodes(); i++) {
        double x = nodes[i]->getX();
        double y = nodes[i]->getY();
        if (prob.is_on_boundary(x,y)) {
            assert(std::abs(R(i) + 1.0) < 1e-12);
        }
    }
    std::cout << "[PASS] test_dirichlet_bc\n";
}

int main() {
    test_mesh_basic();
    test_dirichlet_bc();
    return 0;
}
