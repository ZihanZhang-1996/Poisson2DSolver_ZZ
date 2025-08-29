#ifndef PDE2D_FEM_2D_HPP_INCLUDED
#define PDE2D_FEM_2D_HPP_INCLUDED

#include "pde2d/quadrature_boost.hpp"
#include <Eigen/Sparse>

namespace pde2d {
namespace fem_2d {

///  装配有限元全局矩阵和残量向量（Newton 迭代用）
///
/// - 输入：
///   - mesh : 网格对象（包含节点和单元）
///   - u    : 当前解向量（在所有节点上的值）
///   - prob : PDE 问题定义（包含 f(u,x,y), df/du 以及边界条件）
///   - impose_dirichlet : 是否施加 Dirichlet 边界条件
///
/// - 输出：
///   - K : 全局稀疏刚度矩阵（Jacobian 矩阵）
///   - R : 全局残量向量
///
/// - 说明：
///   1. 该函数会遍历所有单元，积分（quadrature_boost）
///      得到单元刚度矩阵 Ke 和单元残量/载荷向量 Fe；
///   2. 再将 Ke、Fe 按照局部节点到全局节点的映射，叠加到全局 K、R；
///   3. 若 impose_dirichlet=true，则对边界节点强制施加 Dirichlet 条件，
///      即将对应行改为单位行，并修改 R(i)=u(i)-gD。
///
void assembleForIteration(
    const ::pde2d::mesh_2d::Mesh2D& mesh,   // 网格（节点 + 单元）
    const Eigen::VectorXd& u,               // 当前解向量
    Eigen::SparseMatrix<double>& K,         // [输出] 全局刚度矩阵 (Jacobian)
    Eigen::VectorXd& R,                     // [输出] 全局残量向量
    bool impose_dirichlet,                  // 是否施加 Dirichlet 边界条件
    const ::pde2d::problem::ProblemExpr2D& prob); // PDE 问题定义

} // namespace fem_2d
} // namespace pde2d

#endif
