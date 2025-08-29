#ifndef PDE2D_SOLVER_NEWTON_2D_HPP_INCLUDED
#define PDE2D_SOLVER_NEWTON_2D_HPP_INCLUDED

#include <Eigen/Dense>
#include <iostream>
#include <cmath>

#include "pde2d/Mesh2D.hpp"
#include "pde2d/fem_2D.hpp"

namespace pde2d { namespace problem { class ProblemExpr2D; } }

namespace pde2d {
namespace solver {

// 配置参数
struct Newton2DOptions {
  int    max_iter       = 30;     // 最大迭代步
  double tol_R          = 1e-10;  // 残量范数阈值
  double tol_rel_du         = 1e-12;  // 增量范数阈值
  double damping        = 1.0;    // 初始步长（线搜索起点）
  int    max_backtrack  = 0;     // 最大回溯次数
  double c1             = 1e-4;   // Armijo 系数
  bool   verbose        = true;   // 打印迭代日志
};

// 返回是否收敛；u 将被原地更新为解
bool newton_solve_2D(const pde2d::mesh_2d::Mesh2D& mesh,
                     Eigen::VectorXd& u,                  // in/out, size = nnode
                     const Newton2DOptions& opt,
                     const pde2d::problem::ProblemExpr2D& prob);

} // namespace solver
} // namespace pde2d

#endif // PDE2D_SOLVER_NEWTON_2D_HPP_INCLUDED
