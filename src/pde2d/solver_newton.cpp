#include "pde2d/solver_newton.hpp"
#include "pde2d/problem_syb.hpp"
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

namespace pde2d {
namespace solver {

// 使用 Newton 迭代法求解 2D 非线性方程组
bool newton_solve_2D(const pde2d::mesh_2d::Mesh2D& mesh,
                     Eigen::VectorXd& u,                      // 输入: 初始猜测; 输出: 最终解
                     const Newton2DOptions& opt,              // Newton 迭代参数 (容差、迭代数、阻尼系数等)
                     const pde2d::problem::ProblemExpr2D& prob) // PDE 问题定义 (f, df/du, 边界条件)
{
  using SpMat = Eigen::SparseMatrix<double>;

  const int n = static_cast<int>(u.size());
  SpMat K;                   // 全局稀疏雅可比/刚度矩阵
  Eigen::VectorXd R;         // 全局残量向量
  Eigen::VectorXd du(n);     // 本次牛顿步 (校正量)

  // === 第 0 步: 初值装配 ===
  // 根据初始猜测 u，装配全局矩阵 K 和残量 R (施加 Dirichlet 边界条件)
  pde2d::fem_2d::assembleForIteration(mesh, u, K, R, /*impose_dirichlet=*/true, prob);
  double Rnorm = R.norm();   // 初始残量范数
  if (opt.verbose) std::cout << "[Newton2D] iter=0, ||R||=" << Rnorm << "\n";

  // 稀疏矩阵求解器
  Eigen::SparseLU<SpMat> solver;
  // 符号分析 (非零结构分析) 只需做一次，后续迭代可以复用
  solver.analyzePattern(K);

  // === Newton 迭代主循环 ===
  for (int it = 1; it <= opt.max_iter; ++it) {
      // === 1. 解牛顿步 ===
      solver.factorize(K);
      if (solver.info() != Eigen::Success) {
        std::cerr << "[Newton2D] factorize failed at iter " << it << "\n";
        return false;
      }
      du = solver.solve(-R);
      if (solver.info() != Eigen::Success) {
        std::cerr << "[Newton2D] solve failed at iter " << it << "\n";
        return false;
      }

      // 计算牛顿步的相对范数 (pre-line-search)
      // const double rel_du = du.norm() / std::max(u.norm(), 1e-30);

      // === 2. 回溯线搜索 ===
      double alpha = opt.damping;
      const Eigen::VectorXd u_prev = u;
      SpMat K_trial;
      Eigen::VectorXd R_trial;
      bool accepted = false;

      for (int bt = 0; bt <= opt.max_backtrack; ++bt) {
        const Eigen::VectorXd u_trial_vec = u_prev + alpha * du;
        pde2d::fem_2d::assembleForIteration(mesh, u_trial_vec, K_trial, R_trial,
                                            /*impose_dirichlet=*/true, prob);
        const double Rnorm_trial = R_trial.norm();

        if (Rnorm_trial <= (1.0 - opt.c1 * alpha) * Rnorm) {
          u = u_trial_vec;
          K.swap(K_trial);
          R.swap(R_trial);
          Rnorm = Rnorm_trial;
          accepted = true;
          break;
        }
        alpha *= 0.5;
      }

      if (!accepted) {
        u = u_prev + du;
        pde2d::fem_2d::assembleForIteration(mesh, u, K, R, /*impose_dirichlet=*/true, prob);
        Rnorm = R.norm();
      }

      // === 3. 实际步长 (post-line-search) ===
      const double step_norm = (u - u_prev).norm();
      const double unorm     = std::max(u.norm(), 1e-30);
      const double rel_step  = step_norm / unorm;

      if (opt.verbose) {
        std::cout << "[Newton2D] iter=" << it
                  // << "  pre_rel_du=" << rel_du
                  << "  Rel_du=" << rel_step
                  << "  ||R||=" << Rnorm << "\n";
      }

      // === 4. 收敛判据 ===
      if ((Rnorm < opt.tol_R) && (rel_step < opt.tol_rel_du)) {
        if (opt.verbose) std::cout << "[Newton2D] Converged (both criteria).\n";
        return true;
      }
  }

  // 若迭代结束仍未满足收敛条件
  if (opt.verbose) std::cout << "[Newton2D] NOT converged in max_iter.\n";
  return false;
}

} // namespace solver
} // namespace pde2d
