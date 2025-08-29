#ifndef PDE2D_IO_2D_HPP_INCLUDED
#define PDE2D_IO_2D_HPP_INCLUDED

#include <string>
#include <Eigen/Dense>
#include "pde2d/Mesh2D.hpp"
#include "pde2d/problem_syb.hpp"
#include "pde2d/solver_newton.hpp"

namespace pde2d {
namespace io {

///  将数值解 u 按照网格拓扑保存为 CSV 矩阵格式
/// 
/// - 输入:
///   - mesh     : 网格 (包含 Nx, Ny, 节点坐标)
///   - u        : 解向量 (长度应为 (Nx+1)*(Ny+1))
///   - filename : 输出文件名 (默认 "u_matrix.csv")
///   - precision: 输出小数位数 (默认 16)
///   - flip_y   : 是否翻转 y 方向 (默认 false)
/// 
/// - 输出:
///   在 CSV 文件中保存解矩阵 (按网格节点顺序排列)
///
void save_u_matrix_csv(const pde2d::mesh_2d::Mesh2D& mesh,
                       const Eigen::VectorXd& u,
                       const std::string& filename = "u_matrix.csv",
                       int precision = 16,
                       bool flip_y = false);


/// 从 JSON 文件读取问题定义 (f(u,x,y), df/du)
/// 
/// - 输入:
///   - filename : JSON 文件路径
/// - 输出:
///   - prob     : ProblemExpr2D 对象，存储源项 f 和导数 df/du 的表达式
/// 
/// - 返回值:
///   true  = 成功解析
///   false = 出错 (文件不存在或解析失败)
///
bool load_problem_from_json(const std::string& filename,
                            pde2d::problem::ProblemExpr2D& prob);


///  从 JSON 文件加载配置 (网格 + PDE + 初值 + 边界条件 + Newton 参数)
/// 
/// - 输入:
///   - filename   : JSON 文件路径
/// - 输出:
///   - mesh       : 构造完成的网格对象
///   - prob       : PDE 问题定义 (f, df/du, 边界条件 gD)
///   - u0         : 初始解向量 (由 initial_condition 表达式生成，或默认 0)
///   - newton_opt : Newton 迭代选项 (max_iter, tol_R, tol_rel_du 等)
/// 
/// - 返回值:
///   true  = 成功
///   false = 出错 (文件不存在或 JSON 格式不正确)
///
bool load_config_from_json(const std::string& filename,
                           pde2d::mesh_2d::Mesh2D& mesh,
                           pde2d::problem::ProblemExpr2D& prob,
                           Eigen::VectorXd& u0,
                           pde2d::solver::Newton2DOptions& newton_opt);


///  根据用户给定的解析表达式 (ExprTk) 生成初始解
/// 
/// - 输入:
///   - mesh : 网格 (提供节点坐标)
///   - expr : 初值表达式 (例如 "sin(pi*x)*sin(pi*y)")
/// - 输出:
///   - u0   : 在每个网格节点上评估后的数值
/// 
/// - 返回值:
///   true  = 成功 (表达式编译和计算正常)
///   false = 出错 (表达式编译失败)
///
bool evaluate_initial_condition(const pde2d::mesh_2d::Mesh2D& mesh,
                                const std::string& expr,
                                Eigen::VectorXd& u0);

} // namespace io
} // namespace pde2d

#endif // PDE2D_IO_2D_HPP_INCLUDED
