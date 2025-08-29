#include "pde2d/problem_syb.hpp"
#include <iostream>
#include <cmath>

namespace pde2d {
namespace problem {

// 构造函数：初始化 PDE 表达式 f(u,x,y)、df/du(u,x,y)，并编译成 ExprTk 表达式
ProblemExpr2D::ProblemExpr2D(std::string f_expr, std::string dfdu_expr)
  : f_expr_(std::move(f_expr)), dfdu_expr_(std::move(dfdu_expr))
{
    // ====== 设置 f(u,x,y) ======
    // ExprTk 使用符号表(symbol_table)，需要把变量作为引用绑定进去
    st_f_.add_variable("u", vu_);
    st_f_.add_variable("x", vx_);
    st_f_.add_variable("y", vy_);
    st_f_.add_constants();               // 加入内置常数 (如 pi, e)
    ex_f_.register_symbol_table(st_f_);

    // ====== 设置 df/du(u,x,y) ======
    st_df_.add_variable("u", vu_);
    st_df_.add_variable("x", vx_);
    st_df_.add_variable("y", vy_);
    st_df_.add_constants();
    ex_df_.register_symbol_table(st_df_);

    // 编译 f 和 df/du 表达式
    if (!parser_f_.compile(f_expr_, ex_f_)) {
      std::cerr << "[ExprTk] compile failed for f: " << f_expr_ << "\n";
    }
    if (!parser_df_.compile(dfdu_expr_, ex_df_)) {
      std::cerr << "[ExprTk] compile failed for df/du: " << dfdu_expr_ << "\n";
    }

    // ====== 绑定 Dirichlet 边界条件函数 ======
    auto bind = [&](exprtk::symbol_table<double>& st, exprtk::expression<double>& ex){
      st.add_variable("x", vx_);
      st.add_variable("y", vy_);
      st.add_constants();
      ex.register_symbol_table(st);
    };
    // 初始时边界条件都是默认表达式 (0)，但要先绑定 symbol_table
    bind(st_gL_, ex_gL_); parser_gL_.compile(gL_expr_, ex_gL_);
    bind(st_gR_, ex_gR_); parser_gR_.compile(gR_expr_, ex_gR_);
    bind(st_gB_, ex_gB_); parser_gB_.compile(gB_expr_, ex_gB_);
    bind(st_gT_, ex_gT_); parser_gT_.compile(gT_expr_, ex_gT_);
}

// 设置 f(u,x,y) 的表达式
bool ProblemExpr2D::set_f_expression(const std::string& expr) {
    f_expr_ = expr;
    if (!parser_f_.compile(f_expr_, ex_f_)) {
      std::cerr << "[ExprTk] compile failed for f: " << f_expr_ << "\n";
      return false;
    }
    return true;
}

// 设置 df/du(u,x,y) 的表达式
bool ProblemExpr2D::set_dfdu_expression(const std::string& expr) {
    dfdu_expr_ = expr;
    if (!parser_df_.compile(dfdu_expr_, ex_df_)) {
      std::cerr << "[ExprTk] compile failed for df/du: " << dfdu_expr_ << "\n";
      return false;
    }
    return true;
}

// 计算 f(u,x,y)
double ProblemExpr2D::f(double u, double x, double y) const {
    vu_ = u; vx_ = x; vy_ = y;
    return ex_f_.value();
}

// 计算 df/du(u,x,y)
double ProblemExpr2D::df_du(double u, double x, double y) const {
    vu_ = u; vx_ = x; vy_ = y;
    return ex_df_.value();
}

// ====== Dirichlet 边界条件设置 ======
// 允许用户分别为左/右/下/上边界设置表达式
// 例如 "sin(pi*y)"，"x*x+y*y" 等
bool ProblemExpr2D::set_dirichlet_left (const std::string& e){ gL_expr_ = e; return parser_gL_.compile(gL_expr_, ex_gL_); }
bool ProblemExpr2D::set_dirichlet_right(const std::string& e){ gR_expr_ = e; return parser_gR_.compile(gR_expr_, ex_gR_); }
bool ProblemExpr2D::set_dirichlet_bottom(const std::string& e){ gB_expr_ = e; return parser_gB_.compile(gB_expr_, ex_gB_); }
bool ProblemExpr2D::set_dirichlet_top   (const std::string& e){ gT_expr_ = e; return parser_gT_.compile(gT_expr_, ex_gT_); }

//
// 判断一个点 (x,y) 是否在边界上
bool ProblemExpr2D::is_on_boundary(double x, double y, double eps) const {
    const double Lx  = xmax_ - xmin_;
    const double Ly  = ymax_ - ymin_;
    const double tol = std::max({1.0, std::abs(Lx), std::abs(Ly)}) * eps;

    const bool onL = std::abs(x - xmin_) <= tol;
    const bool onR = std::abs(x - xmax_) <= tol;
    const bool onB = std::abs(y - ymin_) <= tol;
    const bool onT = std::abs(y - ymax_) <= tol;

    return onL || onR || onB || onT;
}


// 计算 Dirichlet 边界条件值 gD(x,y)
// 根据点所在边界 (L/R/B/T)，返回对应的边界表达式值
double ProblemExpr2D::gD(double x, double y) const {
  vx_ = x; vy_ = y;
  // 容差 ex：根据网格尺度自动缩放
  const double ex = std::max({1.0, xmax_-xmin_, ymax_-ymin_}) * 1e-12;
  const bool onL = (std::abs(x - xmin_) < ex);
  const bool onR = (std::abs(x - xmax_) < ex);
  const bool onB = (std::abs(y - ymin_) < ex);
  const bool onT = (std::abs(y - ymax_) < ex);

  if (onL) return ex_gL_.value();
  if (onR) return ex_gR_.value();
  if (onB) return ex_gB_.value();
  if (onT) return ex_gT_.value();
  // 理论上不会走到这里（只在边界节点调用）
  return 0.0;
}

} // namespace problem
} // namespace pde2d
