#ifndef PDE2D_PROBLEM_EXPRTK_CLASS_INCLUDED
#define PDE2D_PROBLEM_EXPRTK_CLASS_INCLUDED

#include <string>
#include "exprtk.hpp"

namespace pde2d {
namespace problem {


class ProblemExpr2D {
public:
  explicit ProblemExpr2D(std::string f_expr = "sin(pi*x)*sin(pi*y)",
                         std::string dfdu_expr = "0");

  // 设置PDE的表达式
  bool set_f_expression(const std::string& expr);
  bool set_dfdu_expression(const std::string& expr);
  double a(double /*u*/, double /*x*/, double /*y*/) const { return 1.0; }
  double da_du(double /*u*/, double /*x*/, double /*y*/) const { return 0.0; }
  double f(double u, double x, double y) const;
  double df_du(double u, double x, double y) const;

  // 设置边条件
  bool is_on_boundary(double x, double y, double eps = 1e-12) const;

  //解析边条件表达式
  double gD(double x, double y) const;
  void set_domain_bounds(double xmin, double xmax, double ymin, double ymax) {
    xmin_ = xmin; xmax_ = xmax; ymin_ = ymin; ymax_ = ymax;
  }
  bool set_dirichlet_left (const std::string& expr);  
  bool set_dirichlet_right(const std::string& expr);  
  bool set_dirichlet_bottom(const std::string& expr); 
  bool set_dirichlet_top(const std::string& expr);


private:
  // 变量存储（被 ExprTk 以引用方式读取）
  mutable double vu_ = 0.0, vx_ = 0.0, vy_ = 0.0;

  // f 与 df/du 的符号表/表达式/解析器
  exprtk::symbol_table<double> st_f_, st_df_;
  exprtk::expression<double>   ex_f_,  ex_df_;
  exprtk::parser<double>       parser_f_, parser_df_;

  // 解析前的表达式
  std::string f_expr_;
  std::string dfdu_expr_;
  // 边界条件符号表/表达式/解析器
  exprtk::symbol_table<double> st_gL_, st_gR_, st_gB_, st_gT_;
  exprtk::expression<double>   ex_gL_,  ex_gR_,  ex_gB_,  ex_gT_;
  exprtk::parser<double>       parser_gL_, parser_gR_, parser_gB_, parser_gT_;
  std::string gL_expr_ = "0";
  std::string gR_expr_ = "0";
  std::string gB_expr_ = "0";
  std::string gT_expr_ = "0";

  // 域边界
  double xmin_ = 0.0, xmax_ = 1.0, ymin_ = 0.0, ymax_ = 1.0;
};

} // namespace problem
} // namespace pde2d

#endif // PDE2D_PROBLEM_EXPRTK_CLASS_INCLUDED
