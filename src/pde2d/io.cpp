#include "pde2d/io.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include "json.hpp"

namespace pde2d {
namespace io {

//
// 从 JSON 文件中读取 problem 部分 (仅 f, df/du 表达式)
//
bool load_problem_from_json(const std::string& filename,
                            pde2d::problem::ProblemExpr2D& prob)
{
    try {
        std::ifstream ifs(filename);
        if (!ifs) {
            std::cerr << "[IO] Cannot open " << filename << "\n";
            return false;
        }

        nlohmann::json j;
        ifs >> j;

        // 必须包含 "problem" 节
        if (!j.contains("problem")) {
            std::cerr << "[IO] JSON missing 'problem' section\n";
            return false;
        }

        const auto& jp = j["problem"];
        if (!jp.contains("f_expr") || !jp.contains("dfdu_expr")) {
            std::cerr << "[IO] JSON missing f_expr/dfdu_expr\n";
            return false;
        }

        std::string f_expr    = jp["f_expr"];
        std::string dfdu_expr = jp["dfdu_expr"];

        if (!prob.set_f_expression(f_expr)) return false;
        if (!prob.set_dfdu_expression(dfdu_expr)) return false;

        std::cout << "[IO] Loaded problem from " << filename << "\n";
        return true;
    }
    catch (const std::exception& e) {
        std::cerr << "[IO] Failed to parse JSON: " << e.what() << "\n";
        return false;
    }
}

//
// 从 JSON 读取 Newton 迭代选项 (max_iter, tol_R, tol_rel_du)
//
static void read_newton_options(const nlohmann::json& j,
                                pde2d::solver::Newton2DOptions& opt)
{
    if (!j.contains("newton")) return;
    const auto& jn = j["newton"];
    if (jn.contains("max_iter"))    opt.max_iter   = jn["max_iter"].get<int>();
    if (jn.contains("tol_R"))       opt.tol_R      = jn["tol_R"].get<double>();
    if (jn.contains("tol_rel_du"))  opt.tol_rel_du = jn["tol_rel_du"].get<double>();
}

//
// 根据输入表达式计算初始条件 u0 (在每个节点上求值)
// 如果 expr 为空，默认初始解为全零
//
bool evaluate_initial_condition(const pde2d::mesh_2d::Mesh2D& mesh,
                                const std::string& expr,
                                Eigen::VectorXd& u0)
{
    const int nnode = static_cast<int>(mesh.getNumNodes());
    u0.resize(nnode);

    // 如果没给表达式，直接返回零向量
    if (expr.empty()) {
        u0.setZero();
        return true;
    }

    // 用 ExprTk 编译表达式
    double vx = 0.0, vy = 0.0;
    exprtk::symbol_table<double> sym;
    sym.add_variable("x", vx);
    sym.add_variable("y", vy);
    sym.add_constants();

    exprtk::expression<double> ex;
    ex.register_symbol_table(sym);

    exprtk::parser<double> parser;
    if (!parser.compile(expr, ex)) {
        std::cerr << "[IO] ExprTk compile failed for initial_condition: " << expr << "\n";
        return false;
    }

    // 在所有节点上计算初值
    for (int i = 0; i < nnode; ++i) {
        vx = mesh.getNodes()[i]->getX();
        vy = mesh.getNodes()[i]->getY();
        u0[i] = ex.value();
    }
    return true;
}

//
// 从 JSON 配置文件中读取：
// 1. PDE problem (f, df/du)
// 2. 网格 (nx, ny, Lx, Ly, type)
// 3. 初始条件 (optional)
// 4. Dirichlet 边界条件
// 5. Newton 迭代参数
//
bool load_config_from_json(const std::string& filename,
                           pde2d::mesh_2d::Mesh2D& mesh,
                           pde2d::problem::ProblemExpr2D& prob,
                           Eigen::VectorXd& u0,
                           pde2d::solver::Newton2DOptions& newton_opt)
{
    try {
        std::ifstream ifs(filename);
        if (!ifs) { std::cerr << "[IO] Cannot open " << filename << "\n"; return false; }
        nlohmann::json j; ifs >> j;

        // === 1. problem ===
        if (!j.contains("problem")) { std::cerr << "[IO] JSON missing 'problem'\n"; return false; }
        {
            const auto& jp = j["problem"];
            std::string f_expr    = jp["f_expr"];
            std::string dfdu_expr = jp["dfdu_expr"];
            if (!prob.set_f_expression(f_expr)) return false;
            if (!prob.set_dfdu_expression(dfdu_expr)) return false;
        }

        // === 2. mesh ===
        int nx = j.value("mesh", nlohmann::json{}).value("nx", 16);
        int ny = j.value("mesh", nlohmann::json{}).value("ny", 16);
        double Lx = j.value("mesh", nlohmann::json{}).value("Lx", 1.0);
        double Ly = j.value("mesh", nlohmann::json{}).value("Ly", 1.0);
        std::string type_str = j.value("mesh", nlohmann::json{}).value("type", "triangular");

        pde2d::mesh_2d::ElementType elem_type;
        if (type_str == "triangular") {
            elem_type = pde2d::mesh_2d::ElementType::TRIANGULAR;
        } else if (type_str == "rectangular") {
            elem_type = pde2d::mesh_2d::ElementType::RECTANGULAR;
        } else {
            std::cerr << "[IO] Unknown mesh type: " << type_str
                      << ", fallback to triangular.\n";
            elem_type = pde2d::mesh_2d::ElementType::TRIANGULAR;
        }

        // 创建网格对象
        mesh = pde2d::mesh_2d::Mesh2D(Lx, Ly, nx, ny, elem_type);

        // === 3. initial condition (可选) ===
        std::string ic_expr;
        if (j.contains("initial_condition")) {
            const auto& ji = j["initial_condition"];
            if (ji.is_string()) ic_expr = ji.get<std::string>();
            else if (ji.is_object() && ji.contains("expr")) ic_expr = ji["expr"].get<std::string>();
        }
        if (!evaluate_initial_condition(mesh, ic_expr, u0)) return false;

        std::cout << "[IO] Loaded config from " << filename << "\n";

        // 设置 PDE 的域范围
        prob.set_domain_bounds(0.0, Lx, 0.0, Ly);

        // === 4. Dirichlet 边界条件 ===
        if (j.contains("boundary")) {
            const auto& jb = j["boundary"];
            const std::string type = jb.value("type", "dirichlet");
            if (type != "dirichlet") {
                std::cerr << "[IO] only dirichlet boundary supported currently\n";
                return false;
            }
            if (jb.contains("left"))   prob.set_dirichlet_left  (jb["left"].get<std::string>());
            if (jb.contains("right"))  prob.set_dirichlet_right (jb["right"].get<std::string>());
            if (jb.contains("bottom")) prob.set_dirichlet_bottom(jb["bottom"].get<std::string>());
            if (jb.contains("top"))    prob.set_dirichlet_top   (jb["top"].get<std::string>());
        }

        // === 5. Newton 迭代选项 ===
        read_newton_options(j, newton_opt);
        if (newton_opt.verbose) {
            std::cout << "[IO] Newton opts: max_iter=" << newton_opt.max_iter
                      << ", tol_R=" << newton_opt.tol_R
                      << ", tol_rel_du=" << newton_opt.tol_rel_du << "\n";
        }

        return true;
    }
    catch (const std::exception& e) {
        std::cerr << "[IO] Failed to parse JSON: " << e.what() << "\n"; 
        return false;
    }
}

//
// 保存解向量 u 到 CSV 文件
// 输出为 (ny+1) 行, (nx+1) 列的矩阵形式
// flip_y = true 时，y 方向翻转（通常用于可视化）
//
void save_u_matrix_csv(const pde2d::mesh_2d::Mesh2D& mesh,
                       const Eigen::VectorXd& u,
                       const std::string& filename,
                       int precision,
                       bool flip_y)
{
    const int nx = mesh.getnx() + 1;
    const int ny = mesh.getny() + 1;
    if (u.size() != nx * ny) {
        std::cerr << "save_u_matrix_csv: size mismatch (u.size()="
                  << u.size() << ", expected " << nx*ny << ")\n";
        return;
    }

    std::ofstream ofs(filename);
    if (!ofs) {
        std::cerr << "save_u_matrix_csv: cannot open file '" << filename << "'\n";
        return;
    }

    ofs.setf(std::ios::fixed);
    ofs << std::setprecision(precision);

    if (!flip_y) {
        // 正常顺序
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                ofs << u[j*nx + i];
                if (i + 1 < nx) ofs << ",";
            }
            ofs << "\n";
        }
    } else {
        // 翻转 y 轴顺序
        for (int j = ny - 1; j >= 0; --j) {
            for (int i = 0; i < nx; ++i) {
                ofs << u[j*nx + i];
                if (i + 1 < nx) ofs << ",";
            }
            ofs << "\n";
        }
    }
}

} // namespace io
} // namespace pde2d
