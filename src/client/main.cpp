#include <iostream>
#include <Eigen/Dense>
#include "pde2d/Mesh2D.hpp"
#include "pde2d/problem_syb.hpp"
#include "pde2d/fem_2D.hpp"
#include "pde2d/solver_newton.hpp"
#include "pde2d/io.hpp"
#include <filesystem>
#include <chrono>                  // 用于计时

#if defined(__APPLE__)
  #include <mach-o/dyld.h>
#elif defined(_WIN32)
  #include <windows.h>
#else
  #include <unistd.h>
  #include <limits.h>
#endif

// 获取当前可执行文件所在目录
static std::filesystem::path exe_dir() {
#if defined(__APPLE__)
    uint32_t size = 0;
    _NSGetExecutablePath(nullptr, &size);
    std::string buf(size, '\0');
    _NSGetExecutablePath(buf.data(), &size);
    return std::filesystem::weakly_canonical(std::filesystem::path(buf).parent_path());
#elif defined(_WIN32)
    wchar_t buf[MAX_PATH];
    GetModuleFileNameW(nullptr, buf, MAX_PATH);
    return std::filesystem::path(buf).parent_path();
#else  // Linux/Unix
    char buf[PATH_MAX];
    ssize_t len = readlink("/proc/self/exe", buf, sizeof(buf)-1);
    if (len > 0) { buf[len] = '\0'; }
    return std::filesystem::weakly_canonical(std::filesystem::path(buf).parent_path());
#endif
}

int main(int argc, char* argv[]) {
    using namespace pde2d::mesh_2d;
    using clock_t = std::chrono::steady_clock;

    // --- 命令行输入检查 ---
    if (argc < 2) {
        std::cerr << "用法: " << argv[0] << " <输入文件名.json>\n"
                  << "注意: 输入文件应放在 ./input/ 文件夹中\n";
        return 1;
    }

    // --- 构造 input 文件完整路径 ---
    std::string filename = argv[1];  // 用户传入的文件名（例如 example_input.json）
    std::filesystem::path input_file = std::filesystem::path("input") / filename;

    // 程序开始计时
    const auto t_program_start = clock_t::now();

    // --- 定义求解所需对象 ---
    pde2d::mesh_2d::Mesh2D mesh;                // 网格
    pde2d::problem::ProblemExpr2D prob;         // PDE 问题定义（f, df/du, 边界条件等）
    pde2d::solver::Newton2DOptions opt;         // Newton 迭代选项（容差、最大迭代数等）
    Eigen::VectorXd u;                          // 解向量（节点自由度）
    opt.verbose = true; // 打开迭代日志

    // 从 JSON 文件读取配置（网格、问题定义、初始条件、Newton 参数等）
    if (!pde2d::io::load_config_from_json(input_file.string(), mesh, prob, u, opt)) {
        std::cerr << "配置文件加载失败: " << input_file << "\n";
        return 1;
    }

    // 调用 Newton 迭代求解器
    const auto t_newton_start = clock_t::now();
    bool ok = pde2d::solver::newton_solve_2D(mesh, u, opt, prob);
    const auto t_newton_end   = clock_t::now();
    using ms = std::chrono::milliseconds;
    const auto newton_ms = std::chrono::duration_cast<ms>(t_newton_end - t_newton_start).count();

    // 打印收敛信息
    std::cout << (ok ? "Converged\n" : "Not converged\n");
    std::cout << "[Timing] Newton solve time: " << newton_ms << " ms\n";

    // --- 导出数值解到 CSV ---
    // 输出文件名格式: input文件名的stem + "_u_matrix.csv"
    std::string stem = std::filesystem::path(filename).stem().string();
    auto out_csv = exe_dir() / "output" / (stem + "_u_matrix.csv");

    // 保存矩阵形式的数值解
    std::filesystem::create_directories(out_csv.parent_path());
    pde2d::io::save_u_matrix_csv(mesh, u, out_csv.string());

    // --- 全程序耗时统计 ---
    std::cout << "[Mesh] Mesh: " << mesh.getnx() <<"x"<< mesh.getny()<< "\n";
    const auto t_program_end = clock_t::now();
    const auto program_ms = std::chrono::duration_cast<ms>(t_program_end - t_program_start).count();
    std::cout << "[Timing] Total program time: " << program_ms << " ms\n";

    return 0;
}