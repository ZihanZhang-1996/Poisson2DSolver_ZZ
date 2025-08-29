#include "pde2d/fem_2D.hpp"
#include <Eigen/Sparse>
#include <cmath>

using namespace ::pde2d::mesh_2d;

namespace pde2d {
namespace fem_2d {

//
// 工具函数：把单元刚度矩阵 Ke 和单元残量 Re 累加到全局
// 参数：
//   ids   —— 单元对应的全局自由度编号
//   Ke    —— 单元刚度矩阵
//   Re    —— 单元残量向量
//   trips —— 全局稀疏矩阵的三元组 (I, J, val) 缓存
//   R     —— 全局残量向量
//   is_dirichlet —— 标记哪些节点是 Dirichlet 边界点
//
template<int N>
inline void add_element_to_global_sparse(
    const std::array<int,N>& ids,
    const Eigen::Matrix<double,N,N>& Ke,
    const Eigen::Matrix<double,N,1>& Re,
    std::vector<Eigen::Triplet<double>>& trips,
    Eigen::VectorXd& R,
    const std::vector<char>& is_dirichlet)
{
    for (int a=0; a<N; ++a) {
        const int I = ids[a];
        // 如果节点不是 Dirichlet 点，则累加残量
        if (!is_dirichlet[I]) R(I) += Re(a);

        for (int b=0; b<N; ++b) {
            const int J = ids[b];
            // 如果行/列中有 Dirichlet 自由度，则不加入三元组（后面单独处理）
            if (!is_dirichlet[I] && !is_dirichlet[J]) {
                trips.emplace_back(I, J, Ke(a,b));
            }
        }
    }
}

//
// 主装配函数：根据当前解 u，装配全局稀疏矩阵 K 和残量 R
// 参数：
//   mesh —— 网格
//   u    —— 当前解向量
//   K    —— 输出：稀疏全局刚度矩阵
//   R    —— 输出：全局残量
//   impose_dirichlet —— 是否施加 Dirichlet 边界条件
//   prob —— PDE 问题定义（包括 f, df/du, gD 边界条件）
//
void assembleForIteration(
    const ::pde2d::mesh_2d::Mesh2D& mesh,
    const Eigen::VectorXd& u,
    Eigen::SparseMatrix<double>& K,   // 输出: 稀疏矩阵
    Eigen::VectorXd& R,               // 输出: 残量
    bool impose_dirichlet,
    const ::pde2d::problem::ProblemExpr2D& prob)
{
    using SpMat = Eigen::SparseMatrix<double>;
    const int nnode = static_cast<int>(mesh.getNumNodes());    // 总节点数
    const int nelem = static_cast<int>(mesh.getNumElements()); // 总单元数
    assert(u.size() == nnode);

    R.setZero(nnode); // 初始化全局残量

    // 预分配稀疏矩阵结构
    K.resize(nnode, nnode);
    K.setZero();
    std::vector<Eigen::Triplet<double>> trips;
    // 每个节点大概有 5~7 个非零（二维 Laplace 型问题），预留容量
    trips.reserve(std::max(1, nnode * 7));

    // Dirichlet 掩码与边界值（如果 impose_dirichlet = true 则启用）
    std::vector<char> is_dirichlet(nnode, 0);
    std::vector<double> ubc(nnode, 0.0);

    if (impose_dirichlet) {
        auto nodes = mesh.getNodes();
        for (int i=0; i<nnode; ++i) {
            const double x = nodes[i]->getX();
            const double y = nodes[i]->getY();
            // 如果节点在边界上，标记并存储边界值 gD(x,y)
            if (prob.is_on_boundary(x,y)) {
                is_dirichlet[i] = 1;
                ubc[i] = prob.gD(x,y);
            }
        }
    }

    // 装配局部矩阵/向量到全局 
    auto elems = mesh.getElements();

    if (mesh.getElementType() == ::pde2d::mesh_2d::ElementType::RECTANGULAR) {
        // Q4 四节点矩形单元
        Eigen::Matrix4d Ke;
        Eigen::Vector4d Re;
        for (int e=0; e<nelem; ++e) {
            auto elem = elems[e];
            // 调用数值积分，计算单元刚度矩阵和残量
            quad_boost::computeElementStiffnessAndForcingOnRecElem_boost(elem, u, Ke, Re, prob);

            // 收集该单元的全局节点编号
            auto nodesInElem = elem->getNodes();
            std::array<int,4> ids;
            for (int i=0; i<4; ++i) ids[i] = nodesInElem[i]->getId();

            // 加入全局系统（跳过 Dirichlet 节点）
            add_element_to_global_sparse<4>(ids, Ke, Re, trips, R, is_dirichlet);
        }
    } else {
        // T3 三节点三角形单元
        Eigen::Matrix3d Ke;
        Eigen::Vector3d Re;
        for (int e=0; e<nelem; ++e) {
            auto elem = elems[e];
            quad_boost::computeElementStiffnessAndForcingOnTriElem_boost(elem, u, Ke, Re, prob);

            auto nodesInElem = elem->getNodes();
            std::array<int,3> ids;
            for (int i=0; i<3; ++i) ids[i] = nodesInElem[i]->getId();

            add_element_to_global_sparse<3>(ids, Ke, Re, trips, R, is_dirichlet);
        }
    }

    // Dirichlet 
    // 对于 Dirichlet 自由度：矩阵只保留 (i,i)=1，残量为 (u(i)-边界值)
    if (impose_dirichlet) {
        for (int i=0; i<nnode; ++i) {
            if (is_dirichlet[i]) {
                trips.emplace_back(i, i, 1.0);
                R(i) = u(i) - ubc[i];   // 保持和原稠密版本一致的处理方式
            }
        }
    }

    // 一次性生成稀疏矩阵
    K.setFromTriplets(trips.begin(), trips.end());
    K.makeCompressed();
}

} // namespace fem_2d
} // namespace pde2d
