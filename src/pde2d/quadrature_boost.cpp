#include "pde2d/fem_2D.hpp"
#include "pde2d/quadrature_boost.hpp"
#include <cmath>

using namespace ::pde2d::mesh_2d;

namespace pde2d {
namespace fem_2d {
namespace quad_boost {

// ----------------- 形函数工具 -----------------
namespace ShapeFunctions {
    
    // 三角形线性形函数 φ1=1-ξ-η, φ2=ξ, φ3=η
    inline void triangle(double xi, double eta, double phi[3]) {
        phi[0] = 1.0 - xi - eta;
        phi[1] = xi;
        phi[2] = eta;
    }
    
    // 三角形线性形函数的导数（常数，不依赖 ξ,η）
    inline void triangleDerivatives(double dphi_dxi[3], double dphi_deta[3]) {
        dphi_dxi[0] = -1.0; dphi_dxi[1] = 1.0; dphi_dxi[2] = 0.0;
        dphi_deta[0] = -1.0; dphi_deta[1] = 0.0; dphi_deta[2] = 1.0;
    }
    
    // 四节点矩形（双线性）形函数
    // 顺序: 左下, 右下, 右上, 左上
    inline void rectangle(double xi, double eta, double phi[4]) {
        phi[0] = 0.25 * (1 - xi) * (1 - eta);
        phi[1] = 0.25 * (1 + xi) * (1 - eta);
        phi[2] = 0.25 * (1 + xi) * (1 + eta);
        phi[3] = 0.25 * (1 - xi) * (1 + eta);
    }
    
    // 矩形形函数的导数
    inline void rectangleDerivatives(double xi, double eta, 
                                   double dphi_dxi[4], double dphi_deta[4]) {
        dphi_dxi[0] = -0.25 * (1 - eta);
        dphi_dxi[1] =  0.25 * (1 - eta);
        dphi_dxi[2] =  0.25 * (1 + eta);
        dphi_dxi[3] = -0.25 * (1 + eta);
        
        dphi_deta[0] = -0.25 * (1 - xi);
        dphi_deta[1] = -0.25 * (1 + xi);
        dphi_deta[2] =  0.25 * (1 + xi);
        dphi_deta[3] =  0.25 * (1 - xi);
    }
} // namespace ShapeFunctions


// ----------------- 雅可比矩阵计算 -----------------
// J = [dx/dξ  dx/dη; dy/dξ  dy/dη], detJ = det(J)
// 对于线性三角形和矩形，Jacobian 在整个单元内是常数
void computeJacobian(const std::shared_ptr<::pde2d::mesh_2d::Element> elem, Eigen::Matrix2d& J, double& detJ) {
    double dphi_dxi[4], dphi_deta[4];
    
    if (elem->getType() == pde2d::mesh_2d::ElementType::TRIANGULAR) {
        ShapeFunctions::triangleDerivatives(dphi_dxi, dphi_deta);
    } else {
        ShapeFunctions::rectangleDerivatives(0, 0, dphi_dxi, dphi_deta); // 在中心点
    }
    
    J.setZero();
    for (int i = 0; i < elem->getNumNodes(); i++) {
        const std::shared_ptr<Node> node = elem->getNode(i);
        J(0,0) += dphi_dxi[i]  * node->getX();
        J(0,1) += dphi_dxi[i]  * node->getY();
        J(1,0) += dphi_deta[i] * node->getX();
        J(1,1) += dphi_deta[i] * node->getY();
    }
    
    detJ = J(0,0)*J(1,1) - J(0,1)*J(1,0);
}

// ----------------- 坐标变换 (参考域 → 物理域) -----------------
void transformToPhysical(const std::shared_ptr<::pde2d::mesh_2d::Element> elem, double xi, double eta, 
                         double& x, double& y, double* phi) {
    // 计算参考域上的形函数
    if (elem->getType() == pde2d::mesh_2d::ElementType::TRIANGULAR) {
        ShapeFunctions::triangle(xi, eta, phi);
    } else {
        ShapeFunctions::rectangle(xi, eta, phi);
    }
    
    // 根据形函数插值，得到物理坐标 (x,y)
    x = 0.0; y = 0.0;
    for (int i = 0; i < elem->getNumNodes(); i++) {
        const std::shared_ptr<Node> node = elem->getNode(i);
        x += phi[i] * node->getX();
        y += phi[i] * node->getY();
    }
}

// ----------------- 形函数梯度 (参考坐标系 → 物理坐标系) -----------------
void computeShapeFunctionGradients(const std::shared_ptr<::pde2d::mesh_2d::Element> elem, const Eigen::Matrix2d& invJ, 
                                   Eigen::MatrixXd& grad_phi) {
    double dphi_dxi[4], dphi_deta[4];
    
    if (elem->getType() == pde2d::mesh_2d::ElementType::TRIANGULAR) {
        ShapeFunctions::triangleDerivatives(dphi_dxi, dphi_deta);
    } else {
        ShapeFunctions::rectangleDerivatives(0, 0, dphi_dxi, dphi_deta);
    }
    
    grad_phi.resize(2, elem->getNumNodes());
    for (int i = 0; i < elem->getNumNodes(); i++) {
        grad_phi(0,i) = invJ(0,0)*dphi_dxi[i] + invJ(0,1)*dphi_deta[i];  // dφ/dx
        grad_phi(1,i) = invJ(1,0)*dphi_dxi[i] + invJ(1,1)*dphi_deta[i];  // dφ/dy
    }
}

// ----------------- 三角形单元刚度矩阵 & 残量装配 (高斯积分) -----------------
void computeElementStiffnessAndForcingOnTriElem_boost(
    const std::shared_ptr<::pde2d::mesh_2d::Element> elem,
    const Eigen::VectorXd& u,
    Eigen::Matrix3d& Ke, Eigen::Vector3d& Fe,
    const ::pde2d::problem::ProblemExpr2D& prob
) {
    // Jacobian (三角形常数)
    Eigen::Matrix2d J;
    double detJ;
    computeJacobian(elem, J, detJ);
    double abs_detJ = std::abs(detJ);
    Eigen::Matrix2d invJ = J.inverse();
    
    // 形函数梯度
    Eigen::MatrixXd grad_phi;
    computeShapeFunctionGradients(elem, invJ, grad_phi);
    
    // 当前单元节点解向量
    Eigen::Vector3d u_elem;
    for (int i = 0; i < 3; i++) {
        u_elem(i) = u[elem->getNode(i)->getId()];
    }
    
    Ke.setZero();
    Fe.setZero();
    
    // 使用 Boost.Math 高斯积分 (三角形积分)
    constexpr int order = 3; // 积分阶数
    using rule = ::boost::math::quadrature::gauss<double, order>;
    
    // 将 [-1,1] → [0,1] 映射，用于三角形自然坐标
    for (int i = 0; i < order; i++) {
        double xi_ref = rule::abscissa()[i];
        double w1 = rule::weights()[i];
        double xi = 0.5 * (xi_ref + 1.0); 
        double jacobian1 = 0.5; // dξ/dξ_ref
        
        for (int j = 0; j < order; j++) {
            double eta_ref = rule::abscissa()[j];
            double w2 = rule::weights()[j];
            double eta = 0.5 * (eta_ref + 1.0) * (1.0 - xi);
            double jacobian2 = 0.5 * (1.0 - xi);
            
            double weight = w1 * w2 * jacobian1 * jacobian2 * abs_detJ;
            
            if (xi + eta > 1.0) continue; // 超出三角形范围则跳过
            
            // 形函数 + 物理坐标
            double phi[3];
            double x, y;
            transformToPhysical(elem, xi, eta, x, y, phi);
            
            // 插值解和梯度
            double uq = 0.0;
            Eigen::Vector2d grad_uq = Eigen::Vector2d::Zero();
            for (int k = 0; k < 3; k++) {
                uq        += phi[k] * u_elem(k);
                grad_uq(0) += grad_phi(0,k) * u_elem(k);
                grad_uq(1) += grad_phi(1,k) * u_elem(k);
            }
            
            // 源项 f 和 df/du
            double fq     = prob.f(uq, x, y);
            double dfdu_q = prob.df_du(uq, x, y);
            
            // 载荷向量
            for (int k = 0; k < 3; k++) {
                Fe(k) += (grad_phi(0,k)*grad_uq(0) + grad_phi(1,k)*grad_uq(1) - fq * phi[k]) * weight;
            }
            
            // 刚度矩阵
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    double term = (grad_phi(0,k)*grad_phi(0,l) + grad_phi(1,k)*grad_phi(1,l)) 
                                  - phi[k] * phi[l] * dfdu_q;
                    Ke(k,l) += term * weight;
                }
            }
        }
    }
}

// ----------------- 矩形单元刚度矩阵 & 残量装配 (高斯积分) -----------------
void computeElementStiffnessAndForcingOnRecElem_boost(
    const std::shared_ptr<::pde2d::mesh_2d::Element> elem,
    const Eigen::VectorXd& u,
    Eigen::Matrix4d& Ke,Eigen::Vector4d& Fe,
    const ::pde2d::problem::ProblemExpr2D& prob
) {
    // Jacobian (矩形常数)
    Eigen::Matrix2d J;
    double detJ;
    computeJacobian(elem, J, detJ);
    double abs_detJ = std::abs(detJ);
    Eigen::Matrix2d invJ = J.inverse();
    
    // 形函数梯度
    Eigen::MatrixXd grad_phi;
    computeShapeFunctionGradients(elem, invJ, grad_phi);
    
    // 当前单元节点解
    Eigen::Vector4d u_elem;
    for (int i = 0; i < 4; i++) {
        u_elem(i) = u[elem->getNode(i)->getId()];
    }
    
    Ke.setZero();
    Fe.setZero();
    
    // Boost 高斯积分 (矩形在 [-1,1]x[-1,1])
    const int order = 3;
    using rule = ::boost::math::quadrature::gauss<double, order>;
    
    for (int i = 0; i < order; i++) {
        double xi = rule::abscissa()[i];
        double w1 = rule::weights()[i];
        
        for (int j = 0; j < order; j++) {
            double eta = rule::abscissa()[j];
            double w2 = rule::weights()[j];
            double weight = w1 * w2 * abs_detJ;
            
            // 形函数 + 物理坐标
            double phi[4];
            double x, y;
            transformToPhysical(elem, xi, eta, x, y, phi);
            
            // 插值解和梯度
            double uq = 0.0;
            Eigen::Vector2d grad_uq = Eigen::Vector2d::Zero();
            for (int k = 0; k < 4; k++) {
                uq        += phi[k] * u_elem(k);
                grad_uq(0) += grad_phi(0,k) * u_elem(k);
                grad_uq(1) += grad_phi(1,k) * u_elem(k);
            }
            
            // 源项 f 和 df/du
            double fq     = prob.f(uq, x, y);
            double dfdu_q = prob.df_du(uq, x, y);
            
            // 载荷向量
            for (int k = 0; k < 4; k++) {
                Fe(k) += (grad_phi(0,k)*grad_uq(0) + grad_phi(1,k)*grad_uq(1) - fq * phi[k]) * weight;
            }
            
            // 刚度矩阵
            for (int k = 0; k < 4; k++) {
                for (int l = 0; l < 4; l++) {
                    double term = (grad_phi(0,k)*grad_phi(0,l) + grad_phi(1,k)*grad_phi(1,l)) 
                                 - phi[k] * phi[l] * dfdu_q;
                    Ke(k,l) += term * weight;
                }
            }
        }
    }
}

} // namespace quad_boost
} // namespace fem_2d
} // namespace pde2d
