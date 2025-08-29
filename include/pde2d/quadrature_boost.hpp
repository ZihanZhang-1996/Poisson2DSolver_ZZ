#ifndef PDE2D_QUADRATURE_BOOST_HPP_INCLUDED
#define PDE2D_QUADRATURE_BOOST_HPP_INCLUDED

#include <Eigen/Dense>
#include <array>
#include <vector>
#include <cmath>
#include <cassert>
#include <boost/math/quadrature/gauss.hpp>

#include "pde2d/Mesh2D.hpp"
#include "pde2d/problem_syb.hpp"

namespace pde2d {
namespace fem_2d {
namespace quad_boost {

// -----------------------------
// 单元装配：给定当前解 u，装配 Q4 单元的 Re(4) 和 Ke(4x4)
// e: 单元索引
// -----------------------------

void computeElementStiffnessAndForcingOnTriElem_boost(
    const std::shared_ptr<::pde2d::mesh_2d::Element> elem,
    const Eigen::VectorXd& u,
    Eigen::Matrix3d& Ke,
    Eigen::Vector3d& Fe,
    const ::pde2d::problem::ProblemExpr2D& prob
);

void computeElementStiffnessAndForcingOnRecElem_boost(
    const std::shared_ptr<::pde2d::mesh_2d::Element> elem,
    const Eigen::VectorXd& u,
    Eigen::Matrix4d& Ke,
    Eigen::Vector4d& Re,
    const ::pde2d::problem::ProblemExpr2D& prob);

} // namespace quad_boost
} // namespace fem_2d
} // namespace pde2d

#endif // PDE2D_QUADRATURE_BOOST_HPP_INCLUDED
