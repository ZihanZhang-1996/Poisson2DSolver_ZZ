#ifndef RECTANGULAR_ELEMENT_HPP
#define RECTANGULAR_ELEMENT_HPP

#include "Element.hpp"
#include <cmath>

namespace pde2d {
namespace mesh_2d {

// 矩形单元类，继承自 Element
// 表示有限元网格中的一个四边形单元（Q4 元素）
class RectangularElement : public Element {
public:
    // 构造函数：指定单元 ID，并设置类型为矩形
    RectangularElement(int id) : Element(id, ElementType::RECTANGULAR) {}
    
    // 返回单元的节点数（矩形单元固定为 4）
    int getNumNodes() const override { return 4; }
    
    // 计算矩形单元的面积
    double area() const override {
        if (nodes_.size() != 4) {
            throw std::runtime_error("Rectangular element requires exactly 4 nodes");
        }

        // 使用多边形面积公式（Shoelace Formula）
        double area = 0.0;
        for (int i = 0; i < 4; ++i) {
            int j = (i + 1) % 4;  // 下一个顶点
            area += nodes_[i]->getX() * nodes_[j]->getY();
            area -= nodes_[j]->getX() * nodes_[i]->getY();
        }
        return std::abs(area) / 2.0;
    }
};

} // namespace mesh_2d
} // namespace pde2d

#endif // RECTANGULAR_ELEMENT_HPP
