#ifndef TRIANGULARELEMENT_HPP
#define TRIANGULARELEMENT_HPP

#include "Element.hpp"
#include <cmath>

namespace pde2d {
namespace mesh_2d {

// 三角形单元类，继承自 Element
// 表示有限元网格中的一个三角形单元（T3 元素）
class TriangularElement : public Element {
public:
    // 构造函数：指定单元 ID，并设置类型为三角形
    TriangularElement(int id) : Element(id, ElementType::TRIANGULAR) {}
    
    // 返回单元的节点数（三角形单元固定为 3）
    int getNumNodes() const override { return 3; }
    
    // 计算三角形单元的面积
    double area() const override {
        if (nodes_.size() != 3) {
            throw std::runtime_error("Triangular element requires exactly 3 nodes");
        }
        
        // 提取三个节点的坐标
        double x1 = nodes_[0]->getX(), y1 = nodes_[0]->getY();
        double x2 = nodes_[1]->getX(), y2 = nodes_[1]->getY();
        double x3 = nodes_[2]->getX(), y3 = nodes_[2]->getY();
        
        // 使用三角形面积公式（叉积/行列式公式）
        return std::abs((x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)) / 2.0);
    }
};

} // namespace mesh_2d
} // namespace pde2d

#endif // TRIANGULARELEMENT_HPP
