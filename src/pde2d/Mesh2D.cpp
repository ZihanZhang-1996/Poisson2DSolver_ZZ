#include "pde2d/Mesh2D.hpp"
#include <cassert>

namespace pde2d {
namespace mesh_2d {

//
// 根据 nx, ny 创建结构化网格的节点
// 节点按照 (i,j) 坐标顺序生成，存储在 nodes_ 向量中
//
void Mesh2D::createNodes() {
    int nodeId = 0;
    for (int j = 0; j <= ny_; ++j) {
        for (int i = 0; i <= nx_; ++i) {
            double x = i * dx_; // 节点横坐标
            double y = j * dy_; // 节点纵坐标
            // 创建节点并放入 nodes_ 列表
            nodes_.push_back(std::make_shared<Node>(nodeId++, x, y));
        }
    }
}

//
// 根据 elementType_ 创建单元
// - RECTANGULAR: 四节点矩形单元
// - TRIANGULAR:  三节点三角形单元 (每个矩形切分为两个三角形)
//
void Mesh2D::createElements() {
    int elemId = 0;
    
    if (elementType_ == ElementType::RECTANGULAR) {
        createRectangularElements(elemId);
    } else {
        createTriangularElements(elemId);
    }
}

//
// 创建四节点矩形单元 (Q4)
// 每个矩形单元由 (i,j) 网格的 4 个顶点组成
//
void Mesh2D::createRectangularElements(int& elemId) {
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            auto element = std::make_shared<RectangularElement>(elemId++);
            
            // 计算该单元 4 个角点的全局节点编号
            int n1 = j * (nx_ + 1) + i;            // 左下角
            int n2 = j * (nx_ + 1) + i + 1;        // 右下角
            int n3 = (j + 1) * (nx_ + 1) + i + 1;  // 右上角
            int n4 = (j + 1) * (nx_ + 1) + i;      // 左上角
            
            // 添加节点到单元
            element->addNode(nodes_[n1]);
            element->addNode(nodes_[n2]);
            element->addNode(nodes_[n3]);
            element->addNode(nodes_[n4]);
            
            elements_.push_back(element);
        }
    }
}

//
// 创建三角形单元 (T3)
// 每个矩形单元被切分为两个三角形 (对角线划分)
//
void Mesh2D::createTriangularElements(int& elemId) {
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            // 当前小矩形的 4 个顶点编号
            int n1 = j * (nx_ + 1) + i;
            int n2 = j * (nx_ + 1) + i + 1;
            int n3 = (j + 1) * (nx_ + 1) + i + 1;
            int n4 = (j + 1) * (nx_ + 1) + i;
            
            // 第一个三角形 (n1, n2, n3)
            auto tri1 = std::make_shared<TriangularElement>(elemId++);
            tri1->addNode(nodes_[n1]);
            tri1->addNode(nodes_[n2]);
            tri1->addNode(nodes_[n3]);
            elements_.push_back(tri1);
            
            // 第二个三角形 (n1, n3, n4)
            auto tri2 = std::make_shared<TriangularElement>(elemId++);
            tri2->addNode(nodes_[n1]);
            tri2->addNode(nodes_[n3]);
            tri2->addNode(nodes_[n4]);
            elements_.push_back(tri2);
        }
    }
}

//
// 查找包含点 (x,y) 的单元
// 遍历所有单元，调用 isPointInElement 判断
//
std::shared_ptr<Element> 
Mesh2D::findElementContaining(double x, double y) const {
    for (const auto& element : elements_) {
        if (isPointInElement(element, x, y)) {
            return element;
        }
    }
    return nullptr;
}

//
// 判断点 (x,y) 是否在某个单元内
// 当前实现只是做一个简单的box

bool Mesh2D::isPointInElement(const std::shared_ptr<Element>& element, double x, double y) const {
    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::lowest();
    
    // 遍历单元节点，找出box
    for (const auto& node : element->getNodes()) {
        minX = std::min(minX, node->getX());
        maxX = std::max(maxX, node->getX());
        minY = std::min(minY, node->getY());
        maxY = std::max(maxY, node->getY());
    }
    
    // 判断点是否在box内
    return (x >= minX && x <= maxX && y >= minY && y <= maxY);
}

} // namespace mesh_2d
} // namespace pde2d
