#ifndef MESH2D_HPP_INCLUDED
#define MESH2D_HPP_INCLUDED

#include "Element.hpp"
#include "RectangularElement.hpp"
#include "TriangularElement.hpp"
#include <vector>
#include <memory>
#include <functional>

namespace pde2d {
namespace mesh_2d {

// 二维结构化网格类，支持矩形单元和三角形单元
class Mesh2D {
private:
    // 网格的所有节点
    std::vector<std::shared_ptr<Node>> nodes_;
    // 网格的所有单元
    std::vector<std::shared_ptr<Element>> elements_;
    // 网格在 x 和 y 方向的划分数
    int nx_, ny_;
    // 物理区域的长和宽
    double Lx_, Ly_;
    // 网格间距
    double dx_, dy_;
    // 单元类型（矩形或三角形）
    ElementType elementType_;
    
public:
    Mesh2D() {}

    // 构造函数：根据区域大小、划分数和单元类型生成网格
    Mesh2D(double Lx, double Ly, int nx, int ny, ElementType type)
        : nx_(nx), ny_(ny)
        , Lx_(Lx), Ly_(Ly)
        , dx_(Lx/nx), dy_(Ly/ny)
        , elementType_(type) {
        createNodes();
        createElements();
    }

    const std::vector<std::shared_ptr<Node>>& getNodes() const { return nodes_; }
    const std::vector<std::shared_ptr<Element>>& getElements() const { return elements_; }
    int getnx() const {return nx_;}
    int getny() const {return ny_;}
    double getLx() const {return Lx_;}
    double getLy() const {return Ly_;}
    int getNumNodes() const { return nodes_.size(); }
    int getNumElements() const { return elements_.size(); }
    ElementType getElementType() const { return elementType_; }

    // 创建网格节点
    void createNodes();
    
    // 创建网格单元
    void createElements();
    
    // 生成矩形单元
    void createRectangularElements(int& elemId);
    
    // 生成三角形单元
    void createTriangularElements(int& elemId);
    
    // 查找包含某一点 (x,y) 的单元
    std::shared_ptr<Element> findElementContaining(double x, double y) const;
    
private:
    // 判断点 (x,y) 是否在某个单元内（box判断）
    bool isPointInElement(const std::shared_ptr<Element>& element, double x, double y) const;
};

} // namespace mesh_2d
} // namespace pde2d

#endif // MESH2D_HPP_INCLUDED
