#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include "Node.hpp"
#include <vector>
#include <memory>
#include <stdexcept>

namespace pde2d {
namespace mesh_2d {

/// 元素类型：矩形单元 or 三角形单元
enum class ElementType { RECTANGULAR, TRIANGULAR };

/// 有限元网格中的单元基类
/// - 每个单元有一个 ID
/// - 一个类型（矩形 / 三角形）
/// - 一组节点（shared_ptr<Node>）
class Element {
protected:
    int id_;                                   // 单元 ID
    ElementType type_;                         // 单元类型（矩形/三角形）
    std::vector<std::shared_ptr<Node>> nodes_; // 单元包含的节点

public:
    /// 构造函数
    Element(int id, ElementType type) : id_(id), type_(type) {}

    /// 析构函数
    virtual ~Element() = default;
    
    // === 基本信息访问接口 ===
    int getId() const { return id_; }                       // 获取单元 ID
    ElementType getType() const { return type_; }           // 获取单元类型
    const std::vector<std::shared_ptr<Node>>& getNodes() const { return nodes_; } // 获取节点列表
    
    // === 抽象接口 ===
    // 每个子类（矩形单元/三角形单元）必须实现以下接口：
    virtual double area() const = 0;        // 单元面积
    virtual int getNumNodes() const = 0;    // 单元的节点个数（矩形=4，三角形=3）
    
    // === 节点操作 ===

    /// 向单元中添加一个节点
    /// - 检查是否超过单元允许的节点数
    void addNode(std::shared_ptr<Node> node) {
        if (nodes_.size() >= static_cast<size_t>(getNumNodes())) {
            throw std::runtime_error("Cannot add more nodes to element"); // 超出节点数时报错
        }
        nodes_.push_back(node);
    }
    
    /// 根据局部索引获取节点
    /// - localIndex：0,1,...,(getNumNodes()-1)
    /// - 越界则抛出异常
    std::shared_ptr<Node> getNode(int localIndex) const {
        if (localIndex < 0 || localIndex >= static_cast<int>(nodes_.size())) {
            throw std::out_of_range("Invalid node index");
        }
        return nodes_[localIndex];
    }
};

} // namespace mesh_2d
} // namespace pde2d

#endif // ELEMENT_HPP
