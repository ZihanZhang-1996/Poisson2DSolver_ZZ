#ifndef NODE_H
#define NODE_H

#include <array>
#include <cmath>

namespace pde2d {
namespace mesh_2d {

// 网格中的节点类，表示二维坐标点 (x,y)
class Node {
private:
    int id_;     // 节点编号
    double x_;   // 节点的 x 坐标
    double y_;   // 节点的 y 坐标
    
public:
    // 构造函数：初始化节点编号和坐标
    Node(int id, double x, double y) : id_(id), x_(x), y_(y) {}
    
    // 获取节点编号
    int getId() const { return id_; }
    // 获取节点 x 坐标
    double getX() const { return x_; }
    // 获取节点 y 坐标
    double getY() const { return y_; }
    
    // 设置节点编号
    void setId(int id) { id_ = id; }
    // 设置节点 x 坐标
    void setX(double x) { x_ = x; }
    // 设置节点 y 坐标
    void setY(double y) { y_ = y; }
    
    // 计算当前节点到另一个节点的距离
    double distanceTo(const Node& other) const {
        double dx = x_ - other.x_;
        double dy = y_ - other.y_;
        return std::sqrt(dx*dx + dy*dy);
    }
};

} // namespace mesh_2d
} // namespace pde2d

#endif // NODE_H
