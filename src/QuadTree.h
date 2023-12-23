#ifndef QUAD_TREE_H
#define QUAD_TREE_H

#include <Eigen/Dense>

#include <memory>
#include <vector>

class QuadTree : public std::enable_shared_from_this<QuadTree> {
public:
  typedef std::shared_ptr<QuadTree> ptr;

  QuadTree(Eigen::Vector3d lower_, Eigen::Vector3d upper_, int lvl, int index_);

  void balance();
  void center(Eigen::Vector3d &c);
  double lengthx();
  double lengthy();

  int get_depth() { return level; }
  const std::shared_ptr<QuadTree> &get_parent() { return _parent; };
  void set_parent(std::shared_ptr<QuadTree> &parent) { _parent = parent; };

  const std::shared_ptr<QuadTree> &get_child(int i) { return trees[i]; };

  double get_xmin() { return lower(0); }
  double get_xmax() { return upper(0); }
  double get_ymin() { return lower(1); }
  double get_ymax() { return upper(1); }

  int level;
  Eigen::Vector3d lower;
  Eigen::Vector3d upper;
  //  std::vector<Eigen::Vector3d> points;
  std::vector<std::shared_ptr<QuadTree>> trees;

  int index;

  int get_id() { return _id; }
  void set_id(int id) { _id = id; }

  bool flag = false;

private:
  int _id = -1;

  void insert_leaf(std::shared_ptr<QuadTree> q,
                   std::vector<std::shared_ptr<QuadTree>> &quad_list);

  std::shared_ptr<QuadTree> _parent;
};

#endif // QUAD_TREE_H
