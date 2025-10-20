#ifndef OC_TREE_H
#define OC_TREE_H

#include <Eigen/Dense>
#include <memory>
#include <vector>

class OcTree : public std::enable_shared_from_this<OcTree> {
public:
  enum CHILDREN { SOUTH_WEST_BOTTOM = 0, SOUTH_EAST_BOTTOM = 1, NORTH_WEST_BOTTOM = 2, NORTH_EAST_BOTTOM = 3, SOUTH_WEST_TOP = 4, SOUTH_EAST_TOP = 5, NORTH_WEST_TOP = 6, NORTH_EAST_TOP = 7 };

  typedef std::shared_ptr<OcTree> ptr;

  OcTree(Eigen::Vector3d lower_, Eigen::Vector3d upper_, int level_, int index_) : _lower(lower_), _upper(upper_), _level(level_), _index(index_) {}

  void center(Eigen::Vector3d &c) { c = 0.5 * (_lower + _upper); }
  void length(Eigen::Vector3d &l) { l = _upper - _lower; }

  const std::vector<std::shared_ptr<OcTree>> &childs() { return _childs; };
  const Eigen::Vector3d                      &lower() { return _lower; }
  const Eigen::Vector3d                      &upper() { return _upper; }
  int                                         index() { return _index; }
  int                                         level() { return _level; }

  void set_index(int index) { _index = index; }

private:
  int                                  _index;
  int                                  _level;
  Eigen::Vector3d                      _lower;
  Eigen::Vector3d                      _upper;
  std::vector<std::shared_ptr<OcTree>> _childs;
};

#endif // OC_TREE_H
