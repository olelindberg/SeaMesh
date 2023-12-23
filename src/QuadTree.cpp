#include "QuadTree.h"

#include <cmath>
#include <limits>

QuadTree::QuadTree(Eigen::Vector3d lower_, Eigen::Vector3d upper_, int lvl,
                   int index_)
    : lower(lower_), upper(upper_), index(index_) {
  level = lvl;
}

void QuadTree::center(Eigen::Vector3d &c) {
  c(0) = 0.5 * (lower(0) + upper(0));
  c(1) = 0.5 * (lower(1) + upper(1));
}

double QuadTree::lengthx() { return upper(0) - lower(0); }
double QuadTree::lengthy() { return upper(1) - lower(1); }

//--------------------------------------------------------------------------------
// Insert leaf nodes in list:
//--------------------------------------------------------------------------------
void QuadTree::insert_leaf(std::shared_ptr<QuadTree> q,
                           std::vector<std::shared_ptr<QuadTree>> &quad_list) {
  if (q->trees.size() == 0)
    quad_list.push_back(q);
  else
    for (std::vector<std::shared_ptr<QuadTree>>::iterator qi = q->trees.begin();
         qi != q->trees.end(); qi++)
      q->insert_leaf(*qi, quad_list);
}
