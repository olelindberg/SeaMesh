#ifndef QUAD_TREE_POINTS_INSERTER_H
#define QUAD_TREE_POINTS_INSERTER_H

#include "BoostGeometryTypes.h"
#include "QuadTree.h"
#include <memory>

class QuadTreePointsInserter {
public:
  QuadTreePointsInserter(int levelmax, double lengthmin);
  void create(std::shared_ptr<QuadTree> &tree, point_t point);
  void balance(std::shared_ptr<QuadTree> &tree);
  int _levelmax;

  void subdivide(std::shared_ptr<QuadTree> &tree);

private:
  void create(std::shared_ptr<QuadTree> &tree, int level, int levelmax,
              point_t point);

  double _lengthmin;
};

#endif // QUAD_TREE_POINTS_INSERTER_H