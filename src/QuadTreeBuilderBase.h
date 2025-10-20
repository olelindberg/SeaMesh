#ifndef QUAD_TREE_BUILDER_BASE_H
#define QUAD_TREE_BUILDER_BASE_H

#include "BoostGeometryTypes.h"
#include "QuadTree.h"
#include <memory>

class QuadTreeBuilderBase {
public:
  void balance(std::shared_ptr<QuadTree> &tree);
  void subdivide(std::shared_ptr<QuadTree> &tree);
};

#endif // QUAD_TREE_BUILDER_BASE_H