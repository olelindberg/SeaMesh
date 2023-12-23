#ifndef QUAD_TREE_BUILDER_H
#define QUAD_TREE_BUILDER_H

#include "QuadTree.h"
#include <memory>

class QuadTreeBuilder {
public:
  QuadTreeBuilder(int levelmax, double lengthmin);
  void create(std::shared_ptr<QuadTree> &tree);
  void balance(std::shared_ptr<QuadTree> &tree);
  int _levelmax;

  void subdivide(std::shared_ptr<QuadTree> &tree);

private:
  void create(std::shared_ptr<QuadTree> &tree, int level, int levelmax);

  double _lengthmin;
};

#endif // QUAD_TREE_BUILDER_H