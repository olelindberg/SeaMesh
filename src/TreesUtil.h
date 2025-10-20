#ifndef TREES_UTIL_H
#define TREES_UTIL_H

#include "OcTree.h"
#include "QuadTree.h"

class TreesUtil {
public:
  static bool createSurfaceMesh(const OcTree::ptr &octree, QuadTree::ptr &quadtree);
  static bool createBottomMesh(const OcTree::ptr &octree, QuadTree::ptr &quadtree);

private:
  static bool _createSurfaceMesh(const OcTree::ptr &octree, QuadTree::ptr &quadtree);
  static bool _createBottomMesh(const OcTree::ptr &octree, QuadTree::ptr &quadtree);
};

#endif // TREES_UTIL_H