#ifndef OCTREE_NODE_H
#define OCTREE_NODE_H

#include "Morton3D.h"

#include <fstream>

#include <cstdint>
#include <iostream>
#include <memory>
#include <vector>

// ============================================================================
// Node structure with spatial extent
// ============================================================================
struct OctreeNode
{
  uint64_t                                 morton;
  int                                      level;
  double                                   xmin, xmax;
  double                                   ymin, ymax;
  double                                   zmin, zmax;
  std::vector<std::unique_ptr<OctreeNode>> children;

  OctreeNode(uint64_t m, int l, double x0, double x1, double y0, double y1, double z0, double z1) : morton(m), level(l), xmin(x0), xmax(x1), ymin(y0), ymax(y1), zmin(z0), zmax(z1) {}
};

bool should_refine(const OctreeNode *node)
{
  return node->level < 3 && (node->morton % 3 == 0); // example condition
}

#endif // OCTREE_NODE_H