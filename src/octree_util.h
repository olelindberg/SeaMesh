#ifndef OCTREE_UTIL_H
#define OCTREE_UTIL_H

#include "Morton3D.h"
#include "OctreeNode.h"

#include <cstdint>
#include <iostream>

// ============================================================================
// Printing utility
// ============================================================================
void print_tree(const OctreeNode *node, int indent = 0)
{
  for (int i = 0; i < indent; ++i)
    std::cout << "  ";
  std::cout << "Level " << node->level << " Morton " << node->morton;

  uint32_t ix, iy, iz;
  Morton3D::grid_index(node->morton, ix, iy, iz);
  std::cout << "  Index(" << ix << "," << iy << "," << iz << ")";
  std::cout << "  Box[(" << node->xmin << "," << node->ymin << "," << node->zmin << ")-(" << node->xmax << "," << node->ymax << "," << node->zmax << ")]\n";

  for (const auto &c : node->children)
    print_tree(c.get(), indent + 1);
}

// ============================================================================
// Collect leaf cells
// ============================================================================
void collect_leaves(const OctreeNode *node, std::vector<const OctreeNode *> &leaves)
{
  if (node->children.empty())
  {
    leaves.push_back(node);
    return;
  }
  for (const auto &c : node->children)
    collect_leaves(c.get(), leaves);
}

// ------------------------- Build children according to node->mask -------------------------
auto compute_bounds = [](double minv, double maxv, int n, int idx)
{
  double mid = 0.5 * (minv + maxv);
  if (n == 1)
    return std::make_pair(minv, maxv);
  return idx == 0 ? std::make_pair(minv, mid) : std::make_pair(mid, maxv);
};

#endif // OCTREE_UTIL_H