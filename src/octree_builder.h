#ifndef OCTREE_BUILDER_H
#define OCTREE_BUILDER_H

#include "OctreeNode.h"
#include <iostream>

// ============================================================================
// Recursive builder
// ============================================================================
class OctreeBuilder
{
public:
  OctreeBuilder(int max_level) : _max_level(max_level) {}

  virtual bool should_refine(const OctreeNode *node)
  {
    std::cout << "Checking refinement for node level " << node->level << " morton " << node->morton << "\n";
    // Default refinement criterion (can be overridden in derived classes)
    return node->level < 3 && (node->morton % 3 == 0); // example condition
  }

  virtual void build_tree(OctreeNode *node)
  {

    if (node->level >= _max_level || !should_refine(node))
      return;

    double xm = 0.5 * (node->xmin + node->xmax);
    double ym = 0.5 * (node->ymin + node->ymax);
    double zm = 0.5 * (node->zmin + node->zmax);

    for (int i = 0; i < 8; ++i)
    {
      // child id bits: (xbit,ybit,zbit)
      int xb = (i >> 0) & 1;
      int yb = (i >> 1) & 1;
      int zb = (i >> 2) & 1;

      //-----------------------------------------------------------------------
      // Compute child bounding box:
      //-----------------------------------------------------------------------

      // Lower corner
      double x0 = xb ? xm : node->xmin;
      double y0 = yb ? ym : node->ymin;
      double z0 = zb ? zm : node->zmin;

      // Upper corner
      double x1 = xb ? node->xmax : xm;
      double y1 = yb ? node->ymax : ym;
      double z1 = zb ? node->zmax : zm;

      uint64_t child_code = (node->morton << 3) | i;

      node->children.emplace_back(std::make_unique<OctreeNode>(child_code, node->level + 1, x0, x1, y0, y1, z0, z1));

      build_tree(node->children.back().get());
    }
  }

protected:
  int _max_level;
};
#endif // OCTREE_BUILDER_H
