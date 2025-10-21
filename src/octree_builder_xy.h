#ifndef OCTREE_BUILDER_XY_H
#define OCTREE_BUILDER_XY_H

#include "OctreeNode.h"
#include "octree_builder.h"

class OctreeBuilderXY : public OctreeBuilder
{
public:
  OctreeBuilderXY(int max_level) : OctreeBuilder(max_level) {}

  // Same OctreeNode struct you already have; use this build variant:
  virtual void build_tree(OctreeNode *node) override
  {
    if (node->level >= _max_level || !should_refine(node))
      return;

    double xm = 0.5 * (node->xmin + node->xmax);
    double ym = 0.5 * (node->ymin + node->ymax);
    // keep z split unchanged: use node->zmin/node->zmax for all children

    for (int i = 0; i < 4; ++i)
    {
      int xb = (i >> 0) & 1; // x-bit
      int yb = (i >> 1) & 1; // y-bit
      int zb = 0;            // no z split

      double x0 = xb ? xm : node->xmin;
      double x1 = xb ? node->xmax : xm;
      double y0 = yb ? ym : node->ymin;
      double y1 = yb ? node->ymax : ym;
      double z0 = node->zmin;
      double z1 = node->zmax;

      // use 2 bits per level (child id 0..3)
      uint64_t child_code = (node->morton << 2) | (uint64_t)i;
      node->children.emplace_back(std::make_unique<OctreeNode>(child_code, node->level + 1, x0, x1, y0, y1, z0, z1));

      build_tree(node->children.back().get()); // only recurse on half the children to create quadtree
    }
  }
};

#endif // OCTREE_BUILDER_XY_H