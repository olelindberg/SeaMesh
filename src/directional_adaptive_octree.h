#ifndef DIRECTIONAL_ADAPTIVE_OCTREE_H
#define DIRECTIONAL_ADAPTIVE_OCTREE_H

#include "octree_refine_mask_util.h"
#include <cstdint>
#include <memory>
#include <vector>

// ------------------------- Octree node -------------------------
struct DirectionalAdaptiveOctree
{
  uint64_t                                                morton;
  int                                                     level_x, level_y, level_z;
  double                                                  xmin, xmax, ymin, ymax, zmin, zmax;
  RefineMask                                              mask; // which axes to refine when subdividing this node
  std::vector<std::unique_ptr<DirectionalAdaptiveOctree>> children;
  DirectionalAdaptiveOctree                              *parent = nullptr; // <-- new

  DirectionalAdaptiveOctree(uint64_t m, int xl, int yl, int zl, double x0, double x1, double y0, double y1, double z0, double z1,
                            RefineMask rm = RefineMask::REFINE_NONE, DirectionalAdaptiveOctree *p = nullptr)
      : morton(m), level_x(xl), level_y(yl), level_z(zl), xmin(x0), xmax(x1), ymin(y0), ymax(y1), zmin(z0), zmax(z1), mask(rm), parent(p)
  {
  }
  bool is_leaf() const { return children.empty(); }
};

#endif // DIRECTIONAL_ADAPTIVE_OCTREE_H
