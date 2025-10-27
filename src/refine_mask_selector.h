#ifndef REFINE_MASK_SELECTOR_H
#define REFINE_MASK_SELECTOR_H

#include "directional_adaptive_octree.h"
#include "refine_mask_util.h"

class RefineMaskSelector
{
public:
  virtual RefineMask operator()(const DirectionalAdaptiveOctree *node, double cx, double cy, double cz, int level_x, int level_y, int level_z) const = 0;
};

#endif // REFINE_MASK_SELECTOR_H