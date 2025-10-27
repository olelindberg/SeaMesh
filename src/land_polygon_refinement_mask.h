#ifndef LAND_POLYGON_REFINEMENT_MASK_H
#define LAND_POLYGON_REFINEMENT_MASK_H

#include "directional_adaptive_octree.h"
#include "refine_mask_selector.h"

class LandPolygonRefineMask : public RefineMaskSelector
{
public:
  LandPolygonRefineMask(int levelmax, segment_rtree_t &segment_rtree) : _levelmax(levelmax), _segment_rtree(segment_rtree) {}

  virtual RefineMask operator()(const DirectionalAdaptiveOctree *node, double cx, double cy, double cz, int level_x, int level_y, int level_z) const override
  {
    if (level_x >= _levelmax || level_y >= _levelmax || level_z >= _levelmax)
      return RefineMask::REFINE_NONE;

    box_t box(point_t(node->xmin, node->ymin), point_t(node->xmax, node->ymax));

    double center_x = 0.5 * (node->xmin + node->xmax);

    std::vector<segment_value_t> candidates;
    _segment_rtree.query(bgi::intersects(box), std::back_inserter(candidates));

    if (!candidates.empty())
      return RefineMask::REFINE_X | RefineMask::REFINE_Y;

    return RefineMask::REFINE_NONE;
  }

private:
  int              _levelmax;
  segment_rtree_t &_segment_rtree;
};

#endif // LAND_POLYGON_REFINEMENT_MASK_H