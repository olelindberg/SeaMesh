#ifndef OCTREE_BUILDER_LAND_POLYGON_H
#define OCTREE_BUILDER_LAND_POLYGON_H

#include "BoostGeometryTypes.h"
#include "BoostGeometryUtil.h"
#include "octree_builder_xy.h"
#include <iostream>

class OctreeBuilderLandPolygon : public OctreeBuilderXY
{
public:
  OctreeBuilderLandPolygon(int max_level, const multi_polygon_t &land_polygons) : OctreeBuilderXY(max_level) { BoostGeomtryUtil::segment_rtree(land_polygons, _segment_rtree); }

  virtual bool should_refine(const OctreeNode *node) override
  {
    box_t box(point_t(node->xmin, node->ymin), point_t(node->xmax, node->ymax));

    std::vector<segment_value_t> candidates;
    _segment_rtree.query(bgi::intersects(box), std::back_inserter(candidates));

    return !candidates.empty();
  }

private:
  segment_rtree_t _segment_rtree;
};

#endif // OCTREE_BUILDER_LAND_POLYGON_H