#ifndef QUAD_TREE_SEA_BUILDER_H
#define QUAD_TREE_SEA_BUILDER_H

#include "BoostGeometryRtree.h"
#include "BoostGeometryTypes.h"
#include "QuadTree.h"
#include "QuadTreeBuilderBase.h"
#include <memory>

class QuadTreeSeaBuilder : public QuadTreeBuilderBase
{
public:
  QuadTreeSeaBuilder(int levelmax, double lengthmin, const std::shared_ptr<multi_polygon_t> &land_polygons);
  void create(std::shared_ptr<QuadTree> &tree);

private:
  void create(std::shared_ptr<QuadTree> &tree, int level, int levelmax, point_t point);

  int                              _levelmax;
  double                           _lengthmin;
  std::shared_ptr<multi_polygon_t> _land_polygons;
  polygon_t                        _box_polygon;
  segment_rtree_t                  _segment_rtree;
};

#endif // QUAD_TREE_SEA_BUILDER_H