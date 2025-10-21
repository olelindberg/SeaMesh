#include "QuadTreeSeaBuilder.h"
#include "BoostGeometryUtil.h"
#include "QuadTreeSearch.h"
#include <Eigen/Dense>

#include <iostream>

auto box_to_polygon = [](const box_t &b, polygon_t &box_poly)
{
  auto min = b.min_corner();
  auto max = b.max_corner();

  auto &outer = box_poly.outer();

  outer[0] = point_t(min.get<0>(), min.get<1>());
  outer[1] = point_t(min.get<0>(), max.get<1>());
  outer[2] = point_t(max.get<0>(), max.get<1>());
  outer[3] = point_t(max.get<0>(), min.get<1>());
  outer[4] = outer[0]; // close ring
};

QuadTreeSeaBuilder::QuadTreeSeaBuilder(int levelmax, double lengthmin, const std::shared_ptr<multi_polygon_t> &land_polygons)
    : _levelmax(levelmax), _lengthmin(lengthmin), _land_polygons(land_polygons)
{

  BoostGeomtryUtil::segment_rtree(*_land_polygons, _segment_rtree);

  _box_polygon.outer().resize(5);
}

void QuadTreeSeaBuilder::create(std::shared_ptr<QuadTree> &tree)
{

  if (tree->get_level() < _levelmax && tree->lengthx() > _lengthmin && tree->lengthy() > _lengthmin)
  {

    box_t box(point_t(tree->get_xmin(), tree->get_ymin()), point_t(tree->get_xmax(), tree->get_ymax()));
    box_to_polygon(box, _box_polygon);

    std::vector<segment_value_t> candidates;
    _segment_rtree.query(bgi::intersects(box), std::back_inserter(candidates));
    if (candidates.empty())
      return;

    this->subdivide(tree);

    this->create(tree->trees[0]);
    this->create(tree->trees[1]);
    this->create(tree->trees[2]);
    this->create(tree->trees[3]);
  }
}
