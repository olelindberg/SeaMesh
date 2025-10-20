#include "QuadTreeSeaBuilder.h"
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

  for (size_t p = 0; p < _land_polygons->size(); ++p)
  {
    auto const &poly = (*_land_polygons)[p];

    // Outer ring
    auto const &outer = poly.outer();
    for (size_t i = 0; i + 1 < outer.size(); ++i)
    {
      segment_rtree.insert({segment_t(outer[i], outer[i + 1]), {p, 0, i}});
    }

    // Inner rings
    for (size_t r = 0; r < poly.inners().size(); ++r)
    {
      auto const &inner = poly.inners()[r];
      for (size_t i = 0; i + 1 < inner.size(); ++i)
      {
        segment_rtree.insert({segment_t(inner[i], inner[i + 1]), {p, r + 1, i}});
      }
    }
  }

  _box_polygon.outer().resize(5);

  _search_tree = std::make_shared<BoostRtreeSearch<polygon_t>>((*_land_polygons));
}

void QuadTreeSeaBuilder::create(std::shared_ptr<QuadTree> &tree)
{

  if (tree->get_level() < _levelmax && tree->lengthx() > _lengthmin && tree->lengthy() > _lengthmin)
  {

    box_t box(point_t(tree->get_xmin(), tree->get_ymin()), point_t(tree->get_xmax(), tree->get_ymax()));
    box_to_polygon(box, _box_polygon);

    std::vector<segment_value_t> candidates;
    segment_rtree.query(bgi::intersects(box), std::back_inserter(candidates));
    if (candidates.empty())
      return;
    //    auto results = _search_tree->intersects(box);
    //    if (results.empty())
    //      return;
    //    multi_polygon_t polygons;
    //    for (const auto &res : results)
    //      polygons.push_back((*_land_polygons)[res.second]);
    //
    //    if (bg::disjoint(_box_polygon, polygons) || bg::within(_box_polygon, polygons))
    //      return;

    //    if (bg::disjoint(_box_polygon, *_land_polygons) || bg::within(_box_polygon, *_land_polygons))
    //      return;

    this->subdivide(tree);

    this->create(tree->trees[0]);
    this->create(tree->trees[1]);
    this->create(tree->trees[2]);
    this->create(tree->trees[3]);
  }
}
