#include "../../BoostGeometryRtree.h"

#include <iostream>

int main() {

  // polygons
  auto polygons = std::make_shared<std::vector<polygon_t>>();

  // create some polygons
  for (unsigned i = 0; i < 10; ++i) {
    // create a polygon
    polygon_t p;
    for (float a = 0; a < 6.28316f; a += 1.04720f) {
      float x = i + int(10 * ::cos(a)) * 0.1f;
      float y = i + int(10 * ::sin(a)) * 0.1f;
      p.outer().push_back(point_t(x, y));
    }

    // add polygon
    polygons->push_back(p);
  }

  point_t point(10, 10);

  BoostRtreeSearch<polygon_t> search_tree(*polygons);
  auto result = search_tree.nearest(point);

  // display polygons
  std::cout << "generated polygons:" << std::endl;
  for (auto const &p : *polygons)
    std::cout << boost::geometry::wkt<polygon_t>(p) << std::endl;

  std::cout << "knn query point:" << std::endl;
  std::cout << boost::geometry::wkt<point_t>(point) << std::endl;

  std::cout << "knn query result:" << std::endl;
  std::cout << boost::geometry::wkt<polygon_t>(result) << std::endl;

  return 0;
}