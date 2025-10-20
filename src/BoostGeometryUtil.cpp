#include "BoostGeometryUtil.h"

void BoostGeometryUtil::printMultiPolygon(const std::shared_ptr<multi_polygon_t> &multi_polygon) {
  int poly_idx = 0;
  for (auto const &poly : *multi_polygon) {
    std::cout << "Polygon " << poly_idx++ << ":\n";

    // Exterior ring
    std::cout << "  Exterior ring:\n";
    for (auto const &p : poly.outer()) {
      std::cout << "    (" << bg::get<0>(p) << ", " << bg::get<1>(p) << ")\n";
    }

    // Interior rings (holes)
    int hole_idx = 0;
    for (auto const &hole : poly.inners()) {
      std::cout << "  Hole " << hole_idx++ << ":\n";
      for (auto const &p : hole) {
        std::cout << "    (" << bg::get<0>(p) << ", " << bg::get<1>(p) << ")\n";
      }
    }
  }
}