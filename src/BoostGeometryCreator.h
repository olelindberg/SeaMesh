#ifndef BOOST_GEOMETRY_CREATOR_H
#define BOOST_GEOMETRY_CREATOR_H

#include "BoostGeometryTypes.h"

#include <Eigen/Dense>

class BoostGeomtryCreator {
public:
  static std::shared_ptr<polygon_t> squarePolygon(double x0, double y0,
                                                  double halfwidth) {
    std::cout << "Creating rectangular inner boost polygon ...\n";
    auto polygon = std::make_shared<polygon_t>();
    polygon->inners().resize(1);
    boost::geometry::append(polygon->inners().back(),
                            point_t(x0 - halfwidth, y0 - halfwidth));
    boost::geometry::append(polygon->inners().back(),
                            point_t(x0 + halfwidth, y0 - halfwidth));
    boost::geometry::append(polygon->inners().back(),
                            point_t(x0 + halfwidth, y0 + halfwidth));
    boost::geometry::append(polygon->inners().back(),
                            point_t(x0 - halfwidth, y0 + halfwidth));
    boost::geometry::append(polygon->inners().back(),
                            point_t(x0 - halfwidth, y0 - halfwidth));
    return polygon;
  }

  static std::shared_ptr<polygon_t> circularPolygon(double x0, double y0,
                                                    double radius,
                                                    int numberSegments,
                                                    bool isInner = true) {
    std::cout << "Creating circular inner boost polygon ...\n";
    auto polygon = std::make_shared<polygon_t>();
    if (isInner)
      polygon->inners().resize(1);

    auto dtheta = 2.0 * M_PI / double(numberSegments);

    for (int i = 0; i < numberSegments + 1; ++i) {
      auto x = x0 + radius * std::cos(dtheta * i);
      auto y = y0 + radius * std::sin(dtheta * i);
      if (isInner)
        boost::geometry::append(polygon->inners().back(), point_t(x, y));
      else
        boost::geometry::append(polygon->outer(), point_t(x, y));
    }

    return polygon;
  }
};

#endif // BOOST_GEOMETRY_CREATOR_H
