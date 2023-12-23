#ifndef BOOST_GEOMETRY_RTREE_H
#define BOOST_GEOMETRY_RTREE_H

#include "BoostGeometryTypes.h"

template <typename geometry_t> class BoostRtreeSearch {
public:
  BoostRtreeSearch(std::vector<geometry_t> &geometries)
      : _geometries(geometries) {

    for (unsigned i = 0; i < _geometries.size(); ++i)
      _rtree.insert(std::make_pair(
          boost::geometry::return_envelope<box_t>((_geometries)[i]), i));
  }

  const geometry_t &nearest(const point_t &point) {

    std::vector<value_t> result;
    _rtree.query(boost::geometry::index::nearest(point, 1),
                 std::back_inserter(result));

    return (_geometries)[result[0].second];
  }

private:
  boost::geometry::index::rtree<value_t, boost::geometry::index::rstar<16, 4>>
      _rtree;

  std::vector<geometry_t>
      &_geometries; // TODO, change reference to shared pointer
};
#endif // BOOST_GEOMETRY_RTREE_H