#ifndef TEST_GPKG_TO_BOOST_MULTIPOLYGON_H
#define TEST_GPKG_TO_BOOST_MULTIPOLYGON_H

#include "BoostGeometryTypes.h"

// gpkg_to_boost_multipolygon.cpp
#include <gdal_priv.h>
#include <iostream>
#include <memory>
#include <ogrsf_frmts.h>
#include <vector>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>

// namespace bg = boost::geometry;
//
//// ---- Boost types (adjust if you already have your own) ----
// using point_t         = bg::model::d2::point_xy<double>;
// using ring_t          = bg::model::ring<point_t, /*ClockWise?*/ false, /*Closed?*/ false>;
// using polygon_t       = bg::model::polygon<point_t, /*CCW outer*/ false, /*Open*/ false>;
// using multi_polygon_t = bg::model::multi_polygon<polygon_t>;

// Read an OGRLinearRing into a Boost ring (skip duplicate closing point; we'll let bg::correct close it)
static void ogr_linear_ring_to_boost_ring(const OGRLinearRing *lr, ring_t &ring)
{
  if (!lr)
    return;
  int n = lr->getNumPoints();
  // OGR rings are closed (last == first). Skip the last duplicate point.
  int limit = (n >= 2) ? n - 1 : n;
  ring.clear();
  ring.reserve(limit);
  for (int i = 0; i < limit; ++i)
    ring.emplace_back(lr->getX(i), lr->getY(i));
}

// Convert an OGRPolygon (2D) to a Boost polygon
static bool ogr_polygon_to_boost_polygon(const OGRPolygon *opoly, polygon_t &bp)
{
  if (!opoly)
    return false;

  // Exterior
  const OGRLinearRing *ext = opoly->getExteriorRing();
  if (!ext || ext->getNumPoints() < 4) // minimal closed tri -> 4 vertices incl. duplicate
    return false;

  ring_t outer;
  ogr_linear_ring_to_boost_ring(ext, outer);
  bp.outer().assign(outer.begin(), outer.end());

  // Holes
  bp.inners().clear();
  const int nh = opoly->getNumInteriorRings();
  bp.inners().reserve(nh);
  for (int i = 0; i < nh; ++i)
  {
    const OGRLinearRing *in = opoly->getInteriorRing(i);
    if (!in || in->getNumPoints() < 4)
      continue;
    ring_t inner;
    ogr_linear_ring_to_boost_ring(in, inner);
    bp.inners().emplace_back();
    bp.inners().back().assign(inner.begin(), inner.end());
  }

  // Fix orientation/closure to Boost expectations
  bg::correct(bp);
  return true;
}

// Recursively harvest polygons from any OGRGeometry (handles Polygon, MultiPolygon, GeometryCollection)
static void collect_polygons(const OGRGeometry *g, std::vector<const OGRPolygon *> &out)
{
  if (!g)
    return;
  OGRwkbGeometryType t = wkbFlatten(g->getGeometryType());

  if (t == wkbPolygon)
  {
    out.push_back(g->toPolygon());
    return;
  }
  if (t == wkbMultiPolygon)
  {
    const OGRMultiPolygon *mp = g->toMultiPolygon();
    for (int i = 0; i < mp->getNumGeometries(); ++i)
      out.push_back(mp->getGeometryRef(i)->toPolygon());
    return;
  }
  if (t == wkbGeometryCollection)
  {
    const OGRGeometryCollection *gc = g->toGeometryCollection();
    for (int i = 0; i < gc->getNumGeometries(); ++i)
      collect_polygons(gc->getGeometryRef(i), out);
    return;
  }
  // Ignore lines/points, etc.
}

// Load a layer from a .gpkg and return as Boost multi_polygon
// layer_name: pass nullptr to use the first (0) layer
// target_srs_wkt: optional; if non-null, transform coordinates to this SRS
bool gpkg_to_boost_multipolygon(const std::string &filename, std::string layer_name, std::string target_srs_wkt, multi_polygon_t &out_mpoly)
{
  GDALAllRegister();

  std::unique_ptr<GDALDataset> ds(static_cast<GDALDataset *>(GDALOpenEx(filename.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr)));
  if (!ds)
  {
    std::cerr << "Failed to open: " << filename << "\n";
    return false;
  }

  OGRLayer *layer = nullptr;
  if (layer_name.empty() == false)
  {
    layer = ds->GetLayerByName(layer_name.c_str());
    if (!layer)
    {
      std::cerr << "Layer not found: " << layer_name << "\n";
      return false;
    }
  }
  else
  {
    layer = ds->GetLayer(0);
    if (!layer)
    {
      std::cerr << "No layers in dataset.\n";
      return false;
    }
  }

  // Optional transformation
  std::unique_ptr<OGRCoordinateTransformation> coord_tx;
  if (target_srs_wkt.empty() == false)
  {
    OGRSpatialReference srcSRS = *layer->GetSpatialRef(); // may be null; handle below
    OGRSpatialReference dstSRS;
    if (dstSRS.SetFromUserInput(target_srs_wkt.c_str()) != OGRERR_NONE)
    {
      std::cerr << "Invalid target SRS WKT/def.\n";
      return false;
    }
    if (layer->GetSpatialRef())
    {
      coord_tx.reset(OGRCreateCoordinateTransformation(&srcSRS, &dstSRS));
      if (!coord_tx)
      {
        std::cerr << "Failed to create coordinate transformation.\n";
        return false;
      }
    }
    else
    {
      std::cerr << "Warning: source layer has no SRS; skipping transform.\n";
    }
  }

  out_mpoly.clear();

  layer->ResetReading();
  OGRFeature *feat = nullptr;
  while ((feat = layer->GetNextFeature()) != nullptr)
  {
    std::unique_ptr<OGRFeature> feat_guard(feat);
    OGRGeometry                *geom = feat->GetGeometryRef();
    if (!geom)
      continue;

    // Work on a 2D clone (drop Z/M if present)
    std::unique_ptr<OGRGeometry> g2d(geom->clone());
    g2d->flattenTo2D();

    if (coord_tx)
    {
      if (g2d->transform(coord_tx.get()) != OGRERR_NONE)
      {
        std::cerr << "Transform failed for a feature; skipping.\n";
        continue;
      }
    }

    std::vector<const OGRPolygon *> polys;
    collect_polygons(g2d.get(), polys);
    for (const OGRPolygon *op : polys)
    {
      polygon_t bp;
      if (ogr_polygon_to_boost_polygon(op, bp))
      {
        out_mpoly.push_back(std::move(bp));
      }
    }
  }

  // Normalize/correct the whole multi_polygon
  bg::correct(out_mpoly);
  return true;
}

// ---- Example usage ----
int test_gpkg_to_boost_multipolygon(std::string gpkg, multi_polygon_t &mp, std::string layer_name = "", std::string target_srs = "")
{
  //    std::cerr << "Usage: " << argv[0] << " landpolygon.gpkg [layer_name] [target_srs]\n"
  //              << "Example target_srs: EPSG:4326 or PROJ string or WKT\n";

  if (!gpkg_to_boost_multipolygon(gpkg, layer_name, target_srs, mp))
  {
    std::cerr << "Failed to build multi-polygon.\n";
    return 2;
  }

  // Quick sanity: print bbox and polygon count
  bg::model::box<point_t> box;
  bg::envelope(mp, box);
  std::cout << "Polygons: " << mp.size() << "\n";
  std::cout << "BBox: (" << box.min_corner().get<0>() << "," << box.min_corner().get<1>() << ") - (" << box.max_corner().get<0>() << ","
            << box.max_corner().get<1>() << ")\n";
  return 0;
}

#endif // TEST_GPKG_TO_BOOST_MULTIPOLYGON_H