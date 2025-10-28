#ifndef BOOST_MULTIPOLYGON_VTK_WRITER_H
#define BOOST_MULTIPOLYGON_VTK_WRITER_H

#include "BoostGeometryTypes.h"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>

#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyLine.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>

#include <fstream>

class BoostMultiPolygonVtkWriter
{
public:
  //  namespace bg          = boost::geometry;
  //  using point_t         = bg::model::point<double, 2, bg::cs::cartesian>;
  //  using polygon_t       = bg::model::polygon<point_t, true>;
  //  using multi_polygon_t = bg::model::multi_polygon<polygon_t>;

  static vtkSmartPointer<vtkPolyData> multi_polygon_to_lines(const multi_polygon_t &mp)
  {
    auto vtk_points = vtkSmartPointer<vtkPoints>::New();
    auto vtk_lines  = vtkSmartPointer<vtkCellArray>::New();

    for (const auto &poly : mp)
    {
      // Outer + all hole rings
      std::vector<const polygon_t::ring_type *> rings;
      // Push outer ring
      rings.push_back(&poly.outer());
      // Push inner rings
      for (auto const &r : poly.inners())
        rings.push_back(&r);

      for (const auto *ring : rings)
      {
        if (ring->size() < 2)
          continue;

        auto      polyLine = vtkSmartPointer<vtkPolyLine>::New();
        vtkIdType n        = static_cast<vtkIdType>(ring->size());

        polyLine->GetPointIds()->SetNumberOfIds(n);

        for (vtkIdType i = 0; i < n; ++i)
        {
          const auto &pt = (*ring)[i];
          double      x  = bg::get<0>(pt);
          double      y  = bg::get<1>(pt);

          vtkIdType id = vtk_points->InsertNextPoint(x, y, 0.0);
          polyLine->GetPointIds()->SetId(i, id);
        }

        vtk_lines->InsertNextCell(polyLine);
      }
    }

    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(vtk_points);
    polyData->SetLines(vtk_lines);
    return polyData;
  }

  static void write_vtp_lines(vtkPolyData *pd, const std::string &filename)
  {
    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(pd);     // ✅ correct modern API
    writer->SetDataModeToAscii(); // ✅ easier debugging first
    writer->Write();
  }
};

#endif // BOOST_MULTIPOLYGON_VTK_WRITER_H