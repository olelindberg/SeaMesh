#ifndef GEOTIFF_READER_H
#define GEOTIFF_READER_H

#include "geospatial_util.h"
#include "logger.h"

#include "cpl_conv.h"
#include "gdal_priv.h"

#include <iostream>

class GeoTiffReader
{

public:
  GeoTiffReader() { GDALAllRegister(); }

  ~GeoTiffReader()
  {
    // Cleanup GDAL
    GDALDestroyDriverManager();
  }

  int read()
  {
    GDALAllRegister();

    const char  *filename = "/home/ole/Projects/SeaMesh/data/Klimadatastyrelsen/DanmarksDybdeModel/2024/ddm_50m.dybde.tiff";
    GDALDataset *dataset  = (GDALDataset *)GDALOpen(filename, GA_ReadOnly);
    if (dataset == nullptr)
    {
      std::cerr << "Failed to open GeoTIFF: " << filename << std::endl;
      return 1;
    }
    Logger::info("Geotiff reader: Opened GeoTIFF: " + std::string(filename));

    // Get raster size
    int sizex = dataset->GetRasterXSize();
    int sizey = dataset->GetRasterYSize();
    int bands = dataset->GetRasterCount();

    Logger::info("Geotiff reader: GeoTIFF raster size: " + std::to_string(sizex) + " x " + std::to_string(sizey) + ", bands: " + std::to_string(bands));

    // --- Get projection reference (WKT string)
    auto projection_info = std::string(dataset->GetProjectionRef());
    GeospatialUtil::print_projection("Geotiff reader: ", projection_info);

    // Get GeoTransform (affine transform)
    double geoTransform[6];
    if (dataset->GetGeoTransform(geoTransform) == CE_None)
    {
      Logger::info("Geotiff reader: GeoTransform retrieved successfully.");
      Logger::info("Geotiff reader:   Top Left X:   " + std::to_string(geoTransform[0]));
      Logger::info("Geotiff reader:   Top Left Y:   " + std::to_string(geoTransform[3]));
      Logger::info("Geotiff reader:   Pixel Size X: " + std::to_string(geoTransform[1]));
      Logger::info("Geotiff reader:   Pixel Size Y: " + std::to_string(geoTransform[5]));
      Logger::info("Geotiff reader:   Rotation X:   " + std::to_string(geoTransform[2]));
      Logger::info("Geotiff reader:   Rotation Y:   " + std::to_string(geoTransform[4]));
    }

    // Read first band
    GDALRasterBand *band = dataset->GetRasterBand(1);
    // std::vector<float> buffer(sizex * sizey);
    std::vector<float> pixels(sizex * sizey);

    auto status = band->RasterIO(GF_Read, 0, 0, sizex, sizey, pixels.data(), sizex, sizey, GDT_Float32, 0, 0);

    GDALClose(dataset);

    if (status != CE_None)
    {
      Logger::error("Geotiff reader: Failed to read raster data from GeoTIFF.");
      return 1;
    }

    return 0;
  }
};

#endif // GEOTIFF_READER_H