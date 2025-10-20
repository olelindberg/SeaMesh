#ifndef GEOTIFF_READER_H
#define GEOTIFF_READER_H

#include "cpl_conv.h"
#include "gdal_priv.h"
#include <iostream>
#include <vtkImageActor.h>
#include <vtkImageData.h>
#include <vtkImageMapper3D.h>
#include <vtkSmartPointer.h>

class GeoTiffReader {

public:
  GeoTiffReader() { GDALAllRegister(); }

  ~GeoTiffReader() {
    // Cleanup GDAL
    GDALDestroyDriverManager();
  }
  int read(vtkSmartPointer<vtkImageActor> &actor) {
    GDALAllRegister();

    const char  *filename = "/home/ole/Projects/SeaMesh/data/Klimadatastyrelsen/DanmarksDybdeModel/2024/ddm_50m.dybde.tiff";
    GDALDataset *dataset  = (GDALDataset *)GDALOpen(filename, GA_ReadOnly);
    if (dataset == nullptr) {
      std::cerr << "Failed to open GeoTIFF: " << filename << std::endl;
      return 1;
    }

    // Get raster size
    int sizex = dataset->GetRasterXSize();
    int sizey = dataset->GetRasterYSize();
    int bands = dataset->GetRasterCount();

    std::cout << "Size: " << sizex << " x " << sizey << ", bands: " << bands << std::endl;

    // Get GeoTransform (affine transform)
    double geoTransform[6];
    if (dataset->GetGeoTransform(geoTransform) == CE_None) {
      std::cout << "Origin = (" << geoTransform[0] << ", " << geoTransform[3] << ")\n";
      std::cout << "Pixel Size = (" << geoTransform[1] << ", " << geoTransform[5] << ")\n";
    }

    // Read first band
    GDALRasterBand *band = dataset->GetRasterBand(1);
    // std::vector<float> buffer(sizex * sizey);
    std::vector<float> pixels(sizex * sizey);

    auto status = band->RasterIO(GF_Read, 0, 0, sizex, sizey, pixels.data(), sizex, sizey, GDT_Float32, 0, 0);

    std::cout << "Sample value: " << pixels[0] << std::endl;

    auto image = vtkSmartPointer<vtkImageData>::New();

    image->SetDimensions(sizex, sizey, 1);
    image->SetSpacing(geoTransform[1], std::abs(geoTransform[5]), 1.0);
    image->SetOrigin(geoTransform[0], geoTransform[3], 0.0);
    image->AllocateScalars(VTK_FLOAT, 1);

    float *vtkPtr = static_cast<float *>(image->GetScalarPointer());
    std::copy(pixels.begin(), pixels.end(), vtkPtr);

    actor = vtkSmartPointer<vtkImageActor>::New();
    actor->GetMapper()->SetInputData(image);

    GDALClose(dataset);
    return 1;
  }
};

#endif // GEOTIFF_READER_H