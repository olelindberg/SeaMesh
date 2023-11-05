#include <shapefil.h>

#include <iostream>

int main() {

  const char *pszShapeFile =
      "../data/Topo-bathy/shape_format/data/Kystlinie.shp";

  SHPHandle shphandle = SHPOpen(pszShapeFile, "rb");

  int pnEntities;
  int pnShapeType;
  double padfMinBound[4];
  double padfMaxBound[4];
  SHPGetInfo(shphandle, &pnEntities, &pnShapeType, padfMinBound, padfMaxBound);

  std::cout << std::endl;
  std::cout << pnEntities << std::endl;
  std::cout << pnShapeType << std::endl;
  std::cout << padfMinBound[0] << std::endl;
  std::cout << padfMinBound[1] << std::endl;
  std::cout << padfMinBound[2] << std::endl;
  std::cout << padfMinBound[3] << std::endl;
  std::cout << padfMaxBound[0] << std::endl;
  std::cout << padfMaxBound[1] << std::endl;
  std::cout << padfMaxBound[2] << std::endl;
  std::cout << padfMaxBound[3] << std::endl;

  for (int i = 0; i < pnEntities; ++i) {
    auto shpobj = SHPReadObject(shphandle, i);

    std::cout << std::endl;
    std::cout << "shape type         " << shpobj->nSHPType << std::endl;
    std::cout << "shape id           " << shpobj->nShapeId << std::endl;
    std::cout << "number of parts    " << shpobj->nParts << std::endl;
    std::cout << "number of vertices " << shpobj->nVertices << std::endl;
    std::cout << "lower corner       " << shpobj->dfXMin << ", "
              << shpobj->dfYMin << std::endl;
    std::cout << "upper corner       " << shpobj->dfXMax << ", "
              << shpobj->dfYMax << std::endl;
    double dfXMin;
    double dfYMin;
    double dfXMax;
    double dfYMax;
  }

  return 1;
}