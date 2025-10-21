#ifndef OCTREE_VTK_WRITER_H
#define OCTREE_VTK_WRITER_H

#include "OctreeNode.h"
#include "octree_util.h"
#include <fstream>

// ============================================================================
// Export to VTK (POLYDATA cubes)
// ============================================================================
void export_to_vtk(const OctreeNode *root, const std::string &filename)
{
  std::vector<const OctreeNode *> leaves;
  collect_leaves(root, leaves);

  std::ofstream vtk(filename);
  vtk << "# vtk DataFile Version 3.0\n";
  vtk << "Octree Export\n";
  vtk << "ASCII\n";
  vtk << "DATASET POLYDATA\n";

  // Each cube = 8 vertices
  vtk << "POINTS " << leaves.size() * 8 << " float\n";

  for (const auto *n : leaves)
  {
    double x0 = n->xmin, x1 = n->xmax;
    double y0 = n->ymin, y1 = n->ymax;
    double z0 = n->zmin, z1 = n->zmax;

    vtk << x0 << " " << y0 << " " << z0 << "\n";
    vtk << x1 << " " << y0 << " " << z0 << "\n";
    vtk << x1 << " " << y1 << " " << z0 << "\n";
    vtk << x0 << " " << y1 << " " << z0 << "\n";
    vtk << x0 << " " << y0 << " " << z1 << "\n";
    vtk << x1 << " " << y0 << " " << z1 << "\n";
    vtk << x1 << " " << y1 << " " << z1 << "\n";
    vtk << x0 << " " << y1 << " " << z1 << "\n";
  }

  // Connectivity
  int num_cells = leaves.size();
  vtk << "POLYGONS " << num_cells * 6 << " " << num_cells * (6 * 5) << "\n";

  int point_id = 0;
  for (int c = 0; c < num_cells; ++c)
  {
    int p = point_id;
    // Each cube has 6 faces, 4 vertices per face
    int faces[6][4] = {
        {p, p + 1, p + 2, p + 3},     // bottom
        {p + 4, p + 5, p + 6, p + 7}, // top
        {p, p + 1, p + 5, p + 4},     // front
        {p + 1, p + 2, p + 6, p + 5}, // right
        {p + 2, p + 3, p + 7, p + 6}, // back
        {p + 3, p + 0, p + 4, p + 7}  // left
    };
    for (auto &f : faces)
    {
      vtk << "4 " << f[0] << " " << f[1] << " " << f[2] << " " << f[3] << "\n";
    }
    point_id += 8;
  }

  // Optional: Color by refinement level
  vtk << "CELL_DATA " << num_cells * 6 << "\n";
  vtk << "SCALARS level int 1\nLOOKUP_TABLE default\n";
  for (const auto *n : leaves)
    for (int f = 0; f < 6; ++f)
      vtk << n->level << "\n";

  vtk.close();
  std::cout << "✅ Exported " << num_cells << " leaf cells to " << filename << "\n";
}

// ============================================================================
// Export to VTU (UnstructuredGrid with HEXAHEDRON cells)
// ============================================================================
void export_to_vtu(const OctreeNode *root, const std::string &filename)
{
  std::vector<const OctreeNode *> leaves;
  collect_leaves(root, leaves);
  int num_cells  = leaves.size();
  int num_points = num_cells * 8;

  std::ofstream vtu(filename);
  vtu << "<?xml version=\"1.0\"?>\n";
  vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  vtu << "  <UnstructuredGrid>\n";
  vtu << "    <Piece NumberOfPoints=\"" << num_points << "\" NumberOfCells=\"" << num_cells << "\">\n";

  // Points
  vtu << "      <Points>\n";
  vtu << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (const auto *n : leaves)
  {
    double x0 = n->xmin, x1 = n->xmax;
    double y0 = n->ymin, y1 = n->ymax;
    double z0 = n->zmin, z1 = n->zmax;
    vtu << x0 << " " << y0 << " " << z0 << "\n";
    vtu << x1 << " " << y0 << " " << z0 << "\n";
    vtu << x1 << " " << y1 << " " << z0 << "\n";
    vtu << x0 << " " << y1 << " " << z0 << "\n";
    vtu << x0 << " " << y0 << " " << z1 << "\n";
    vtu << x1 << " " << y0 << " " << z1 << "\n";
    vtu << x1 << " " << y1 << " " << z1 << "\n";
    vtu << x0 << " " << y1 << " " << z1 << "\n";
  }
  vtu << "        </DataArray>\n";
  vtu << "      </Points>\n";

  // Cells: connectivity, offsets, types
  vtu << "      <Cells>\n";

  vtu << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  int pid = 0;
  for (int i = 0; i < num_cells; ++i)
  {
    for (int j = 0; j < 8; ++j)
      vtu << pid + j << " ";
    vtu << "\n";
    pid += 8;
  }
  vtu << "        </DataArray>\n";

  vtu << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  for (int i = 1; i <= num_cells; ++i)
    vtu << i * 8 << "\n";
  vtu << "        </DataArray>\n";

  vtu << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  for (int i = 0; i < num_cells; ++i)
    vtu << "12\n"; // VTK_HEXAHEDRON
  vtu << "        </DataArray>\n";

  vtu << "      </Cells>\n";

  // CellData: refinement level
  vtu << "      <CellData Scalars=\"level\">\n";
  vtu << "        <DataArray type=\"Int32\" Name=\"level\" format=\"ascii\">\n";
  for (const auto *n : leaves)
    vtu << n->level << "\n";
  vtu << "        </DataArray>\n";
  vtu << "      </CellData>\n";

  vtu << "    </Piece>\n";
  vtu << "  </UnstructuredGrid>\n";
  vtu << "</VTKFile>\n";

  vtu.close();
  std::cout << "✅ Exported " << num_cells << " hexahedra to " << filename << "\n";
}

#endif // OCTREE_VTK_WRITER_H