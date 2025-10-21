#ifndef TEST_OCTREE_H
#define TEST_OCTREE_H

#include "OctreeNode.h"
#include "octree_builder.h"
#include "octree_builder_xy.h"
#include "octree_util.h"
#include "octree_vtk_writer.h"

#include <filesystem>

namespace fs = std::filesystem;

bool test_octree(const fs::path &filepath)
{

  // ============================================================================
  // Demo
  // ============================================================================
  // Root domain = [0,1]^3
  OctreeNode root(0, 0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

  OctreeBuilder builder(5); // max level = 5
  builder.build_tree(&root);

  print_tree(&root);

  // Export to VTK and VTU:
  auto vtk_filename = filepath / "octree.vtk";
  export_to_vtk(&root, vtk_filename);

  auto vtu_filename = filepath / "octree.vtu";
  export_to_vtu(&root, vtu_filename);

  // ============================================================================
  // Octree in XY plane demo
  // ============================================================================
  OctreeBuilderXY octree_builder_xy(5); // max level = 5
  OctreeNode      root_xy(0, 0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

  octree_builder_xy.build_tree(&root_xy);
  print_tree(&root_xy);

  // VTU export:
  auto vtu_filename_xy = filepath / "quadtree_xy.vtu";
  export_to_vtu(&root_xy, vtu_filename_xy);

  return true;
}

#endif // TEST_OCTREE_H