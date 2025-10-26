#ifndef TEST_DIRECTIONAL_ADAPTIVE_OCTREE_H
#define TEST_DIRECTIONAL_ADAPTIVE_OCTREE_H

// Directional (anisotropic) octree with Morton codes and balancing.

#include "directional_adaptive_octree.h"
#include "directional_adaptive_octree_util.h"

int test_directional_adaptive_octree()
{
  // Root domain: [0,1]^3, morton=0, level=0
  std::unique_ptr<DirectionalAdaptiveOctree> root =
      std::make_unique<DirectionalAdaptiveOctree>(0ULL, 0, 0, 0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, RefineMask::REFINE_NONE);

  RefineMaskSelectorExample choose_refinement_pattern_example;

  // Initial directional build: recursively decide per-node mask and create children
  const int initial_max_level = 4;
  DirectionalAdaptiveOctreeUtil::build_directional_recursive(root.get(), initial_max_level, choose_refinement_pattern_example);

  std::cout << "After initial directional build:\n";
  // print_leaves(root.get());
  DirectionalAdaptiveOctreeUtil::write_vtu(root.get(), "/home/ole/tmp/directional_octree_before_balance.vtu");
  // Now balance the octree so neighboring leaves differ by at most 1 level
  std::cout << "\nBalancing...\n";
  DirectionalAdaptiveOctreeUtil::balance_octree(root.get());

  std::cout << "After balancing:\n";
  // print_leaves(root.get());
  DirectionalAdaptiveOctreeUtil::write_vtu(root.get(), "/home/ole/tmp/directional_octree_after_balance.vtu");

  return 0;
}

#endif // TEST_DIRECTIONAL_ADAPTIVE_OCTREE_H