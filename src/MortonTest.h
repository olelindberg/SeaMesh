#ifndef MORTON_TEST_H
#define MORTON_TEST_H

#include "Morton2D.h"
#include "MortonUtil.h"
//#include "Morton3D.h"

#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

class MortonTest {

public:
  static int run() {

    int     example = 3;
    CellMap mesh;

    if (example == 0) {
      // Coarse level 1
      mesh[MortonUtil::make_cell(0, 0, 0, 1)] = true;
      mesh[MortonUtil::make_cell(1, 0, 0, 1)] = true;

      // Finer level 2 neighbors of the cell (1,0,0)
      mesh[MortonUtil::make_cell(4, 0, 0, 2)] = true;
      mesh[MortonUtil::make_cell(4, 1, 0, 2)] = true;
      mesh[MortonUtil::make_cell(5, 0, 0, 2)] = true;
      mesh[MortonUtil::make_cell(5, 1, 0, 2)] = true;

      Cell query = MortonUtil::make_cell(1, 0, 0, 1);
      MortonUtil::find_neighbors(query, mesh);
    } else if (example == 1) {

      // 3x3 grid

      mesh[MortonUtil::make_cell(0, 0, 0, 1)] = true;
      mesh[MortonUtil::make_cell(1, 0, 0, 1)] = true;
      mesh[MortonUtil::make_cell(2, 0, 0, 1)] = true;

      mesh[MortonUtil::make_cell(0, 1, 0, 1)] = true;
      mesh[MortonUtil::make_cell(1, 1, 0, 1)] = true;
      mesh[MortonUtil::make_cell(2, 1, 0, 1)] = true;

      mesh[MortonUtil::make_cell(0, 2, 0, 1)] = true;
      mesh[MortonUtil::make_cell(1, 2, 0, 1)] = true;
      mesh[MortonUtil::make_cell(2, 2, 0, 1)] = true;

      // Finer level 2 neighbors of the cell (1,0,0)
      //    mesh[MortonUtil::make_cell(4, 0, 0, 2)] = true;
      //    mesh[MortonUtil::make_cell(4, 1, 0, 2)] = true;
      //    mesh[MortonUtil::make_cell(5, 0, 0, 2)] = true;
      //    mesh[MortonUtil::make_cell(5, 1, 0, 2)] = true;

      MortonUtil::find_neighbors(MortonUtil::make_cell(0, 0, 0, 1), mesh);
      MortonUtil::find_neighbors(MortonUtil::make_cell(1, 0, 0, 1), mesh);
      MortonUtil::find_neighbors(MortonUtil::make_cell(2, 0, 0, 1), mesh);
      MortonUtil::find_neighbors(MortonUtil::make_cell(0, 1, 0, 1), mesh);
      MortonUtil::find_neighbors(MortonUtil::make_cell(1, 1, 0, 1), mesh);
      MortonUtil::find_neighbors(MortonUtil::make_cell(2, 1, 0, 1), mesh);
      MortonUtil::find_neighbors(MortonUtil::make_cell(0, 2, 0, 1), mesh);
      MortonUtil::find_neighbors(MortonUtil::make_cell(1, 2, 0, 1), mesh);
      MortonUtil::find_neighbors(MortonUtil::make_cell(2, 2, 0, 1), mesh);
    } else if (example == 2) {

      std::cout << Morton2D::linear_index(0, 0) << std::endl;
      std::cout << Morton2D::linear_index(1, 0) << std::endl;
      std::cout << Morton2D::linear_index(2, 0) << std::endl;
      std::cout << Morton2D::linear_index(0, 1) << std::endl;
      std::cout << Morton2D::linear_index(1, 1) << std::endl;
      std::cout << Morton2D::linear_index(2, 1) << std::endl;
      std::cout << Morton2D::linear_index(0, 2) << std::endl;
      std::cout << Morton2D::linear_index(1, 2) << std::endl;
      std::cout << Morton2D::linear_index(2, 2) << std::endl;
    } else if (example == 3) {

      // Level 0:
      mesh[MortonUtil::make_cell(0, 0, 0, 0)] = true;

      // Level 1:
      mesh[MortonUtil::make_cell(2, 0, 0, 1)] = true;
      mesh[MortonUtil::make_cell(2, 1, 0, 1)] = true;
      mesh[MortonUtil::make_cell(0, 2, 0, 1)] = true;
      mesh[MortonUtil::make_cell(1, 2, 0, 1)] = true;
      mesh[MortonUtil::make_cell(2, 2, 0, 1)] = true;

      MortonUtil::find_neighbors(MortonUtil::make_cell(2, 0, 0, 1), mesh);
      MortonUtil::find_neighbors(MortonUtil::make_cell(2, 1, 0, 1), mesh);
      MortonUtil::find_neighbors(MortonUtil::make_cell(0, 2, 0, 1), mesh);
      MortonUtil::find_neighbors(MortonUtil::make_cell(1, 2, 0, 1), mesh);
      MortonUtil::find_neighbors(MortonUtil::make_cell(2, 2, 0, 1), mesh);

      MortonUtil::find_neighbors(MortonUtil::make_cell(2, 0, 1, 1), mesh);
      MortonUtil::find_neighbors(MortonUtil::make_cell(2, 1, 1, 1), mesh);
      MortonUtil::find_neighbors(MortonUtil::make_cell(0, 2, 1, 1), mesh);
      MortonUtil::find_neighbors(MortonUtil::make_cell(1, 2, 1, 1), mesh);
      MortonUtil::find_neighbors(MortonUtil::make_cell(2, 2, 1, 1), mesh);
    }

    return 0;
  }
};

#endif // MORTON_TEST_H