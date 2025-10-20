#include "MortonUtil.h"
//#include "Morton2D.h"
//#include "Morton3D.h"
//
//#include <cstdint>
//#include <iostream>
//#include <string>
//#include <unordered_map>
//#include <vector>

// -----------------------------------------------------------------------------
// Face/Quadrant mapping utilities
// -----------------------------------------------------------------------------
enum Face { XPLUS, XMINUS, YPLUS, YMINUS, ZPLUS, ZMINUS };

int MortonUtil::subface_index_relative(const Cell &fine, const Cell &coarse, Face dir) {
  uint32_t fi, fj, fk;
  Morton3D::grid_index(fine.morton, fi, fj, fk);

  uint32_t ci, cj, ck;
  Morton3D::grid_index(coarse.morton, ci, cj, ck);

  int delta = fine.level - coarse.level;
  if (delta <= 0)
    return -1; // invalid: fine not actually finer

  // Extract the bit identifying which child quadrant (within the coarse cell)
  int      shift = delta - 1;
  uint32_t sub_i = (fi >> shift) & 1;
  uint32_t sub_j = (fj >> shift) & 1;
  uint32_t sub_k = (fk >> shift) & 1;

  // Combine into octant index (x-major order)
  return (sub_k << 2) | (sub_j << 1) | sub_i;

  switch (dir) {
  case XPLUS:
  case XMINUS:
    return (sub_j << 1) | sub_k; // face in Y–Z plane
  case YPLUS:
  case YMINUS:
    return (sub_i << 1) | sub_k; // face in X–Z plane
  case ZPLUS:
  case ZMINUS:
    return (sub_i << 1) | sub_j; // face in X–Y plane
  default:
    return -1;
  }
}

// -----------------------------------------------------------------------------
// Neighbor finder with coarse–fine face mapping
// -----------------------------------------------------------------------------
void MortonUtil::find_neighbors(const Cell &query_cell, const CellMap &cells) {

  const int dirs[6][3] = {
      {+1, 0, 0}, {-1, 0, 0}, // east - west
      {0, +1, 0}, {0, -1, 0}, // north - south
      {0, 0, +1}, {0, 0, -1}  // top - bottom
  };

  const Face faces[6] = {XPLUS, XMINUS, YPLUS, YMINUS, ZPLUS, ZMINUS};

  uint32_t i, j, k;
  Morton3D::grid_index(query_cell.morton, i, j, k);

  std::cout << "Neighbors of (" << i << "," << j << "," << k << ") level " << query_cell.level << ":\n";

  //-------------------------------------------------------------------------//
  // Loop over cell faces::
  //-------------------------------------------------------------------------//
  for (int d = 0; d < 6; ++d) {

    uint32_t ni     = i + dirs[d][0];
    uint32_t nj     = j + dirs[d][1];
    uint32_t nk     = k + dirs[d][2];
    int      nl     = query_cell.level;
    int      fine   = nl + 1;
    int      coarse = nl - 1;

    //-----------------------------------------------------------------------//
    // 1️) Same-level neighbor
    //-----------------------------------------------------------------------//
    Cell n = make_cell(ni, nj, nk, nl);
    if (cells.count(n)) {

      std::cout << "  Neighbor (" << ni << "," << nj << "," << nk << ") level " << nl << "\n";
      continue;
    }

    //-----------------------------------------------------------------------//
    // 2️) Coarser neighbor (parent of neighbor)
    //-----------------------------------------------------------------------//
    if (nl > 0) {

      Cell parent = make_cell(i >> 1, j >> 1, k >> 1, coarse);

      uint32_t pi, pj, pk;
      Morton3D::grid_index(parent.morton, pi, pj, pk);

      Cell parent_neigh = make_cell(pi + dirs[d][0], pj + dirs[d][1], pk + dirs[d][2], coarse);

      if (cells.count(parent_neigh)) {
        int suboct = subface_index_relative(query_cell, parent, faces[d]);
        std::cout << "  Neighbor (" << (pi + dirs[d][0]) << "," << (pj + dirs[d][1]) << "," << (pk + dirs[d][2]) << ") level " << (coarse) << " touches coarse face ";
        std::cout << ((dirs[d][0] > 0) ? "+X" : (dirs[d][0] < 0) ? "-X" : (dirs[d][1] > 0) ? "+Y" : (dirs[d][1] < 0) ? "-Y" : (dirs[d][2] > 0) ? "+Z" : "-Z");
        std::cout << " at suboctant " << suboct << " (" << MortonUtil::octant_name(suboct) << ")\n";
        continue;
      }
    }

    // 3️) Finer neighbors (children)
    uint32_t base_i    = ni << 1;
    uint32_t base_j    = nj << 1;
    uint32_t base_k    = nk << 1;
    bool     foundFine = false;
    for (int di = 0; di < 2; ++di)
      for (int dj = 0; dj < 2; ++dj)
        for (int dk = 0; dk < 2; ++dk) {
          Cell child = make_cell(base_i + di, base_j + dj, base_k + dk, fine);
          if (cells.count(child)) {
            if (!foundFine) {
              std::cout << "  Finer neighbors on face "
                        << ((dirs[d][0] > 0)   ? "+X"
                            : (dirs[d][0] < 0) ? "-X"
                            : (dirs[d][1] > 0) ? "+Y"
                            : (dirs[d][1] < 0) ? "-Y"
                            : (dirs[d][2] > 0) ? "+Z"
                                               : "-Z")
                        << ":\n";
              foundFine = true;
            }
            uint32_t fi, fj, fk;
            Morton3D::grid_index(child.morton, fi, fj, fk);
            int subq = subface_index_relative(child, query_cell, faces[d]);
            std::cout << "    -> (" << fi << "," << fj << "," << fk << ") level " << fine << " touches subquadrant " << subq << " (" << quadrant_name(subq) << ")\n";
          }
        }
  }
}

std::string MortonUtil::quadrant_name(int q) {
  static const char *names[4] = {"lower-left", "lower-right", "upper-left", "upper-right"};
  return names[q];
}

std::string MortonUtil::octant_name(int id) {
  static const char *names[8] = {"west-south-bottom", "east-south-bottom", "west-north-bottom", "east-north-bottom", "west-south-top", "east-south-top", "west-north-top", "east-north-top"};
  return names[id];
}
