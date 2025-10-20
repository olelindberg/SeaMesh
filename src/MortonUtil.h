#ifndef MORTON_UTIL_H
#define MORTON_UTIL_H

//#include "Morton2D.h"
#include "Morton3D.h"
//
#include <cstdint>
//#include <iostream>
//#include <string>
#include <unordered_map>
//#include <vector>

// -----------------------------------------------------------------------------
// Cell structure
// -----------------------------------------------------------------------------
struct Cell {
  uint64_t morton;
  int      level;
  bool     operator==(const Cell &o) const { return morton == o.morton && level == o.level; }
};

struct CellHash {
  std::size_t operator()(const Cell &c) const noexcept { return std::hash<uint64_t>()(c.morton ^ (c.level * 2654435761ULL)); }
};

using CellMap = std::unordered_map<Cell, bool, CellHash>;

class MortonUtil {

public:
  // Helper functions
  static inline Cell make_cell(uint32_t i, uint32_t j, uint32_t k, int level) {
    uint64_t morton = Morton3D::linear_index(i, j, k);
    return {morton, level};
  }

  // -----------------------------------------------------------------------------
  // Face/Quadrant mapping utilities
  // -----------------------------------------------------------------------------
  enum Face { XPLUS, XMINUS, YPLUS, YMINUS, ZPLUS, ZMINUS };

  static int  subface_index_relative(const Cell &fine, const Cell &coarse, Face dir);
  static void find_neighbors(const Cell &query_cell, const CellMap &cells);

private:
  static std::string quadrant_name(int q);
  static std::string octant_name(int id);
};

#endif // MORTON_UTIL_H