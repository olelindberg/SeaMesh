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
struct Cell
{
  uint64_t morton;
  int      level;
  bool     operator==(const Cell &o) const { return morton == o.morton && level == o.level; }
};

struct CellHash
{
  std::size_t operator()(const Cell &c) const noexcept { return std::hash<uint64_t>()(c.morton ^ (c.level * 2654435761ULL)); }
};

using CellMap = std::unordered_map<Cell, bool, CellHash>;

class MortonUtil
{

public:
  // Helper functions
  static inline Cell make_cell(uint32_t i, uint32_t j, uint32_t k, int level)
  {
    uint64_t morton = Morton3D::linear_index(i, j, k);
    return {morton, level};
  }

  // -----------------------------------------------------------------------------
  // Face/Quadrant mapping utilities
  // -----------------------------------------------------------------------------
  enum Face
  {
    XPLUS,
    XMINUS,
    YPLUS,
    YMINUS,
    ZPLUS,
    ZMINUS
  };

  static int  subface_index_relative(const Cell &fine, const Cell &coarse, Face dir);
  static void find_neighbors(const Cell &query_cell, const CellMap &cells);

  // Parent / child code helpers for 3-bit-per-level encoding
  // static inline uint64_t parent(uint64_t code) { return code >> 3; }

  static inline uint64_t parent_x(uint64_t code)
  {
    // Mask pattern for x bits: 0b001001001... (0x1249249249249249)
    constexpr uint64_t xmask = 0x1249249249249249ULL;
    // Extract x bits
    uint64_t xbits = code & xmask;
    // Drop lowest x bit
    xbits >>= 3; // remove one x bit (each is spaced by 3 bits)
    // Reinsert back into the Morton stream
    code &= ~xmask;          // clear x bits
    code |= (xbits & xmask); // put shifted bits back
    return code;
  }

  static inline uint64_t parent_y(uint64_t code)
  {
    constexpr uint64_t ymask = 0x2492492492492492ULL; // pattern shifted by 1
    uint64_t           ybits = code & ymask;
    ybits >>= 3;
    code &= ~ymask;
    code |= (ybits & ymask);
    return code;
  }

  static inline uint64_t parent_z(uint64_t code)
  {
    constexpr uint64_t zmask = 0x4924924924924924ULL; // pattern shifted by 2
    uint64_t           zbits = code & zmask;
    zbits >>= 3;
    code &= ~zmask;
    code |= (zbits & zmask);
    return code;
  }

  static inline uint64_t child_x(uint64_t code, int bit)
  {
    // bit = 0 for "left" child, 1 for "right" child in X
    constexpr uint64_t xmask = 0x1249249249249249ULL;
    uint64_t           xbits = code & xmask;
    xbits                    = (xbits << 3) | (bit ? 0x1 : 0x0);
    code &= ~xmask;
    code |= (xbits & xmask);
    return code;
  }

  static inline uint64_t child_y(uint64_t code, int bit)
  {
    constexpr uint64_t ymask = 0x2492492492492492ULL;
    uint64_t           ybits = code & ymask;
    ybits                    = (ybits << 3) | (bit ? 0x2 : 0x0);
    code &= ~ymask;
    code |= (ybits & ymask);
    return code;
  }

  static inline uint64_t child_z(uint64_t code, int bit)
  {
    constexpr uint64_t zmask = 0x4924924924924924ULL;
    uint64_t           zbits = code & zmask;
    zbits                    = (zbits << 3) | (bit ? 0x4 : 0x0);
    code &= ~zmask;
    code |= (zbits & zmask);
    return code;
  }
  // static inline uint64_t child(uint64_t parent, int child_id) { return (parent << 3) | (child_id & 0x7); }

  static inline uint64_t refine(uint64_t code, int ix, int iy, int iz, bool refine_x, bool refine_y, bool refine_z)
  {
    if (refine_x)
      code = MortonUtil::child_x(code, ix);
    if (refine_y)
      code = MortonUtil::child_y(code, iy);
    if (refine_z)
      code = MortonUtil::child_z(code, iz);
    return code;
  }

private:
  static std::string quadrant_name(int q);
  static std::string octant_name(int id);
};

#endif // MORTON_UTIL_H