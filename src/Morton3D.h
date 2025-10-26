#ifndef MORTON3D_H
#define MORTON3D_H

#include <cstdint>
#include <iostream>

/// @brief Morton3D provides static methods for encoding and decoding 3D Morton codes (Z-order curve).
/// @details Morton codes interleave the bits of the x, y, and z coordinates to produce a single linear index.
/// This is useful for spatial indexing and hierarchical data structures like octrees.
/// The implementation supports up to 21 bits per coordinate, allowing for grids of size up to 2^21 in each dimension.
/// Usage example:
//// ```cpp
/// uint32_t x = 5, y = 10, z = 3;
/// uint64_t morton_code = Morton3D::linear_index(x, y, z);
/// uint32_t ix, iy, iz;
/// Morton3D::grid_index(morton_code, ix, iy, iz);
/// ```
/// This will encode the (x,y,z) coordinates into a Morton code and then decode it back to (ix,iy,iz).
///

///-----------------------------------------------------------------------------
/// Note on bitwise operations:
///
/// Bitwise operations used are: AND, OR, XOR, shifts.
///
/// AND operator: &
/// logical AND operation on each pair of the corresponding bits.
/// Thus, if both bits in the compared position are 1, the bit in the resulting binary
/// representation is 1 (1 × 1 = 1); otherwise, the result is 0 (1 × 0 = 0 and 0 × 0 = 0)
///
/// OR operator: |
/// logical inclusive OR operation on each pair of corresponding bits.
/// The result in each position is 0 if both bits are 0, while otherwise the result is 1.
///
/// XOR operator: ^
/// logical exclusive OR operation on each pair of corresponding bits.
/// The result in each position is 1 if only one of the bits is 1, but will be 0 if both are 0 or both are 1.
///
/// left shift: <<
/// Moves bits to the left by the specified number of positions.
/// Each left shift effectively multiplies the number by 2.
///
/// right shift: >>
/// Moves bits to the right by the specified number of positions.
/// Each right shift effectively divides the number by 2, discarding any fractional part.
///
//// -----------------------------------------------------------------------------

class Morton3D
{

public:
  // -----------------------------------------------------------------------------
  // Morton encoding utilities (bit interleaving for 3D indices)
  // -----------------------------------------------------------------------------

  static inline uint64_t linear_index(uint32_t x, uint32_t y, uint32_t z) { return _part1by2(x) | (_part1by2(y) << 1) | (_part1by2(z) << 2); }

  static inline void grid_index(uint64_t code, uint32_t &i, uint32_t &j, uint32_t &k)
  {
    i = _compact1by2(code);
    j = _compact1by2(code >> 1);
    k = _compact1by2(code >> 2);
  }

private:
  static inline uint64_t _part1by2(uint64_t n)
  {
    n &= 0x1fffff; // 21 bits
    n = (n | (n << 32)) & 0x1f00000000ffff;
    n = (n | (n << 16)) & 0x1f0000ff0000ff;
    n = (n | (n << 8)) & 0x100f00f00f00f00f;
    n = (n | (n << 4)) & 0x10c30c30c30c30c3;
    n = (n | (n << 2)) & 0x1249249249249249;
    return n;
  }
  static inline uint32_t _compact1by2(uint64_t n)
  {
    n &= 0x1249249249249249;
    n = (n ^ (n >> 2)) & 0x10c30c30c30c30c3;
    n = (n ^ (n >> 4)) & 0x100f00f00f00f00f;
    n = (n ^ (n >> 8)) & 0x1f0000ff0000ff;
    n = (n ^ (n >> 16)) & 0x1f00000000ffff;
    n = (n ^ (n >> 32)) & 0x1fffff;
    return (uint32_t)n;
  }
};

#endif // MORTON3D_H