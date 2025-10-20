#include <cstdint>
#include <iostream>

class Morton3D {

public:
  // -----------------------------------------------------------------------------
  // Morton encoding utilities (bit interleaving for 3D indices)
  // -----------------------------------------------------------------------------

  static inline uint64_t linear_index(uint32_t x, uint32_t y, uint32_t z) { return _part1by2(x) | (_part1by2(y) << 1) | (_part1by2(z) << 2); }

  static inline void grid_index(uint64_t code, uint32_t &i, uint32_t &j, uint32_t &k) {
    i = _compact1by2(code);
    j = _compact1by2(code >> 1);
    k = _compact1by2(code >> 2);
  }

private:
  static inline uint64_t _part1by2(uint64_t n) {
    n &= 0x1fffff; // 21 bits
    n = (n | (n << 32)) & 0x1f00000000ffff;
    n = (n | (n << 16)) & 0x1f0000ff0000ff;
    n = (n | (n << 8)) & 0x100f00f00f00f00f;
    n = (n | (n << 4)) & 0x10c30c30c30c30c3;
    n = (n | (n << 2)) & 0x1249249249249249;
    return n;
  }
  static inline uint32_t _compact1by2(uint64_t n) {
    n &= 0x1249249249249249;
    n = (n ^ (n >> 2)) & 0x10c30c30c30c30c3;
    n = (n ^ (n >> 4)) & 0x100f00f00f00f00f;
    n = (n ^ (n >> 8)) & 0x1f0000ff0000ff;
    n = (n ^ (n >> 16)) & 0x1f00000000ffff;
    n = (n ^ (n >> 32)) & 0x1fffff;
    return (uint32_t)n;
  }
};