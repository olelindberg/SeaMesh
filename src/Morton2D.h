#include <cstdint>
#include <iostream>

class Morton2D {

public:
  // Morton linear_index 2D (interleave bits of x and y)
  static inline uint32_t linear_index(uint16_t x, uint16_t y) { return (_part1by1(y) << 1) | _part1by1(x); }

  // Decode Morton number into (x,y)
  static inline void grid_index(uint32_t code, uint16_t &x, uint16_t &y) {
    x = _compact1by1(code);
    y = _compact1by1(code >> 1);
  }

private:
  // Spread bits of a 16-bit integer so that there is one zero bit between each
  static inline uint32_t _part1by1(uint32_t n) {
    n &= 0x0000ffff; // 16 bits
    n = (n | (n << 8)) & 0x00FF00FF;
    n = (n | (n << 4)) & 0x0F0F0F0F;
    n = (n | (n << 2)) & 0x33333333;
    n = (n | (n << 1)) & 0x55555555;
    return n;
  }

  // Extract bits back out (reverse operation)
  static inline uint32_t _compact1by1(uint32_t n) {
    n &= 0x55555555;
    n = (n ^ (n >> 1)) & 0x33333333;
    n = (n ^ (n >> 2)) & 0x0F0F0F0F;
    n = (n ^ (n >> 4)) & 0x00FF00FF;
    n = (n ^ (n >> 8)) & 0x0000FFFF;
    return n;
  }
};