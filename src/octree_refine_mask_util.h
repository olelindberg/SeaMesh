#ifndef OCTREE_REFINE_MASK_UTIL_H
#define OCTREE_REFINE_MASK_UTIL_H

#include <cstdint>
#include <string>

enum class RefineMask : uint8_t
{
  REFINE_NONE = 0,
  REFINE_X    = 1 << 0,
  REFINE_Y    = 1 << 1,
  REFINE_Z    = 1 << 2
};

inline RefineMask operator|(RefineMask a, RefineMask b) { return static_cast<RefineMask>(static_cast<uint8_t>(a) | static_cast<uint8_t>(b)); }

inline RefineMask operator&(RefineMask a, RefineMask b) { return static_cast<RefineMask>(static_cast<uint8_t>(a) & static_cast<uint8_t>(b)); }

inline RefineMask &operator|=(RefineMask &a, RefineMask b)
{
  a = a | b;
  return a;
}

inline bool has_flag(RefineMask value, RefineMask flag) { return static_cast<uint8_t>(value & flag) != 0; }

class OctreeRefineMaskUtil
{
public:
  static inline RefineMask mask_or(RefineMask a, RefineMask b) { return static_cast<RefineMask>(static_cast<int>(a) | static_cast<int>(b)); }

  static inline std::string mask_to_string(RefineMask m)
  {
    std::string s;
    if (has_flag(m, RefineMask::REFINE_X))
      s += "X";
    if (has_flag(m, RefineMask::REFINE_Y))
      s += "Y";
    if (has_flag(m, RefineMask::REFINE_Z))
      s += "Z";
    if (s.empty())
      s = "NONE";
    return s;
  }
};

#endif // OCTREE_REFINE_MASK_UTIL_H