#ifndef DIRECTIONAL_ADAPTIVE_OCTREE_UTIL_H
#define DIRECTIONAL_ADAPTIVE_OCTREE_UTIL_H

#include "MortonUtil.h"
#include "directional_adaptive_octree.h"
#include "logger.h"
#include "octree_util.h"
#include "refine_mask_selector.h"
#include "refine_mask_util.h"

#include <cmath>
#include <deque>
#include <functional>
#include <iomanip>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class RefineMaskSelectorExample : public RefineMaskSelector
{
public:
  virtual RefineMask operator()(const DirectionalAdaptiveOctree *node, double cx, double cy, double cz, int level_x, int level_y, int level_z) const override
  {
    // Example: refine in XY near center, refine in XZ for z small, and full refinement near very center

    int level = std::max({level_x, level_y, level_z});

    std::cout << "Choosing refinement at level " << level << " center=(" << cx << "," << cy << "," << cz << ")\n";
    RefineMask m = RefineMask::REFINE_NONE;

    std::cout << "  preliminary mask: " << RefineMaskUtil::mask_to_string(m) << "\n";
    std::cout << "  cz check: " << std::fabs(cz - 0.1) << "\n";
    if (std::fabs(cz - 0.1) < 0.2 && level < 4)
    {
      m = RefineMaskUtil::mask_or(m, RefineMask::REFINE_X);
      m = RefineMaskUtil::mask_or(m, RefineMask::REFINE_Z);
    }

    std::cout << "  preliminary mask: " << RefineMaskUtil::mask_to_string(m) << "\n";
    std::cout << "  radial2D check: " << std::sqrt((cx - 0.5) * (cx - 0.5) + (cy - 0.5) * (cy - 0.5)) << "\n";
    //  if (std::sqrt((cx - 0.5) * (cx - 0.5) + (cy - 0.5) * (cy - 0.5)) < 0.2 && level < 3)
    //  {
    //    m = RefineMaskUtil::mask_or(m, RefineMask::REFINE_X);
    //    m = RefineMaskUtil::mask_or(m, RefineMask::REFINE_Y);
    //  }

    double r2 = (cx - 0.5) * (cx - 0.5) + (cy - 0.5) * (cy - 0.5) + (cz - 0.5) * (cz - 0.5);
    std::cout << "  radial3D check: " << r2 << "\n";
    if (r2 < 0.2 && level < 5)
    {
      m = RefineMaskUtil::mask_or(m, RefineMask::REFINE_X);
      m = RefineMaskUtil::mask_or(m, RefineMask::REFINE_Y);
      m = RefineMaskUtil::mask_or(m, RefineMask::REFINE_Z);
    }

    std::cout << "  chosen mask: " << RefineMaskUtil::mask_to_string(m) << "\n";
    return m;
  }
};

class DirectionalAdaptiveOctreeUtil
{
public:
  static void create_children_for_node(DirectionalAdaptiveOctree *node)
  {
    // If already has children, do nothing
    if (!node->children.empty())
      return;

    bool refine_x = has_flag(node->mask, RefineMask::REFINE_X);
    bool refine_y = has_flag(node->mask, RefineMask::REFINE_Y);
    bool refine_z = has_flag(node->mask, RefineMask::REFINE_Z);

    int nx = refine_x ? 2 : 1;
    int ny = refine_y ? 2 : 1;
    int nz = refine_z ? 2 : 1;

    struct ChildTemp
    {
      uint64_t                   code;
      double                     xmin, xmax, ymin, ymax, zmin, zmax;
      int                        lx, ly, lz;
      uint64_t                   local_morton;
      DirectionalAdaptiveOctree *parent;
    };

    std::vector<ChildTemp> temp;
    temp.reserve(nx * ny * nz);

    // Create all children
    for (int iz = 0; iz < nz; ++iz)
    {
      for (int iy = 0; iy < ny; ++iy)
      {
        for (int ix = 0; ix < nx; ++ix)
        {
          auto [x0, x1] = compute_bounds(node->xmin, node->xmax, nx, ix);
          auto [y0, y1] = compute_bounds(node->ymin, node->ymax, ny, iy);
          auto [z0, z1] = compute_bounds(node->zmin, node->zmax, nz, iz);

          uint64_t child_code = MortonUtil::refine(node->morton, ix, iy, iz, refine_x, refine_y, refine_z);

          int child_level_x = node->level_x + (refine_x ? 1 : 0);
          int child_level_y = node->level_y + (refine_y ? 1 : 0);
          int child_level_z = node->level_z + (refine_z ? 1 : 0);

          uint64_t local_morton = Morton3D::linear_index(ix, iy, iz);

          temp.push_back({child_code, x0, x1, y0, y1, z0, z1, child_level_x, child_level_y, child_level_z, local_morton, node});
        }
      }
    }

    // Sort by Morton order to ensure consistent spatial ordering
    std::sort(temp.begin(), temp.end(), [](const ChildTemp &a, const ChildTemp &b) { return a.local_morton < b.local_morton; });

    // Move into node->children
    node->children.reserve(temp.size());
    for (auto &c : temp)
    {
      node->children.emplace_back(
          std::make_unique<DirectionalAdaptiveOctree>(c.code, c.lx, c.ly, c.lz, c.xmin, c.xmax, c.ymin, c.ymax, c.zmin, c.zmax, RefineMask::REFINE_NONE, node));
    }
  }

  // ------------------------- Initial directional build (recursive) -------------------------
  // choose_mask_fn(centerx, centery, centerz, level) -> RefineMask
  static void build_directional_recursive(DirectionalAdaptiveOctree *node, int max_level, RefineMaskSelector &choose_mask_fn)
  {
    if (node->level_x >= max_level || node->level_y >= max_level || node->level_z >= max_level)
      return;

    // decide for this node whether to refine and in which axes
    double cx  = 0.5 * (node->xmin + node->xmax);
    double cy  = 0.5 * (node->ymin + node->ymax);
    double cz  = 0.5 * (node->zmin + node->zmax);
    node->mask = choose_mask_fn(node, cx, cy, cz, node->level_x, node->level_y, node->level_z);

    if (node->mask == RefineMask::REFINE_NONE)
      return;

    // create children according to mask and recurse
    create_children_for_node(node);
    for (auto &chptr : node->children)
    {
      build_directional_recursive(chptr.get(), max_level, choose_mask_fn);
    }
  }

  // ------------------------- Collect leaves and build lookup map -------------------------
  struct Key
  {
    uint64_t mort;
    int      level_x;
    int      level_y;
    int      level_z;
    bool     operator==(Key const &o) const noexcept { return mort == o.mort && level_x == o.level_x && level_y == o.level_y && level_z == o.level_z; }
  };

  struct KeyHash
  {
    std::size_t operator()(Key const &k) const noexcept
    {
      uint64_t h   = std::hash<uint64_t>()(k.mort);
      uint64_t mix = ((uint64_t)k.level_x << 16) | ((uint64_t)k.level_y << 8) | (uint64_t)k.level_z;
      // Combine with a 64-bit mix constant
      h ^= mix + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
      return h;
    }
  };

  static void collect_leaves(DirectionalAdaptiveOctree *node, std::vector<DirectionalAdaptiveOctree *> &leaves)
  {
    if (node->is_leaf())
    {
      leaves.push_back(node);
      return;
    }
    for (auto &c : node->children)
      collect_leaves(c.get(), leaves);
  }

  static void print_octrees(const std::vector<DirectionalAdaptiveOctree *> &nodes)
  {
    for (const auto *n : nodes)
    {
      uint32_t ix, iy, iz;
      Morton3D::grid_index(n->morton, ix, iy, iz);
      std::cout << "  level_x=" << n->level_x << " level_y=" << n->level_y << " level_z=" << n->level_z << " morton=" << n->morton << " idx=(" << ix << ","
                << iy << "," << iz << ")"
                << " bbox=[(" << n->xmin << "," << n->ymin << "," << n->zmin << ")-(" << n->xmax << "," << n->ymax << "," << n->zmax << ")]"
                << " mask=" << RefineMaskUtil::mask_to_string(n->mask) << "\n";
    }
  }

  static std::unordered_map<Key, DirectionalAdaptiveOctree *, KeyHash> build_leaf_lookup(const std::vector<DirectionalAdaptiveOctree *> &leaves)
  {
    std::unordered_map<Key, DirectionalAdaptiveOctree *, KeyHash> map;
    map.reserve(leaves.size() * 2);
    for (DirectionalAdaptiveOctree *n : leaves)
    {
      Key k{n->morton, n->level_x, n->level_y, n->level_z};
      map[k] = n;
    }
    return map;
  }

  // ------------------------- Neighbor probing -------------------------------------------------
  // Given a node, find a neighbor in direction (dx,dy,dz) at same or coarser/finer levels.
  // The function returns the found node pointer and its level, or nullptr if none found (e.g., outside domain).
  // We limit search to find the neighbor by walking up (coarser) until present in lookup.
  // Assumes domain is [0,1]^3 discretized at integer coords per level.
  static DirectionalAdaptiveOctree *find_neighbor_in_lookup(const std::unordered_map<Key, DirectionalAdaptiveOctree *, KeyHash> &lookup,
                                                            DirectionalAdaptiveOctree *node, int dx, int dy, int dz)
  {
    // Decode integer coordinates at node's per-axis resolution
    uint32_t ix, iy, iz;
    Morton3D::grid_index(node->morton, ix, iy, iz);

    // Per-axis levels and grid sizes
    int Lx = node->level_x;
    int Ly = node->level_y;
    int Lz = node->level_z;

    const int gx = 1 << Lx;
    const int gy = 1 << Ly;
    const int gz = 1 << Lz;

    // Same-resolution neighbor integer coords
    int nix = (int)ix + dx;
    int niy = (int)iy + dy;
    int niz = (int)iz + dz;

    // If outside domain at this resolution, there is no face neighbor.
    if (nix < 0 || niy < 0 || niz < 0)
      return nullptr;
    if (nix >= gx || niy >= gy || niz >= gz)
      return nullptr;

    // Build same-level Morton candidate
    uint64_t candidate = Morton3D::linear_index((uint32_t)nix, (uint32_t)niy, (uint32_t)niz);

    // 1) Try exact same-level neighbor
    {
      Key  kk{candidate, Lx, Ly, Lz};
      auto it = lookup.find(kk);
      if (it != lookup.end())
        return it->second;
    }

    // 2) Breadth-first search over coarser combinations of axes.
    // We prefer coarsening along the motion axis first, but also allow Y/Z to coarsen.
    struct State
    {
      uint64_t code;
      int      lx, ly, lz;
    };
    std::deque<State> q;
    q.push_back({candidate, Lx, Ly, Lz});

    // (Optional) small visited filter to avoid duplicate checks when coarsening in different orders.
    struct SeenKey
    {
      uint64_t code;
      int      lx, ly, lz;
      bool     operator==(const SeenKey &o) const noexcept { return code == o.code && lx == o.lx && ly == o.ly && lz == o.lz; }
    };
    struct SeenHash
    {
      size_t operator()(const SeenKey &s) const noexcept
      {
        // simple mix
        size_t   h   = std::hash<uint64_t>()(s.code);
        uint64_t mix = ((uint64_t)(uint32_t)s.lx << 16) ^ ((uint64_t)(uint32_t)s.ly << 8) ^ (uint64_t)(uint32_t)s.lz;
        h ^= mix + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
      }
    };
    std::unordered_set<SeenKey, SeenHash> seen;
    seen.insert({candidate, Lx, Ly, Lz});

    auto push_if_new = [&](uint64_t code, int lx, int ly, int lz)
    {
      SeenKey sk{code, lx, ly, lz};
      if (seen.insert(sk).second)
        q.push_back({code, lx, ly, lz});
    };

    while (!q.empty())
    {
      State s = q.front();
      q.pop_front();

      // Try this state in the lookup (skip the initial state which we already checked,
      // but harmless to check again).
      {
        Key  k{s.code, s.lx, s.ly, s.lz};
        auto it = lookup.find(k);
        if (it != lookup.end())
          return it->second;
      }

      // Generate coarser states (prefer coarsening along the motion axis first)
      // Motion-priority order:
      if (dx != 0 && s.lx > 0)
      {
        push_if_new(MortonUtil::parent_x(s.code), s.lx - 1, s.ly, s.lz);
      }
      if (dy != 0 && s.ly > 0)
      {
        push_if_new(MortonUtil::parent_y(s.code), s.lx, s.ly - 1, s.lz);
      }
      if (dz != 0 && s.lz > 0)
      {
        push_if_new(MortonUtil::parent_z(s.code), s.lx, s.ly, s.lz - 1);
      }

      // Also allow coarsening in orthogonal axes, because neighbors can be coarser there too.
      // (These are pushed after motion-axis coarsening to preserve priority.)
      if (dx == 0 && s.lx > 0)
      {
        push_if_new(MortonUtil::parent_x(s.code), s.lx - 1, s.ly, s.lz);
      }
      if (dy == 0 && s.ly > 0)
      {
        push_if_new(MortonUtil::parent_y(s.code), s.lx, s.ly - 1, s.lz);
      }
      if (dz == 0 && s.lz > 0)
      {
        push_if_new(MortonUtil::parent_z(s.code), s.lx, s.ly, s.lz - 1);
      }
    }

    // Not found
    return nullptr;
  }

  // ------------------------- Balancing algorithm -------------------------------------------
  // Ensure that for every face neighbor, |level(node) - level(neighbor)| <= 1.
  // If the neighbor is too coarse, mark it for directional refinement (refine only axes necessary).
  static void balance_octree(DirectionalAdaptiveOctree *root)
  {
    Logger::info("Balancing octree ...");

    bool changed = true;
    int  iter    = 0;

    while (changed)
    {
      ++iter;
      changed = false;

      std::vector<DirectionalAdaptiveOctree *> leaves;
      collect_leaves(root, leaves);
      auto lookup = build_leaf_lookup(leaves);

      std::unordered_map<DirectionalAdaptiveOctree *, RefineMask> refine_requests;

      const std::array<std::array<int, 3>, 6> face_dirs = {{{+1, 0, 0}, {-1, 0, 0}, {0, +1, 0}, {0, -1, 0}, {0, 0, +1}, {0, 0, -1}}};

      for (DirectionalAdaptiveOctree *n : leaves)
      {
        for (auto d : face_dirs)
        {
          const int                  sx = d[0], sy = d[1], sz = d[2];
          DirectionalAdaptiveOctree *nbr = find_neighbor_in_lookup(lookup, n, sx, sy, sz);
          if (!nbr)
            continue;

          // diffs: positive => neighbor is coarser on that axis
          const int dx = n->level_x - nbr->level_x;
          const int dy = n->level_y - nbr->level_y;
          const int dz = n->level_z - nbr->level_z;

          // already balanced?
          if (std::abs(dx) <= 1 && std::abs(dy) <= 1 && std::abs(dz) <= 1)
            continue;

          DirectionalAdaptiveOctree *coarse = nullptr;
          RefineMask                 req    = RefineMask::REFINE_NONE;

          if (dx > 1 || dy > 1 || dz > 1)
          {
            // neighbor is coarser in at least one axis: refine neighbor
            coarse = nbr;
            if (dx > 1)
              req = RefineMaskUtil::mask_or(req, RefineMask::REFINE_X);
            if (dy > 1)
              req = RefineMaskUtil::mask_or(req, RefineMask::REFINE_Y);
            if (dz > 1)
              req = RefineMaskUtil::mask_or(req, RefineMask::REFINE_Z);
          }
          else if (dx < -1 || dy < -1 || dz < -1)
          {
            // current node is coarser in at least one axis: refine current node
            coarse = n;
            if (dx < -1)
              req = RefineMaskUtil::mask_or(req, RefineMask::REFINE_X);
            if (dy < -1)
              req = RefineMaskUtil::mask_or(req, RefineMask::REFINE_Y);
            if (dz < -1)
              req = RefineMaskUtil::mask_or(req, RefineMask::REFINE_Z);
          }

          if (coarse && req != RefineMask::REFINE_NONE)
          {
            auto it = refine_requests.find(coarse);
            if (it == refine_requests.end())
              refine_requests[coarse] = req;
            else
              it->second = RefineMaskUtil::mask_or(it->second, req);
          }
        }
      }

      // Apply refinement requests
      if (!refine_requests.empty())
      {
        changed = true;
        for (auto &[node, addmask] : refine_requests)
        {
          // Apply only new directions
          RefineMask new_mask = RefineMaskUtil::mask_or(node->mask, addmask);
          if (new_mask == node->mask)
            continue; // already refined

          node->mask = new_mask;

          if (node->is_leaf())
            create_children_for_node(node);

          // Upward propagation for consistency (but stop if parent already has it)
          DirectionalAdaptiveOctree *p = node->parent;
          while (p)
          {
            RefineMask before = p->mask;
            p->mask           = RefineMaskUtil::mask_or(p->mask, addmask);
            if (p->mask == before)
              break;
            p = p->parent;
          }

          // Downward refinement only if child's level diff will still violate |Δ|>1
          for (auto &cptr : node->children)
          {
            auto *child = cptr.get();
            if (!child->is_leaf())
              continue;

            bool need_refine = false;
            if (has_flag(node->mask, RefineMask::REFINE_X))
              need_refine |= (child->level_x < node->level_x + 1);
            if (has_flag(node->mask, RefineMask::REFINE_Y))
              need_refine |= (child->level_y < node->level_y + 1);
            if (has_flag(node->mask, RefineMask::REFINE_Z))
              need_refine |= (child->level_z < node->level_z + 1);

            if (need_refine)
              create_children_for_node(child);
          }
        }
      }
    }

    Logger::info("Balancing complete.");
  }

  // ------------------------- Simple function to print leaves info -------------------------
  static void print_leaves(DirectionalAdaptiveOctree *root)
  {
    std::vector<DirectionalAdaptiveOctree *> leaves;
    collect_leaves(root, leaves);
    std::cout << "Leaves: " << leaves.size() << "\n";
    for (DirectionalAdaptiveOctree *n : leaves)
    {
      uint32_t ix, iy, iz;
      Morton3D::grid_index(n->morton, ix, iy, iz);
      std::cout << "  level_x=" << n->level_x << " level_y=" << n->level_y << " level_z=" << n->level_z << " morton=" << n->morton << " idx=(" << ix << ","
                << iy << "," << iz << ")"
                << " bbox=[(" << n->xmin << "," << n->ymin << "," << n->zmin << ")-(" << n->xmax << "," << n->ymax << "," << n->zmax << ")]"
                << " mask=" << RefineMaskUtil::mask_to_string(n->mask) << "\n";
    }
  }

  // ------------------------- Example refinement rule (user replaces with physics) -------------------------
  static RefineMask choose_refinement_pattern_example(double cx, double cy, double cz, int level)
  {
    // Example: refine in XY near center, refine in XZ for z small, and full refinement near very center

    std::cout << "Choosing refinement at level " << level << " center=(" << cx << "," << cy << "," << cz << ")\n";
    RefineMask m = RefineMask::REFINE_NONE;

    std::cout << "  preliminary mask: " << RefineMaskUtil::mask_to_string(m) << "\n";
    std::cout << "  cz check: " << std::fabs(cz - 0.1) << "\n";
    if (std::fabs(cz - 0.1) < 0.2 && level < 4)
    {
      m = RefineMaskUtil::mask_or(m, RefineMask::REFINE_X);
      m = RefineMaskUtil::mask_or(m, RefineMask::REFINE_Z);
    }

    std::cout << "  preliminary mask: " << RefineMaskUtil::mask_to_string(m) << "\n";
    std::cout << "  radial2D check: " << std::sqrt((cx - 0.5) * (cx - 0.5) + (cy - 0.5) * (cy - 0.5)) << "\n";
    //  if (std::sqrt((cx - 0.5) * (cx - 0.5) + (cy - 0.5) * (cy - 0.5)) < 0.2 && level < 3)
    //  {
    //    m = RefineMaskUtil::mask_or(m, RefineMask::REFINE_X);
    //    m = RefineMaskUtil::mask_or(m, RefineMask::REFINE_Y);
    //  }

    double r2 = (cx - 0.5) * (cx - 0.5) + (cy - 0.5) * (cy - 0.5) + (cz - 0.5) * (cz - 0.5);
    std::cout << "  radial3D check: " << r2 << "\n";
    if (r2 < 0.2 && level < 5)
    {
      m = RefineMaskUtil::mask_or(m, RefineMask::REFINE_X);
      m = RefineMaskUtil::mask_or(m, RefineMask::REFINE_Y);
      m = RefineMaskUtil::mask_or(m, RefineMask::REFINE_Z);
    }

    std::cout << "  chosen mask: " << RefineMaskUtil::mask_to_string(m) << "\n";
    return m;
  }

  static void collect_leaves(DirectionalAdaptiveOctree *node, std::vector<const DirectionalAdaptiveOctree *> &leaves)
  {
    if (node->is_leaf())
    {
      leaves.push_back(node);
      return;
    }
    for (const auto &c : node->children)
      collect_leaves(c.get(), leaves);
  }

  static void collect_leaves(const DirectionalAdaptiveOctree *node, std::vector<const DirectionalAdaptiveOctree *> &leaves)
  {
    if (node->is_leaf())
    {
      leaves.push_back(node);
      return;
    }
    for (const auto &c : node->children)
      collect_leaves(c.get(), leaves);
  }

  static void write_vtu(const DirectionalAdaptiveOctree *root, const std::string &filename)
  {
    std::vector<const DirectionalAdaptiveOctree *> leaves;
    DirectionalAdaptiveOctreeUtil::collect_leaves(root, leaves);

    std::ofstream ofs(filename);
    ofs << std::fixed << std::setprecision(6);
    ofs << "<?xml version=\"1.0\"?>\n";
    ofs << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    ofs << "  <UnstructuredGrid>\n";
    ofs << "    <Piece NumberOfPoints=\"" << leaves.size() * 8 << "\" NumberOfCells=\"" << leaves.size() << "\">\n";

    // --- Points ---
    ofs << "      <Points>\n";
    ofs << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    for (const auto *n : leaves)
    {
      const double x[2] = {n->xmin, n->xmax};
      const double y[2] = {n->ymin, n->ymax};
      const double z[2] = {n->zmin, n->zmax};
      for (int k = 0; k < 2; ++k)
        for (int j = 0; j < 2; ++j)
          for (int i = 0; i < 2; ++i)
            ofs << y[j] << " " << x[i] << " " << z[k] << "\n";
    }
    ofs << "        </DataArray>\n";
    ofs << "      </Points>\n";

    // --- Cells ---
    ofs << "      <Cells>\n";
    ofs << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (std::size_t c = 0; c < leaves.size(); ++c)
    {
      int base = static_cast<int>(c * 8);
      ofs << base + 0 << " " << base + 1 << " " << base + 3 << " " << base + 2 << " " << base + 4 << " " << base + 5 << " " << base + 7 << " " << base + 6
          << "\n";
    }
    ofs << "        </DataArray>\n";

    ofs << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (std::size_t c = 1; c <= leaves.size(); ++c)
      ofs << c * 8 << "\n";
    ofs << "        </DataArray>\n";

    ofs << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (std::size_t c = 0; c < leaves.size(); ++c)
      ofs << 12 << "\n"; // VTK_HEXAHEDRON = 12
    ofs << "        </DataArray>\n";
    ofs << "      </Cells>\n";

    // --- Cell Data: Level or Morton code ---
    ofs << "      <CellData Scalars=\"Level\">\n";

    ofs << "        <DataArray type=\"Int32\" Name=\"Levelx\" format=\"ascii\">\n";
    for (const auto *n : leaves)
      ofs << n->level_x << "\n";
    ofs << "        </DataArray>\n";

    ofs << "        <DataArray type=\"Int32\" Name=\"Levely\" format=\"ascii\">\n";
    for (const auto *n : leaves)
      ofs << n->level_y << "\n";
    ofs << "        </DataArray>\n";

    ofs << "        <DataArray type=\"Int32\" Name=\"Levelz\" format=\"ascii\">\n";
    for (const auto *n : leaves)
      ofs << n->level_z << "\n";
    ofs << "        </DataArray>\n";

    ofs << "        <DataArray type=\"UInt64\" Name=\"Morton\" format=\"ascii\">\n";
    for (const auto *n : leaves)
      ofs << n->morton << "\n";
    ofs << "        </DataArray>\n";
    ofs << "      </CellData>\n";

    ofs << "    </Piece>\n";
    ofs << "  </UnstructuredGrid>\n";
    ofs << "</VTKFile>\n";

    ofs.close();

    Logger::info("✅ Wrote " + filename + " with " + std::to_string(leaves.size()) + " cells.");
  }
};

static void set_parent_links(DirectionalAdaptiveOctree *node)
{
  for (auto &ch : node->children)
  {
    ch->parent = node;
    set_parent_links(ch.get());
  }
}

#endif // DIRECTIONAL_ADAPTIVE_OCTREE_UTIL_H