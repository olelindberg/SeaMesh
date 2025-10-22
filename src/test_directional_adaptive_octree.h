#ifndef TEST_DIRECTIONAL_ADAPTIVE_OCTREE_H
#define TEST_DIRECTIONAL_ADAPTIVE_OCTREE_H

// directional_octree_balance.cpp
// Directional (anisotropic) octree with Morton codes and balancing.
// Do NOT use "using namespace std".

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// ------------------------- Morton 3D utilities (3 bits per level) -------------------------
// Interleave helpers (64-bit)
static inline uint64_t part1by2(uint64_t n)
{
  n &= 0x1fffffULL; // 21 bits
  n = (n | (n << 32)) & 0x1f00000000ffffULL;
  n = (n | (n << 16)) & 0x1f0000ff0000ffULL;
  n = (n | (n << 8)) & 0x100f00f00f00f00fULL;
  n = (n | (n << 4)) & 0x10c30c30c30c30c3ULL;
  n = (n | (n << 2)) & 0x1249249249249249ULL;
  return n;
}
static inline uint64_t morton3D_encode(uint32_t x, uint32_t y, uint32_t z) { return part1by2(x) | (part1by2(y) << 1) | (part1by2(z) << 2); }
static inline uint32_t compact1by2(uint64_t n)
{
  n &= 0x1249249249249249ULL;
  n = (n ^ (n >> 2)) & 0x10c30c30c30c30c3ULL;
  n = (n ^ (n >> 4)) & 0x100f00f00f00f00fULL;
  n = (n ^ (n >> 8)) & 0x1f0000ff0000ffULL;
  n = (n ^ (n >> 16)) & 0x1f00000000ffffULL;
  n = (n ^ (n >> 32)) & 0x1fffffULL;
  return (uint32_t)n;
}
static inline void morton3D_decode(uint64_t m, uint32_t &x, uint32_t &y, uint32_t &z)
{
  x = compact1by2(m);
  y = compact1by2(m >> 1);
  z = compact1by2(m >> 2);
}

// Parent / child code helpers for 3-bit-per-level encoding
static inline uint64_t morton_parent(uint64_t code) { return code >> 3; }
static inline uint64_t morton_child(uint64_t parent, int child_id) { return (parent << 3) | (child_id & 0x7); }

// ------------------------- Directional refinement mask -------------------------
enum RefineMask : uint8_t
{
  REFINE_NONE = 0,
  REFINE_X    = 1 << 0,
  REFINE_Y    = 1 << 1,
  REFINE_Z    = 1 << 2
};
static inline RefineMask mask_or(RefineMask a, RefineMask b) { return static_cast<RefineMask>(static_cast<int>(a) | static_cast<int>(b)); }

// ------------------------- Octree node -------------------------
struct OctreeNodeDirAdapt
{
  uint64_t                                         morton;
  int                                              level;
  double                                           xmin, xmax, ymin, ymax, zmin, zmax;
  RefineMask                                       mask; // which axes to refine when subdividing this node
  std::vector<std::unique_ptr<OctreeNodeDirAdapt>> children;

  OctreeNodeDirAdapt(uint64_t m, int l, double x0, double x1, double y0, double y1, double z0, double z1, RefineMask rm = REFINE_NONE)
      : morton(m), level(l), xmin(x0), xmax(x1), ymin(y0), ymax(y1), zmin(z0), zmax(z1), mask(rm)
  {
  }
  bool is_leaf() const { return children.empty(); }
};

// ------------------------- Utility: pretty-print mask -------------------------
static inline std::string mask_to_string(RefineMask m)
{
  std::string s;
  if (m & REFINE_X)
    s += "X";
  if (m & REFINE_Y)
    s += "Y";
  if (m & REFINE_Z)
    s += "Z";
  if (s.empty())
    s = "NONE";
  return s;
}

// ------------------------- Build children according to node->mask -------------------------
auto compute_bounds = [](double minv, double maxv, int n, int idx)
{
  double mid = 0.5 * (minv + maxv);
  if (n == 1)
    return std::make_pair(minv, maxv);
  return idx == 0 ? std::make_pair(minv, mid) : std::make_pair(mid, maxv);
};

static void create_children_for_node(OctreeNodeDirAdapt *node)
{
  // If already has children, do nothing
  if (!node->children.empty())
    return;

  int nx = (node->mask & REFINE_X) ? 2 : 1;
  int ny = (node->mask & REFINE_Y) ? 2 : 1;
  int nz = (node->mask & REFINE_Z) ? 2 : 1;

  double xm = 0.5 * (node->xmin + node->xmax);
  double ym = 0.5 * (node->ymin + node->ymax);
  double zm = 0.5 * (node->zmin + node->zmax);

  // loop over child grid (ix in [0..nx-1], etc.)
  for (int iz = 0; iz < nz; ++iz)
  {
    for (int iy = 0; iy < ny; ++iy)
    {
      for (int ix = 0; ix < nx; ++ix)
      {

        auto [x0, x1] = compute_bounds(node->xmin, node->xmax, nx, ix);
        auto [y0, y1] = compute_bounds(node->ymin, node->ymax, ny, iy);
        auto [z0, z1] = compute_bounds(node->zmin, node->zmax, nz, iz);

        int      child_local_id = (ix) | (iy << 1) | (iz << 2); // 0..7
        uint64_t child_code     = morton_child(node->morton, child_local_id);

        // child initially has no further refinement (mask NONE) - user can set later
        node->children.emplace_back(std::make_unique<OctreeNodeDirAdapt>(child_code, node->level + 1, x0, x1, y0, y1, z0, z1, REFINE_NONE));
      }
    }
  }
}

// ------------------------- Initial directional build (recursive) -------------------------
// choose_mask_fn(centerx, centery, centerz, level) -> RefineMask
static void build_directional_recursive(OctreeNodeDirAdapt *node, int max_level, const std::function<RefineMask(double, double, double, int)> &choose_mask_fn)
{
  if (node->level >= max_level)
    return;

  // decide for this node whether to refine and in which axes
  double cx  = 0.5 * (node->xmin + node->xmax);
  double cy  = 0.5 * (node->ymin + node->ymax);
  double cz  = 0.5 * (node->zmin + node->zmax);
  node->mask = choose_mask_fn(cx, cy, cz, node->level);

  if (node->mask == REFINE_NONE)
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
  int      level;
  bool     operator==(Key const &o) const noexcept { return mort == o.mort && level == o.level; }
};
struct KeyHash
{
  std::size_t operator()(Key const &k) const noexcept { return std::hash<uint64_t>()((k.mort) ^ ((uint64_t)k.level << 48)); }
};

static void collect_leaves(OctreeNodeDirAdapt *node, std::vector<OctreeNodeDirAdapt *> &leaves)
{
  if (node->is_leaf())
  {
    leaves.push_back(node);
    return;
  }
  for (auto &c : node->children)
    collect_leaves(c.get(), leaves);
}

static std::unordered_map<Key, OctreeNodeDirAdapt *, KeyHash> build_leaf_lookup(const std::vector<OctreeNodeDirAdapt *> &leaves)
{
  std::unordered_map<Key, OctreeNodeDirAdapt *, KeyHash> map;
  map.reserve(leaves.size() * 2);
  for (OctreeNodeDirAdapt *n : leaves)
  {
    Key k{n->morton, n->level};
    map[k] = n;
  }
  return map;
}

// ------------------------- Neighbor probing -------------------------------------------------
// Given a node, find a neighbor in direction (dx,dy,dz) at same or coarser/finer levels.
// The function returns the found node pointer and its level, or nullptr if none found (e.g., outside domain).
// We limit search to find the neighbor by walking up (coarser) until present in lookup.
// Assumes domain is [0,1]^3 discretized at integer coords per level.
static OctreeNodeDirAdapt *find_neighbor_in_lookup(const std::unordered_map<Key, OctreeNodeDirAdapt *, KeyHash> &lookup, OctreeNodeDirAdapt *node, int dx, int dy, int dz)
{
  // decode integer coordinates at node.level
  uint32_t ix, iy, iz;
  morton3D_decode(node->morton, ix, iy, iz);
  int L = node->level;
  // indices live in [0, 2^L - 1]
  int grid = (1 << L);
  int nix  = (int)ix + dx;
  int niy  = (int)iy + dy;
  int niz  = (int)iz + dz;
  // domain wrap or clamp? Here we consider boundary as no neighbor (nullptr).
  if (nix < 0 || niy < 0 || niz < 0)
    return nullptr;
  if (nix >= grid || niy >= grid || niz >= grid)
    return nullptr;

  uint32_t nxi = (uint32_t)nix;
  uint32_t nyi = (uint32_t)niy;
  uint32_t nzi = (uint32_t)niz;

  // candidate Morton at same level
  uint64_t candidate = morton3D_encode(nxi, nyi, nzi);
  int      clevel    = L;

  // if exact same-level neighbor present
  Key  kk{candidate, clevel};
  auto it = lookup.find(kk);
  if (it != lookup.end())
    return it->second;

  // otherwise, try coarser neighbors by moving up parents of candidate
  uint64_t pc = candidate;
  int      pl = clevel;
  while (pl > 0)
  {
    pc = morton_parent(pc);
    --pl;
    Key  pk{pc, pl};
    auto it2 = lookup.find(pk);
    if (it2 != lookup.end())
      return it2->second;
  }

  // optionally, could search for finer neighbors (children of candidate at higher level),
  // but balancing wants to refine coarse neighbors, so finding coarser suffices.
  return nullptr;
}

// ------------------------- Balancing algorithm -------------------------------------------
// Ensure that for every face neighbor, |level(node) - level(neighbor)| <= 1.
// If the neighbor is too coarse, mark it for directional refinement (refine only axes necessary).
static void balance_octree(OctreeNodeDirAdapt *root)
{
  bool changed = true;
  int  iter    = 0;
  while (changed)
  {
    ++iter;
    changed = false;

    // collect all leaves and build lookup
    std::vector<OctreeNodeDirAdapt *> leaves;
    collect_leaves(root, leaves);
    auto lookup = build_leaf_lookup(leaves);

    // To avoid modifying tree while iterating, collect nodes to refine and masks to OR
    std::unordered_map<OctreeNodeDirAdapt *, RefineMask> refine_requests;

    // face neighbor directions (6 faces)
    const std::array<std::array<int, 3>, 6> face_dirs = {{{+1, 0, 0}, {-1, 0, 0}, {0, +1, 0}, {0, -1, 0}, {0, 0, +1}, {0, 0, -1}}};

    for (OctreeNodeDirAdapt *n : leaves)
    {
      for (auto d : face_dirs)
      {
        OctreeNodeDirAdapt *nbr = find_neighbor_in_lookup(lookup, n, d[0], d[1], d[2]);
        if (!nbr)
          continue;
        int dl = n->level - nbr->level;
        if (dl > 1)
        {
          // neighbor is too coarse: request refinement of neighbor in directions that affect the shared face
          RefineMask req = REFINE_NONE;
          if (d[0] != 0)
            req = mask_or(req, REFINE_X);
          if (d[1] != 0)
            req = mask_or(req, REFINE_Y);
          if (d[2] != 0)
            req = mask_or(req, REFINE_Z);
          // accumulate requests
          auto it = refine_requests.find(nbr);
          if (it == refine_requests.end())
          {
            refine_requests[nbr] = req;
          }
          else
          {
            it->second = mask_or(it->second, req);
          }
        }
      }
    }

    if (!refine_requests.empty())
    {
      changed = true;
      // perform refinements: for each requested node, OR its mask and create children
      for (auto &pr : refine_requests)
      {
        OctreeNodeDirAdapt *node       = pr.first;
        RefineMask          additional = pr.second;
        RefineMask          before     = node->mask;
        node->mask                     = mask_or(node->mask, additional);
        // create children according to new mask if not already created
        if (node->is_leaf())
        {
          create_children_for_node(node);
        }
        else
        {
          // if node already had children, ensure children exist for the newly requested directions:
          // If node had previously less splitting, we won't re-create; for simplicity we do not handle mixed partial children here.
          // A robust implementation would refine children further or restructure; this demo assumes create_children_for_node will handle current mask.
        }
      }
      // After creating children, continue loop so new leaves/lookup reflect children
    }
  } // end while

  // done balancing
}

// ------------------------- Simple function to print leaves info -------------------------
static void print_leaves(OctreeNodeDirAdapt *root)
{
  std::vector<OctreeNodeDirAdapt *> leaves;
  collect_leaves(root, leaves);
  std::cout << "Leaves: " << leaves.size() << "\n";
  for (OctreeNodeDirAdapt *n : leaves)
  {
    uint32_t ix, iy, iz;
    morton3D_decode(n->morton, ix, iy, iz);
    std::cout << "  level=" << n->level << " morton=" << n->morton << " idx=(" << ix << "," << iy << "," << iz << ")"
              << " bbox=[(" << n->xmin << "," << n->ymin << "," << n->zmin << ")-(" << n->xmax << "," << n->ymax << "," << n->zmax << ")]"
              << " mask=" << mask_to_string(n->mask) << "\n";
  }
}

// ------------------------- Example refinement rule (user replaces with physics) -------------------------
static RefineMask choose_refinement_pattern_example(double cx, double cy, double cz, int level)
{
  // Example: refine in XY near center, refine in XZ for z small, and full refinement near very center

  std::cout << "Choosing refinement at level " << level << " center=(" << cx << "," << cy << "," << cz << ")\n";
  RefineMask m = REFINE_NONE;

  std::cout << "  preliminary mask: " << mask_to_string(m) << "\n";
  std::cout << "  cz check: " << std::fabs(cz - 0.1) << "\n";
  if (std::fabs(cz - 0.1) < 0.2 && level < 4)
  {
    m = mask_or(m, REFINE_X);
    m = mask_or(m, REFINE_Z);
  }

  std::cout << "  preliminary mask: " << mask_to_string(m) << "\n";
  std::cout << "  radial2D check: " << std::sqrt((cx - 0.5) * (cx - 0.5) + (cy - 0.5) * (cy - 0.5)) << "\n";
  //  if (std::sqrt((cx - 0.5) * (cx - 0.5) + (cy - 0.5) * (cy - 0.5)) < 0.2 && level < 3)
  //  {
  //    m = mask_or(m, REFINE_X);
  //    m = mask_or(m, REFINE_Y);
  //  }

  double r2 = (cx - 0.5) * (cx - 0.5) + (cy - 0.5) * (cy - 0.5) + (cz - 0.5) * (cz - 0.5);
  std::cout << "  radial3D check: " << r2 << "\n";
  if (r2 < 0.2 && level < 5)
  {
    m = mask_or(m, REFINE_X);
    m = mask_or(m, REFINE_Y);
    m = mask_or(m, REFINE_Z);
  }

  std::cout << "  chosen mask: " << mask_to_string(m) << "\n";
  return m;
}

void collect_leaves(const OctreeNodeDirAdapt *node, std::vector<const OctreeNodeDirAdapt *> &leaves)
{
  if (node->is_leaf())
  {
    leaves.push_back(node);
  }
  else
  {
    for (const auto &c : node->children)
      collect_leaves(c.get(), leaves);
  }
}

void write_vtu(const OctreeNodeDirAdapt *root, const std::string &filename)
{
  std::vector<const OctreeNodeDirAdapt *> leaves;
  collect_leaves(root, leaves);

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
          ofs << x[i] << " " << y[j] << " " << z[k] << "\n";
  }
  ofs << "        </DataArray>\n";
  ofs << "      </Points>\n";

  // --- Cells ---
  ofs << "      <Cells>\n";
  ofs << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  for (std::size_t c = 0; c < leaves.size(); ++c)
  {
    int base = static_cast<int>(c * 8);
    ofs << base + 0 << " " << base + 1 << " " << base + 3 << " " << base + 2 << " " << base + 4 << " " << base + 5 << " " << base + 7 << " " << base + 6 << "\n";
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
  ofs << "        <DataArray type=\"Int32\" Name=\"Level\" format=\"ascii\">\n";
  for (const auto *n : leaves)
    ofs << n->level << "\n";
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
  std::cout << "✅ Wrote " << filename << " with " << leaves.size() << " cells.\n";
}

// ------------------------- Demo main ---------------------------------------------------------
int test_directional_adaptive_octree()
{
  // Root domain: [0,1]^3, morton=0, level=0
  std::unique_ptr<OctreeNodeDirAdapt> root = std::make_unique<OctreeNodeDirAdapt>(0ULL, 0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, REFINE_NONE);

  // Initial directional build: recursively decide per-node mask and create children
  const int initial_max_level = 4;
  build_directional_recursive(root.get(), initial_max_level, choose_refinement_pattern_example);

  std::cout << "After initial directional build:\n";
  // print_leaves(root.get());
  write_vtu(root.get(), "/home/ole/tmp/directional_octree_before_balance.vtu");
  // Now balance the octree so neighboring leaves differ by at most 1 level
  std::cout << "\nBalancing...\n";
  balance_octree(root.get());

  std::cout << "After balancing:\n";
  // print_leaves(root.get());
  write_vtu(root.get(), "/home/ole/tmp/directional_octree_after_balance.vtu");

  return 0;
}

#endif // TEST_DIRECTIONAL_ADAPTIVE_OCTREE_H