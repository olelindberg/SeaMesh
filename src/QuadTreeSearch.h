#ifndef QUAD_TREE_UTIL_H
#define QUAD_TREE_UTIL_H

#include "QuadTree.h"

#include <memory>
#include <vector>

#include <limits>

class QuadTreeSearchEvent {
public:
  virtual void south_neighbour(const std::shared_ptr<QuadTree> &q,
                               const std::shared_ptr<QuadTree> &p){};
  virtual void east_neighbour(const std::shared_ptr<QuadTree> &q,
                              const std::shared_ptr<QuadTree> &p){};
  virtual void north_neighbour(const std::shared_ptr<QuadTree> &q,
                               const std::shared_ptr<QuadTree> &p){};
  virtual void west_neighbour(const std::shared_ptr<QuadTree> &q,
                              const std::shared_ptr<QuadTree> &p){};
  virtual void leaf(const std::shared_ptr<QuadTree> &q){};
};

class QuadTreeNeighbourSearch : public QuadTreeSearchEvent {
public:
  QuadTreeNeighbourSearch() { _edge_neighbours.resize(4); }

  virtual void south_neighbour(const std::shared_ptr<QuadTree> &q,
                               const std::shared_ptr<QuadTree> &p) override {
    _neighbours.push_back(p);
    _edge_neighbours[0].push_back(p);
  };

  virtual void east_neighbour(const std::shared_ptr<QuadTree> &q,
                              const std::shared_ptr<QuadTree> &p) override {
    _neighbours.push_back(p);
    _edge_neighbours[1].push_back(p);
  };

  virtual void north_neighbour(const std::shared_ptr<QuadTree> &q,
                               const std::shared_ptr<QuadTree> &p) override {
    _neighbours.push_back(p);
    _edge_neighbours[2].push_back(p);
  };

  virtual void west_neighbour(const std::shared_ptr<QuadTree> &q,
                              const std::shared_ptr<QuadTree> &p) override {
    _neighbours.push_back(p);
    _edge_neighbours[3].push_back(p);
  };

  const std::vector<std::shared_ptr<QuadTree>> &get_neighbours() {
    return _neighbours;
  }
  const std::vector<std::vector<std::shared_ptr<QuadTree>>> &
  get_edge_neighbours() {
    return _edge_neighbours;
  }

  void clear() {
    _neighbours.clear();
    for (auto &edge_neighbours : _edge_neighbours)
      edge_neighbours.clear();
  }

private:
  std::vector<std::shared_ptr<QuadTree>> _neighbours;
  std::vector<std::vector<std::shared_ptr<QuadTree>>> _edge_neighbours;
};

class QuadTreeLeafSearch : public QuadTreeSearchEvent {
public:
  virtual void leaf(const std::shared_ptr<QuadTree> &q) { _leafs.push_back(q); }
  std::vector<std::shared_ptr<QuadTree>> &get_leafs() { return _leafs; }

  void clear() { _leafs.clear(); }

private:
  std::vector<std::shared_ptr<QuadTree>> _leafs;
};

class QuadTreeIndexAssigner : public QuadTreeSearchEvent {
public:
  virtual void leaf(const std::shared_ptr<QuadTree> &q) {
    q->set_id(_id);
    ++_id;
  };

  int get_id() { return _id; }

private:
  int _id = 0;
};

class QuadTreeSearch {
public:
  QuadTreeSearch(std::shared_ptr<QuadTreeSearchEvent> search_event =
                     std::shared_ptr<QuadTreeSearchEvent>())
      : _search_event(search_event) {}

  void set_search_event(std::shared_ptr<QuadTreeSearchEvent> search_event) {
    _search_event = search_event;
  }

  void find_neighbours(const std::shared_ptr<QuadTree> &q);
  void find_leafs(const std::shared_ptr<QuadTree> &q);
  void find_leaf_neighbours(const std::shared_ptr<QuadTree> &q);

private:
  double tol = std::numeric_limits<double>::epsilon();

  bool compare_horizontal(const std::shared_ptr<QuadTree> &q,
                          const std::shared_ptr<QuadTree> &p);
  bool compare_vertical(const std::shared_ptr<QuadTree> &q,
                        const std::shared_ptr<QuadTree> &p);
  void search_south(const std::shared_ptr<QuadTree> &q,
                    const std::shared_ptr<QuadTree> &p);
  void search_east(const std::shared_ptr<QuadTree> &q,
                   const std::shared_ptr<QuadTree> &p);
  void search_north(const std::shared_ptr<QuadTree> &q,
                    const std::shared_ptr<QuadTree> &p);
  void search_west(const std::shared_ptr<QuadTree> &q,
                   const std::shared_ptr<QuadTree> &p);
  void find_neighbours_south(const std::shared_ptr<QuadTree> &q,
                             const std::shared_ptr<QuadTree> &p);
  void find_neighbours_east(const std::shared_ptr<QuadTree> &q,
                            const std::shared_ptr<QuadTree> &p);
  void find_neighbours_north(const std::shared_ptr<QuadTree> &q,
                             const std::shared_ptr<QuadTree> &p);
  void find_neighbours_west(const std::shared_ptr<QuadTree> &q,
                            const std::shared_ptr<QuadTree> &p);

  std::shared_ptr<QuadTreeSearchEvent> _search_event;
};

#endif // QUAD_TREE_UTIL_H
