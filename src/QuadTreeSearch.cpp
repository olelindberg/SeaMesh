#include "QuadTreeSearch.h"

#include <cmath>
#include <memory>
#include <vector>

void QuadTreeSearch::find_neighbours(const std::shared_ptr<QuadTree> &q) {
  find_neighbours_south(q, q);
  find_neighbours_east(q, q);
  find_neighbours_north(q, q);
  find_neighbours_west(q, q);
}

void QuadTreeSearch::find_leafs(const std::shared_ptr<QuadTree> &q) {
  if (q->trees.size() == 0)
    _search_event->leaf(q);
  else
    for (auto &qi : q->trees)
      this->find_leafs(qi);
}

void QuadTreeSearch::find_leaf_neighbours(const std::shared_ptr<QuadTree> &q) {
  if (q->trees.size() == 0)
    this->find_neighbours(q);
  else
    for (auto &qi : q->trees)
      this->find_leaf_neighbours(qi);
}

bool QuadTreeSearch::compare_horizontal(const std::shared_ptr<QuadTree> &q,
                                        const std::shared_ptr<QuadTree> &p) {
  bool is_neighbour = false;

  if (p->get_depth() == q->get_depth()) // Same depth:
  {
    if (std::fabs(q->get_xmin() - p->get_xmin()) < tol) {

      is_neighbour = true;
    }
  } else if (p->get_depth() < q->get_depth()) // This is depper than q:
  {
    // This is smaller than q:
    if ((p->get_xmin() < q->get_xmin() + tol &&
         q->get_xmin() + tol < p->get_xmax()) ||
        (p->get_xmin() < q->get_xmax() - tol &&
         q->get_xmax() - tol < p->get_xmax())) {

      is_neighbour = true;
    }
  } else if (p->get_depth() > q->get_depth()) // q is depper than this:
  {
    if ((q->get_xmin() < p->get_xmin() + tol &&
         p->get_xmin() + tol < q->get_xmax()) ||
        (q->get_xmin() < p->get_xmax() - tol &&
         p->get_xmax() - tol < q->get_xmax())) {

      is_neighbour = true;
    }
  }
  return is_neighbour;
}

bool QuadTreeSearch::compare_vertical(const std::shared_ptr<QuadTree> &q,
                                      const std::shared_ptr<QuadTree> &p) {
  bool is_neighbour = false;

  if (p->get_depth() == q->get_depth()) // Same depth:
  {
    if (std::fabs(q->get_ymin() - p->get_ymin()) < tol) {

      is_neighbour = true;
    }
  } else if (p->get_depth() < q->get_depth()) // This is depper than q:
  {
    if ((p->get_ymin() < q->get_ymin() + tol &&
         q->get_ymin() + tol < p->get_ymax()) ||
        (p->get_ymin() < q->get_ymax() - tol &&
         q->get_ymax() - tol < p->get_ymax())) {

      is_neighbour = true;
    }
  } else if (p->get_depth() > q->get_depth()) // q is depper than this:
  {
    if ((q->get_ymin() < p->get_ymin() + tol &&
         p->get_ymin() + tol < q->get_ymax()) ||
        (q->get_ymin() < p->get_ymax() - tol &&
         p->get_ymax() - tol < q->get_ymax())) {

      is_neighbour = true;
    }
  }
  return is_neighbour;
}

void QuadTreeSearch::search_south(const std::shared_ptr<QuadTree> &q,
                                  const std::shared_ptr<QuadTree> &p) {
  if (p->trees.size() == 0) {
    if (compare_horizontal(q, p))
      _search_event->north_neighbour(q, p);
  } else {
    search_south(q, p->trees[0]);
    search_south(q, p->trees[1]);
  }
}

void QuadTreeSearch::search_east(const std::shared_ptr<QuadTree> &q,
                                 const std::shared_ptr<QuadTree> &p) {
  if (p->trees.size() == 0) {
    if (compare_vertical(q, p))
      _search_event->west_neighbour(q, p);
  } else {
    search_east(q, p->trees[1]);
    search_east(q, p->trees[3]);
  }
}

void QuadTreeSearch::search_north(const std::shared_ptr<QuadTree> &q,
                                  const std::shared_ptr<QuadTree> &p) {
  if (p->trees.size() == 0) {
    if (compare_horizontal(q, p))
      _search_event->south_neighbour(q, p);
  } else {
    search_north(q, p->trees[2]);
    search_north(q, p->trees[3]);
  }
}

void QuadTreeSearch::search_west(const std::shared_ptr<QuadTree> &q,
                                 const std::shared_ptr<QuadTree> &p) {
  if (p->trees.size() == 0) {
    if (compare_vertical(q, p))
      _search_event->east_neighbour(q, p);
  } else {
    search_west(q, p->trees[0]);
    search_west(q, p->trees[2]);
  }
}

void QuadTreeSearch::find_neighbours_south(const std::shared_ptr<QuadTree> &q,
                                           const std::shared_ptr<QuadTree> &p) {
  // If root:
  if (p->get_depth() == 0)
    return;

  if (p->index == 2) // NW-child:
    search_north(q, p->get_parent()->get_child(0));
  else if (p->index == 3) // NE-child:
    search_north(q, p->get_parent()->get_child(1));
  else // SW or SE childs:
    find_neighbours_south(q, p->get_parent());
}

void QuadTreeSearch::find_neighbours_east(const std::shared_ptr<QuadTree> &q,
                                          const std::shared_ptr<QuadTree> &p) {
  // If root:
  if (p->get_depth() == 0)
    return;

  if (p->index == 0) // SW-child:
    search_west(q, p->get_parent()->get_child(1));
  else if (p->index == 2) // NW-child:
    search_west(q, p->get_parent()->get_child(3));
  else // SE or NE childs:
    find_neighbours_east(q, p->get_parent());
}

void QuadTreeSearch::find_neighbours_north(const std::shared_ptr<QuadTree> &q,
                                           const std::shared_ptr<QuadTree> &p) {
  // If root:
  if (p->get_depth() == 0)
    return;

  if (p->index == 0) // SW-child:
    search_south(q, p->get_parent()->get_child(2));
  else if (p->index == 1) // SE-child:
    search_south(q, p->get_parent()->get_child(3));
  else // NW or NE childs:
    find_neighbours_north(q, p->get_parent());
}

void QuadTreeSearch::find_neighbours_west(const std::shared_ptr<QuadTree> &q,
                                          const std::shared_ptr<QuadTree> &p) {
  // If root:
  if (p->get_depth() == 0)
    return;

  if (p->index == 1) // SE-child:
    search_east(q, p->get_parent()->get_child(0));
  else if (p->index == 3) // NE-child:
    search_east(q, p->get_parent()->get_child(2));
  else // SW or NW childs:
    find_neighbours_west(q, p->get_parent());
}
