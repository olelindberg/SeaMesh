#include "QuadTreePointsInserter.h"
#include "QuadTreeSearch.h"

#include <Eigen/Dense>

#include <iostream>

QuadTreePointsInserter::QuadTreePointsInserter(int levelmax, double lengthmin)
    : _levelmax(levelmax), _lengthmin(lengthmin) {}

void QuadTreePointsInserter::create(std::shared_ptr<QuadTree> &tree,
                                    point_t point) {
  this->create(tree, tree->level, _levelmax, point);
}

void QuadTreePointsInserter::create(std::shared_ptr<QuadTree> &tree, int level,
                                    int levelmax, point_t point) {

  if (level < levelmax && tree->lengthx() > _lengthmin &&
      tree->lengthy() > _lengthmin) {

    Eigen::Vector3d center;
    tree->center(center);

    double x = point.get<0>();
    double y = point.get<1>();

    if (tree->get_xmin() < x && x < tree->get_xmax() && tree->get_ymin() < y &&
        y < tree->get_ymax()) {

      if (tree->trees.empty()) {

        tree->trees.resize(4);

        // Sub tree 0:
        tree->trees[0] =
            std::make_shared<QuadTree>(tree->lower, center, level + 1, 0);
        tree->trees[0]->set_parent(tree);

        // Sub tree 1:
        tree->trees[1] = std::make_shared<QuadTree>(
            Eigen::Vector3d(center(0), tree->lower(1), 0.0),
            Eigen::Vector3d(tree->upper(0), center(1), 0.0), level + 1, 1);
        tree->trees[1]->set_parent(tree);

        // Sub tree 2:
        tree->trees[2] = std::make_shared<QuadTree>(
            Eigen::Vector3d(tree->lower(0), center(1), 0.0),
            Eigen::Vector3d(center(0), tree->upper(1), 0.0), level + 1, 2);
        tree->trees[2]->set_parent(tree);

        // Sub tree 3:
        tree->trees[3] =
            std::make_shared<QuadTree>(center, tree->upper, level + 1, 3);
        tree->trees[3]->set_parent(tree);
      }

      this->create(tree->trees[0], level + 1, levelmax, point);
      this->create(tree->trees[1], level + 1, levelmax, point);
      this->create(tree->trees[2], level + 1, levelmax, point);
      this->create(tree->trees[3], level + 1, levelmax, point);
    }
  }
}

void QuadTreePointsInserter::subdivide(std::shared_ptr<QuadTree> &tree) {

  auto level = tree->level;

  Eigen::Vector3d center;
  tree->center(center);

  tree->trees.resize(4);

  // Sub tree 0:
  tree->trees[0] =
      std::make_shared<QuadTree>(tree->lower, center, level + 1, 0);
  tree->trees[0]->set_parent(tree);

  // Sub tree 1:
  tree->trees[1] = std::make_shared<QuadTree>(
      Eigen::Vector3d(center(0), tree->lower(1), 0.0),
      Eigen::Vector3d(tree->upper(0), center(1), 0.0), level + 1, 1);
  tree->trees[1]->set_parent(tree);

  // Sub tree 2:
  tree->trees[2] = std::make_shared<QuadTree>(
      Eigen::Vector3d(tree->lower(0), center(1), 0.0),
      Eigen::Vector3d(center(0), tree->upper(1), 0.0), level + 1, 2);
  tree->trees[2]->set_parent(tree);

  // Sub tree 3:
  tree->trees[3] =
      std::make_shared<QuadTree>(center, tree->upper, level + 1, 3);
  tree->trees[3]->set_parent(tree);
}

void QuadTreePointsInserter::balance(std::shared_ptr<QuadTree> &tree) {

  // Insert leaf nodes/quads into list:
  auto leaf_search = std::make_shared<QuadTreeLeafSearch>();
  QuadTreeSearch search(leaf_search);
  search.find_leafs(tree);
  auto &quad_list = leaf_search->get_leafs();

  while (quad_list.size() > 0) {
    // Get last quad in list:
    std::shared_ptr<QuadTree> q = quad_list.back();
    quad_list.pop_back();

    // Find neighbours:
    std::shared_ptr<QuadTreeNeighbourSearch> neighbour_search =
        std::make_shared<QuadTreeNeighbourSearch>();
    QuadTreeSearch search(neighbour_search);
    search.find_neighbours(q);
    auto &neighbours = neighbour_search->get_neighbours();

    // Get max depth of neighbours:
    int max_depth = 0;
    for (auto &qi : neighbours)
      max_depth = std::max(max_depth, qi->get_depth());

    // If max depth is to large subdivide this quad and search new quads:
    if (max_depth > q->get_depth() + 1) {
      Eigen::Vector3d center;
      q->center(center);

      q->trees.resize(4);

      q->trees[0] =
          std::make_shared<QuadTree>(q->lower, center, q->get_depth() + 1, 0);
      q->trees[0]->set_parent(q);

      q->trees[1] = std::make_shared<QuadTree>(
          Eigen::Vector3d(center(0), q->lower(1), 0.0),
          Eigen::Vector3d(q->upper(0), center(1), 0.0), q->get_depth() + 1, 1);
      q->trees[1]->set_parent(q);

      q->trees[2] = std::make_shared<QuadTree>(
          Eigen::Vector3d(q->lower(0), center(1), 0.0),
          Eigen::Vector3d(center(0), q->upper(1), 0.0), q->get_depth() + 1, 2);
      q->trees[2]->set_parent(q);

      q->trees[3] =
          std::make_shared<QuadTree>(center, q->upper, q->get_depth() + 1, 3);
      q->trees[3]->set_parent(q);

      // Insert neighbours in list:
      for (auto &qi : neighbours)
        quad_list.push_back(qi);
    }

    // Insert new quads in list:
    for (auto &qi : q->trees)
      quad_list.push_back(qi);
  }
}
