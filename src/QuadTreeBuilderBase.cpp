#include "QuadTreeBuilderBase.h"
#include "QuadTreeSearch.h"

#include <Eigen/Dense>

#include <iostream>

void QuadTreeBuilderBase::subdivide(std::shared_ptr<QuadTree> &tree) {

  auto level = tree->level;

  Eigen::Vector3d center;
  tree->center(center);

  tree->trees.resize(4);

  // Sub tree 0:
  tree->trees[0] = std::make_shared<QuadTree>(tree->lower, center, level + 1, 0);
  tree->trees[0]->set_parent(tree);

  // Sub tree 1:
  tree->trees[1] = std::make_shared<QuadTree>(Eigen::Vector3d(center(0), tree->lower(1), 0.0), Eigen::Vector3d(tree->upper(0), center(1), 0.0), level + 1, 1);
  tree->trees[1]->set_parent(tree);

  // Sub tree 2:
  tree->trees[2] = std::make_shared<QuadTree>(Eigen::Vector3d(tree->lower(0), center(1), 0.0), Eigen::Vector3d(center(0), tree->upper(1), 0.0), level + 1, 2);
  tree->trees[2]->set_parent(tree);

  // Sub tree 3:
  tree->trees[3] = std::make_shared<QuadTree>(center, tree->upper, level + 1, 3);
  tree->trees[3]->set_parent(tree);
}

void QuadTreeBuilderBase::balance(std::shared_ptr<QuadTree> &tree) {

  // Insert leaf nodes/quads into list:
  auto           leaf_search = std::make_shared<QuadTreeLeafSearch>();
  QuadTreeSearch search(leaf_search);
  search.find_leafs(tree);
  auto &quad_list = leaf_search->get_leafs();

  while (quad_list.size() > 0) {
    // Get last quad in list:
    std::shared_ptr<QuadTree> q = quad_list.back();
    quad_list.pop_back();

    // Find neighbours:
    std::shared_ptr<QuadTreeNeighbourSearch> neighbour_search = std::make_shared<QuadTreeNeighbourSearch>();
    QuadTreeSearch                           search(neighbour_search);
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

      q->trees[0] = std::make_shared<QuadTree>(q->lower, center, q->get_depth() + 1, 0);
      q->trees[0]->set_parent(q);

      q->trees[1] = std::make_shared<QuadTree>(Eigen::Vector3d(center(0), q->lower(1), 0.0), Eigen::Vector3d(q->upper(0), center(1), 0.0), q->get_depth() + 1, 1);
      q->trees[1]->set_parent(q);

      q->trees[2] = std::make_shared<QuadTree>(Eigen::Vector3d(q->lower(0), center(1), 0.0), Eigen::Vector3d(center(0), q->upper(1), 0.0), q->get_depth() + 1, 2);
      q->trees[2]->set_parent(q);

      q->trees[3] = std::make_shared<QuadTree>(center, q->upper, q->get_depth() + 1, 3);
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
