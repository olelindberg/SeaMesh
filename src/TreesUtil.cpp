#include "TreesUtil.h"

bool TreesUtil::createSurfaceMesh(const OcTree::ptr &octree, QuadTree::ptr &quadtree) {
  quadtree = std::make_shared<QuadTree>(octree->lower(), octree->upper(), 0, 0);
  return _createSurfaceMesh(octree, quadtree);
}

bool TreesUtil::createBottomMesh(const std::shared_ptr<OcTree> &octree, QuadTree::ptr &quadtree) {
  quadtree = std::make_shared<QuadTree>(octree->lower(), octree->upper(), 0, 0);
  return _createBottomMesh(octree, quadtree);
}

bool TreesUtil::_createSurfaceMesh(const OcTree::ptr &octree, QuadTree::ptr &quadtree) {
  if (octree->childs().empty()) {
    return false;
  } else {
    auto center = 0.5 * (octree->lower() + octree->upper());
    auto sw     = std::make_shared<QuadTree>(octree->lower(), center, octree->level() + 1, 0);
    auto se =
        std::make_shared<QuadTree>(Eigen::Vector3d(center(0), octree->lower()(1), octree->lower()(2)), Eigen::Vector3d(octree->upper()(0), center(1), octree->upper()(2)), octree->level() + 1, 1);
    auto nw =
        std::make_shared<QuadTree>(Eigen::Vector3d(octree->lower()(0), center(1), octree->lower()(2)), Eigen::Vector3d(center(0), octree->upper()(1), octree->upper()(2)), octree->level() + 1, 2);
    auto ne = std::make_shared<QuadTree>(center, octree->upper(), octree->level() + 1, 3);

    quadtree->trees.push_back(sw);
    quadtree->trees.push_back(se);
    quadtree->trees.push_back(nw);
    quadtree->trees.push_back(ne);

    _createSurfaceMesh(octree->childs()[OcTree::CHILDREN::SOUTH_WEST_TOP], sw);
    _createSurfaceMesh(octree->childs()[OcTree::CHILDREN::SOUTH_EAST_TOP], se);
    _createSurfaceMesh(octree->childs()[OcTree::CHILDREN::NORTH_WEST_TOP], nw);
    _createSurfaceMesh(octree->childs()[OcTree::CHILDREN::NORTH_EAST_TOP], ne);
  }
  return true;
}

bool TreesUtil::_createBottomMesh(const std::shared_ptr<OcTree> &octree, QuadTree::ptr &quadtree) {
  if (octree->childs().empty()) {
    return false;
  } else {
    auto center = 0.5 * (octree->lower() + octree->upper());
    auto sw     = std::make_shared<QuadTree>(octree->lower(), center, octree->level() + 1, 0);
    auto se =
        std::make_shared<QuadTree>(Eigen::Vector3d(center(0), octree->lower()(1), octree->lower()(2)), Eigen::Vector3d(octree->upper()(0), center(1), octree->upper()(2)), octree->level() + 1, 1);
    auto nw =
        std::make_shared<QuadTree>(Eigen::Vector3d(octree->lower()(0), center(1), octree->lower()(2)), Eigen::Vector3d(center(0), octree->upper()(1), octree->upper()(2)), octree->level() + 1, 2);
    auto ne = std::make_shared<QuadTree>(center, octree->upper(), octree->level() + 1, 3);

    quadtree->trees.push_back(sw);
    quadtree->trees.push_back(se);
    quadtree->trees.push_back(nw);
    quadtree->trees.push_back(ne);

    _createBottomMesh(octree->childs()[OcTree::CHILDREN::SOUTH_WEST_BOTTOM], sw);
    _createBottomMesh(octree->childs()[OcTree::CHILDREN::SOUTH_EAST_BOTTOM], se);
    _createBottomMesh(octree->childs()[OcTree::CHILDREN::NORTH_WEST_BOTTOM], nw);
    _createBottomMesh(octree->childs()[OcTree::CHILDREN::NORTH_EAST_BOTTOM], ne);
  }
  return true;
}