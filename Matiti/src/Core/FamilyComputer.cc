/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#include <Core/FamilyComputer.h> 
#include <Core/Exception.h>

using namespace Matiti;

FamilyComputer::FamilyComputer()
{
}

FamilyComputer::~FamilyComputer()
{
}

// Find which cells the points sit in and create a unordered map that maps nodes to cells
void 
FamilyComputer::createCellNodeMap(const Domain& domain,
                                  const NodePArray& nodeList)
{
  if (!d_map.empty()) {
    d_map.clear();
  }

  for (constNodePIterator iter = nodeList.begin(); iter != nodeList.end(); iter++) {
    NodeP node = *iter;
    IntArray3 cell({{0,0,0}});
    domain.findCellIndex(node->position(), cell);
    long64 cellID = ((long64)cell[0] << 16) | ((long64)cell[1] << 32) | ((long64)cell[2] << 48);
    d_map.insert(CellNodePPair(cellID, node));
  }
}

// Update the cell-node map.
// **WARNING** Will have to be improved later to update only those cells for which particles have
//             entered or left.
void 
FamilyComputer::updateCellNodeMap(const Domain& domain,
                                  const NodePArray& nodeList)
{
  if (!d_map.empty()) {
    d_map.clear();
  }

  for (constNodePIterator iter = nodeList.begin(); iter != nodeList.end(); iter++) {
    NodeP node = *iter;
    const Point3D& position = node->position();
    const Vector3D& displacement = node->displacement();
    Point3D position_new = position+displacement;
    IntArray3 cell({{0,0,0}});
    domain.findCellIndex(position_new, cell);
    long64 cellID = ((long64)cell[0] << 16) | ((long64)cell[1] << 32) | ((long64)cell[2] << 48);
    d_map.insert(CellNodePPair(cellID, node));
  }
}

// print the map
void
FamilyComputer::printCellNodeMap() const
{
  std::cout << "CELL-NODE map: " << std::endl;
  // Print out all the data
  for (auto it = d_map.begin(); it != d_map.end(); ++it) {
    int kk_up = (it->first >> 48);
    int jj_up = (it->first >> 32) & 0xffff;
    int ii_up = (it->first >> 16) & 0xffff;
    std::cout << "    cell key = " << it->first  << " [" << ii_up << "," << jj_up << "," << kk_up <<"] "
              << " node id = " << (it->second)->getID() << std::endl;
  }

  // Check the buckets
  std::cout << "    Bucket_count = " << d_map.bucket_count() << std::endl;
  std::cout << "    cell_node_map's buckets contain:\n";
  for ( unsigned i = 0; i < d_map.bucket_count(); ++i) {
    std::cout << "        bucket #" << i << " contains:";
    for ( auto local_it = d_map.begin(i); local_it!= d_map.end(i); ++local_it ) {
      std::cout << " " << local_it->first << ":" << (local_it->second)->getID();
    }
    std::cout << std::endl;
  }
}

// print the map for one cell
void
FamilyComputer::printCellNodeMap(const IntArray3& cell) const
{
  std::cout << "CELL-NODE map: " << std::endl;
  long64 cellID = ((long64)cell[0] << 16) | ((long64)cell[1] << 32) | ((long64)cell[2] << 48);
  auto nodes = d_map.equal_range(cellID);
  for (auto it = nodes.first; it != nodes.second; ++it) {
    std::cout << "    Cell (" << cell[0] << "," << cell[1] << "," << cell[2]
              <<": key = " << it->first << " value = " << *(it->second) << std::endl;
  }
}

// Finds the family of node m: Find all the nodes inside the horizon of node m
void
FamilyComputer::getInitialFamily(NodeP node,
                                 const Domain& domain,
                                 NodePArray& family) const
{
  // If the node is outside the domain omit it. TODO: Update the family after the node enters the domain.
  if (!domain.inside(node->position())) {
    node->omit(true);
    return;
  }

  // Get horizon size
  double horizon = node->horizonSize();

  // Find cell range within horizon of the node
  Point3D min_pos = node->position() - horizon;
  Point3D max_pos = node->position() + horizon;

  IntArray3 min_cell, max_cell;
  domain.findCellIndex(min_pos, min_cell);
  domain.findCellIndex(max_pos, max_cell);

  IntArray3 num_cells = domain.numCells();
  int iimin = std::max(1, min_cell[0]-1);
  int jjmin = std::max(1, min_cell[1]-1);
  int kkmin = std::max(1, min_cell[2]-1);
  int iimax = std::min(max_cell[0]+1, num_cells[0]);
  int jjmax = std::min(max_cell[1]+1, num_cells[1]);
  int kkmax = std::min(max_cell[2]+1, num_cells[2]);

  IntArray3 cur_cell;
  domain.findCellIndex(node->position(), cur_cell);

  //std::cout << "GET INITIAL FAMILY : " << std::endl;
  //std::cout << " num cells = " << num_cells[0] << "," << num_cells[1] << "," << num_cells[2] << std::endl;
  //std::cout << "  Node = " << node->getID() << " horizon = " << horizon
  //          << " Current cell = [" << cur_cell[0] << "," << cur_cell[1] << "," << cur_cell[2] << "]"
  //          << " Domain cell min = [" << iimin << "," << jjmin << "," << kkmin << "] "
  //          << " max = [" << iimax << "," << jjmax << "," << kkmax << "]" << std::endl;

  //std::cout << "    Near = ";

  // Find the nodes inside the cells within the range
  for (int ii=iimin; ii <= iimax; ++ii) {
    for (int jj=jjmin; jj <= jjmax; ++jj) {
      for (int kk=kkmin; kk <= kkmax; ++kk) {
        long64 cellID = ((long64)ii << 16) | ((long64)jj << 32) | ((long64)kk << 48);
        auto nodes = d_map.equal_range(cellID);
        for (auto it = nodes.first; it != nodes.second; ++it) {
          NodeP near_node = it->second;
          //std::cout << near_node->getID();
          if (node == near_node) {
            //std::cout << ", ";
            continue;
          }
         // ***********for 2D simulation
         // if ((node->distance(*near_node) < horizon) && (node->z() == near_node->z())) {
         // ***********
            if (node->distance(*near_node) < horizon) {
          
            family.push_back(near_node);
            //std::cout << "*, ";
          } 
        } // it loop
      } // kk loop
    } // jj loop
  } // ii loop

  // Throw expecption if the number of nodes in the family is zero
  if (!(family.size() > 0)) {
    std::ostringstream out;
    out << "**ERROR** Cannot have free-hanging nodes.  Family of node " << node->getID() << " is zero " ;
    throw Exception(out.str(), __FILE__, __LINE__);
  }
  //std::cout << std::endl;
}
 
// Finds the family of node m: Find all the nodes inside the horizon of node m
// **WARNING** Need a better way in future (perhaps a private getFamily with position method)
void
FamilyComputer::getCurrentFamily(NodeP node,
                                 const Domain& domain,
                                 NodePArray& family) const
{
  // Find cell range within horizon of the node
  //double horizon = domain.horizon();
  IntArray3 num_cells = domain.numCells();
  const Point3D& position = node->position();
  const Vector3D& displacement = node->displacement();
  Point3D position_new = position+displacement;
  IntArray3 cur_cell;
  domain.findCellIndex(position_new, cur_cell);
  int iimin = std::max(1, cur_cell[0]-1);
  int jjmin = std::max(1, cur_cell[1]-1);
  int kkmin = std::max(1, cur_cell[2]-1);
  int iimax = std::min(cur_cell[0]+1, num_cells[0]);
  int jjmax = std::min(cur_cell[1]+1, num_cells[1]);
  int kkmax = std::min(cur_cell[2]+1, num_cells[2]);

  // Find the nodes inside the cells within the range
  for (int ii=iimin; ii <= iimax; ++ii) {
    for (int jj=jjmin; jj <= jjmax; ++jj) {
      for (int kk=kkmin; kk <= kkmax; ++kk) {
        long64 cellID = ((long64)ii << 16) | ((long64)jj << 32) | ((long64)kk << 48);
        auto nodes = d_map.equal_range(cellID);
        for (auto it = nodes.first; it != nodes.second; ++it) {
          NodeP near_node = it->second;
          if (node == near_node) continue;
          if (node->distance(*near_node) < node->horizonSize()) {
            family.push_back(near_node);
          } 
        } // it loop
      } // kk loop
    } // jj loop
  } // ii loop

  // Throw expecption if the number of nodes in the family is zero
  if (!(family.size() > 0)) {
    std::ostringstream out;
    out << "**ERROR** Cannot have free-hannging nodes.  Family of node " << node->getID() << " is zero " ;
    throw Exception(out.str(), __FILE__, __LINE__);
  }
}

