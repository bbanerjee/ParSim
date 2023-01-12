/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#include <Core/Grid/BoundaryConditions/BCDataArray.h>

#include <Core/Geometry/Point.h>
#include <Core/Grid/BoundaryConditions/BoundCondFactory.h>
#include <Core/Grid/BoundaryConditions/DifferenceBCData.h>
#include <Core/Grid/BoundaryConditions/SideBCData.h>
#include <Core/Grid/Box.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Util/DebugStream.h>
#include <algorithm>
#include <functional>
#include <iostream>
#include <set>
#include <vector>

namespace {
// usage: export SCI_DEBUG="BCDA_DBG:+"
Uintah::DebugStream BCDA_dbg("BCDA_DBG", "Grid_BoundaryConditions", "", false);
}

namespace Uintah {

BCDataArray::BCDataArray() {}

BCDataArray::~BCDataArray()
{
  // cout << "Calling BCDataArray destructor" << std::endl;
  bcDataArrayType::const_iterator mat_id_itr;
  for (mat_id_itr = d_BCDataArray.begin(); mat_id_itr != d_BCDataArray.end();
       ++mat_id_itr) {
    int mat_id               = mat_id_itr->first;
    vector<BCGeomBase*>& vec = d_BCDataArray[mat_id];
    vector<BCGeomBase*>::const_iterator bcd_itr;
    for (bcd_itr = vec.begin(); bcd_itr != vec.end(); ++bcd_itr) {
      delete *bcd_itr;
    }
    vec.clear();
  }
  d_BCDataArray.clear();
}

BCDataArray::BCDataArray(const BCDataArray& mybc)
{
  bcDataArrayType::const_iterator mat_id_itr;
  for (mat_id_itr = mybc.d_BCDataArray.begin();
       mat_id_itr != mybc.d_BCDataArray.end();
       ++mat_id_itr) {
    int mat_id                             = mat_id_itr->first;
    const vector<BCGeomBase*>& mybc_vec    = mat_id_itr->second;
    vector<BCGeomBase*>& d_BCDataArray_vec = d_BCDataArray[mat_id];
    vector<BCGeomBase*>::const_iterator vec_itr;
    for (vec_itr = mybc_vec.begin(); vec_itr != mybc_vec.end(); ++vec_itr) {
      d_BCDataArray_vec.push_back((*vec_itr)->clone());
    }
  }
}

BCDataArray&
BCDataArray::operator=(const BCDataArray& rhs)
{
  if (this == &rhs) {
    return *this;
  }

  // Delete the lhs
  bcDataArrayType::const_iterator mat_id_itr;
  for (mat_id_itr = d_BCDataArray.begin(); mat_id_itr != d_BCDataArray.end();
       ++mat_id_itr) {
    int mat_id               = mat_id_itr->first;
    vector<BCGeomBase*>& vec = d_BCDataArray[mat_id];
    vector<BCGeomBase*>::const_iterator bcd_itr;
    for (bcd_itr = vec.begin(); bcd_itr != vec.end(); ++bcd_itr) {
      delete *bcd_itr;
    }
    vec.clear();
  }
  d_BCDataArray.clear();
  // Copy the rhs to the lhs
  for (mat_id_itr = rhs.d_BCDataArray.begin();
       mat_id_itr != rhs.d_BCDataArray.end();
       ++mat_id_itr) {
    int mat_id                             = mat_id_itr->first;
    vector<BCGeomBase*>& d_BCDataArray_vec = d_BCDataArray[mat_id];
    const vector<BCGeomBase*>& rhs_vec     = mat_id_itr->second;
    vector<BCGeomBase*>::const_iterator vec_itr;
    for (vec_itr = rhs_vec.begin(); vec_itr != rhs_vec.end(); ++vec_itr) {
      d_BCDataArray_vec.push_back((*vec_itr)->clone());
    }
  }
  return *this;
}

BCDataArray*
BCDataArray::clone()
{
  return scinew BCDataArray(*this);
}

void
BCDataArray::determineIteratorLimits(Patch::FaceType face, const Patch* patch)
{
  IntVector lpts, hpts;
  patch->getFaceCells(face, -1, lpts, hpts);
  vector<Point> test_pts;

  for (CellIterator candidatePoints(lpts, hpts); !candidatePoints.done();
       candidatePoints++) {
    IntVector nodes[8];
    patch->findNodesFromCell(*candidatePoints, nodes);
    Point pts[8];
    Vector p(0.0, 0.0, 0.0);
    for (int i = 0; i < 8; i++) {
      pts[i] = patch->getLevel()->getNodePosition(nodes[i]);
    }
    if (face == Patch::xminus) {
      p = (pts[0].asVector() + pts[1].asVector() + pts[2].asVector() +
           pts[3].asVector()) /
          4.;
    }
    if (face == Patch::xplus) {
      p = (pts[4].asVector() + pts[5].asVector() + pts[6].asVector() +
           pts[7].asVector()) /
          4.;
    }
    if (face == Patch::yminus) {
      p = (pts[0].asVector() + pts[1].asVector() + pts[4].asVector() +
           pts[5].asVector()) /
          4.;
    }
    if (face == Patch::yplus) {
      p = (pts[2].asVector() + pts[3].asVector() + pts[6].asVector() +
           pts[7].asVector()) /
          4.;
    }
    if (face == Patch::zminus) {
      p = (pts[0].asVector() + pts[2].asVector() + pts[4].asVector() +
           pts[6].asVector()) /
          4.;
    }
    if (face == Patch::zplus) {
      p = (pts[1].asVector() + pts[3].asVector() + pts[5].asVector() +
           pts[7].asVector()) /
          4.;
    }

    test_pts.push_back(Point(p.x(), p.y(), p.z()));
  }

  BCDataArray::bcDataArrayType::iterator mat_id_itr;
  for (mat_id_itr = d_BCDataArray.begin(); mat_id_itr != d_BCDataArray.end();
       ++mat_id_itr) {
    vector<BCGeomBase*>& bc_objects = mat_id_itr->second;
    for (auto& bc_object : bc_objects) {
      bc_object->determineIteratorLimits(face, patch, test_pts);
#if 0
      (*obj)->printLimits();
#endif
    }
  }

  // A BCDataArry contains a bunch of geometry objects. Here, we remove objects
  // with empty iterators. The reason that we get geometry objects with zero
  // iterators is that, for a given boundary face that is shared across several
  // patches, geometry objects are created on ALL patches. Later, a geometry
  // object is assigned an iterator depending on whether it lives on that patch
  // or not.
  for (mat_id_itr = d_BCDataArray.begin(); mat_id_itr != d_BCDataArray.end();
       ++mat_id_itr) {
    std::vector<BCGeomBase*>& bc_objects = mat_id_itr->second;
    for (std::vector<BCGeomBase*>::iterator obj = bc_objects.begin();
         obj < bc_objects.end();) {
      if (!((*obj)->hasIterator())) {
        delete *obj;
        obj = bc_objects.erase(obj); // point the iterator to the next element
                                     // that was after the one we deleted
      } else {
        ++obj;
      }
    }
  }
}

void
BCDataArray::determineInteriorBndIteratorLimits(Patch::FaceType face,
                                                const Patch* patch)
{
  BCDataArray::bcDataArrayType::iterator mat_id_itr;
  for (mat_id_itr = d_BCDataArray.begin(); mat_id_itr != d_BCDataArray.end();
       ++mat_id_itr) {
    vector<BCGeomBase*>& bc_objects = mat_id_itr->second;
    for (vector<BCGeomBase*>::iterator obj = bc_objects.begin();
         obj != bc_objects.end();
         ++obj) {
      (*obj)->determineInteriorBndIteratorLimits(face, patch);
#if 0
      (*obj)->printLimits();
#endif
    }
  }

  // A BCDataArry contains a bunch of geometry objects. Here, we remove objects
  // with empty iterators. The reason that we get geometry objects with zero
  // iterators is that, for a given boundary face that is shared across several
  // patches, geometry objects are created on ALL patches. Later, a geometry
  // object is assigned an iterator depending on whether it lives on that patch
  // or not.
  for (mat_id_itr = d_BCDataArray.begin(); mat_id_itr != d_BCDataArray.end();
       ++mat_id_itr) {
    vector<BCGeomBase*>& bc_objects = mat_id_itr->second;
    for (vector<BCGeomBase*>::iterator obj = bc_objects.begin();
         obj < bc_objects.end();) {
      if (!((*obj)->hasIterator())) {
        delete *obj;
        obj = bc_objects.erase(obj); // point the iterator to the next element
                                     // that was after the one we deleted
      } else {
        ++obj;
      }
    }
  }
}

void
BCDataArray::addBCData(int mat_id, BCGeomBase* bc)
{
  vector<BCGeomBase*>& d_BCDataArray_vec = d_BCDataArray[mat_id];
  d_BCDataArray_vec.push_back(bc);
}

void
BCDataArray::combineBCGeometryTypes(int mat_id)
{

  vector<BCGeomBase*>& d_BCDataArray_vec = d_BCDataArray[mat_id];

  vector<BCGeomBase*> new_bcdata_array;
  // Look to see if there are duplicate SideBCData types, if so, then
  // combine them into one (i.e. copy the BCData from the duplicate into
  // the one that will actually be stored).

  if (count_if(d_BCDataArray_vec.begin(),
               d_BCDataArray_vec.end(),
               cmp_type<SideBCData>()) > 1) {
    BCDA_dbg << "Have duplicates Before . . ." << std::endl;
    for (vector<BCGeomBase*>::const_iterator v_itr = d_BCDataArray_vec.begin();
         v_itr != d_BCDataArray_vec.end();
         ++v_itr) {
      (*v_itr)->print();
    }
    BCDA_dbg << std::endl << std::endl;
  }

  if (count_if(d_BCDataArray_vec.begin(),
               d_BCDataArray_vec.end(),
               cmp_type<SideBCData>()) > 1) {

    auto side_bc = scinew SideBCData();
    for (const auto& bcDataP : d_BCDataArray_vec) {
      if (typeid(*bcDataP) == typeid(SideBCData)) {
        BCDA_dbg << "Found SideBCData" << std::endl;
        BCData bcd, s_bcd;
        bcDataP->getBCData(bcd);
        side_bc->getBCData(s_bcd);
        s_bcd.combine(bcd);
        side_bc->addBCData(s_bcd);
        side_bc->print();
      } else {
        new_bcdata_array.push_back(bcDataP->clone());
      }
    }
    side_bc->print();
    new_bcdata_array.push_back(side_bc->clone());
    delete side_bc;
    new_bcdata_array.back()->print();

    BCDA_dbg << "Have duplicates After . . ." << std::endl;
    for (vector<BCGeomBase*>::const_iterator v_itr = new_bcdata_array.begin();
         v_itr != new_bcdata_array.end();
         ++v_itr) {
      (*v_itr)->print();
    }
    for_each(d_BCDataArray_vec.begin(),
             d_BCDataArray_vec.end(),
             delete_object<BCGeomBase>());
    d_BCDataArray_vec.clear();
    d_BCDataArray_vec = new_bcdata_array;
  }
}

void
BCDataArray::combineBCGeometryTypes_NEW(int mat_id)
{
  if (d_BCDataArray[mat_id].size() <= 1) {
    BCDA_dbg << "One or fewer elements in BCDataArray" << std::endl << std::endl;
    return;
  }

  vector<BCGeomBase*>& d_BCDataArray_vec = d_BCDataArray[mat_id];

  vector<BCGeomBase*> new_bcdata_array;
  // Look to see if there are duplicate SideBCData types, if so, then
  // combine them into one (i.e. copy the BCData from the duplicate into
  // the one that will actually be stored).

  //  count the number of unique geometry types

  vector<BCGeomBase*>::iterator v_itr, nv_itr;

  for (v_itr = d_BCDataArray_vec.begin(); v_itr != d_BCDataArray_vec.end();
       ++v_itr) {
    BCDA_dbg << "number of SideBCData = "
             << count_if(d_BCDataArray_vec.begin(),
                         d_BCDataArray_vec.end(),
                         cmp_type<SideBCData>())
             << std::endl;
  }

  if (count_if(d_BCDataArray_vec.begin(),
               d_BCDataArray_vec.end(),
               cmp_type<SideBCData>()) > 1) {
    BCDA_dbg << "Have duplicates Before . . ." << std::endl;
    for (const auto& bcDataP : d_BCDataArray_vec) {
      BCDA_dbg << "type of element = " << typeid(*bcDataP).name() << std::endl;
      bcDataP->print();
    }
    BCDA_dbg << std::endl << std::endl;
  }

  // Put the last element in the d_BCDataArray_vec into the new_bcdata_array
  // and delete this element

  BCGeomBase* element       = d_BCDataArray_vec.back();
  BCGeomBase* clone_element = element->clone();

  new_bcdata_array.push_back(clone_element);
  delete element;
  d_BCDataArray_vec.pop_back();

  while (!d_BCDataArray_vec.empty()) {
    element      = d_BCDataArray_vec.back();
    bool foundit = false;
    for (nv_itr = new_bcdata_array.begin(); nv_itr != new_bcdata_array.end();
         ++nv_itr) {
      if (*(*nv_itr) == *element) {
        foundit = true;
        break;
      }
    }
    if (foundit) {
      d_BCDataArray_vec.pop_back();
      BCData bcd, n_bcd;
      element->getBCData(bcd);
      (*nv_itr)->getBCData(n_bcd);
      n_bcd.combine(bcd);
      (*nv_itr)->addBCData(n_bcd);
      delete element;
    } else {
      new_bcdata_array.push_back(element->clone());
      d_BCDataArray_vec.pop_back();
      delete element;
    }
  }

  BCDA_dbg << "size of new_bcdata_array = " << new_bcdata_array.size() << std::endl;
  BCDA_dbg << "size of d_BCDataArray_vec = " << d_BCDataArray_vec.size()
           << std::endl;

  for (nv_itr = new_bcdata_array.begin(); nv_itr != new_bcdata_array.end();
       ++nv_itr) {
    (*nv_itr)->print();
    BCDA_dbg << std::endl << std::endl;
  }

  for_each(d_BCDataArray_vec.begin(),
           d_BCDataArray_vec.end(),
           delete_object<BCGeomBase>());

  d_BCDataArray_vec.clear();
#if 1
  d_BCDataArray_vec = new_bcdata_array;
#endif
}

const BoundCondBaseP
BCDataArray::getBoundCondData(int mat_id, const string type, int ichild) const
{
  //  cout << "type = " << type << std::endl;
  BCData new_bc, new_bc_all;
  // Need to check two scenarios -- the given mat_id and the all mat_id (-1)
  // Check the given mat_id
  auto itr = d_BCDataArray.find(mat_id);

  if (itr != d_BCDataArray.end()) {
    itr->second[ichild]->getBCData(new_bc);
    bool found_it = new_bc.find(type);
    if (found_it == true) {
      return new_bc.getBCValues(type);
    } else {
      found_it = new_bc.find("Auxiliary");
      if (found_it) {
        return new_bc.getBCValues("Auxiliary");
      }
    }
  }
  // Check the mat_id = "all" case
  itr = d_BCDataArray.find(-1);
  if (itr != d_BCDataArray.end()) {
    if (ichild < (int)itr->second.size()) {
      itr->second[ichild]->getBCData(new_bc_all);
      bool found_it = new_bc_all.find(type);
      if (found_it == true) {
        return new_bc_all.getBCValues(type);
      } else {
        found_it = new_bc_all.find("Auxiliary");
        if (found_it) {
          return new_bc_all.getBCValues("Auxiliary");
        }
      }
      return nullptr;
    }
  }
  return nullptr;
}

bool
BCDataArray::checkForBoundCondData(int& mat_id, const string& type, int ichild)
{
  BCData new_bc, new_bc_all;
  // Need to check two scenarios -- the given mat_id and the all mat_id (-1)
  // will update mat_id, to represent the applicable material.
  // Check the given mat_id
  bcDataArrayType::const_iterator itr = d_BCDataArray.find(mat_id);
  if (itr != d_BCDataArray.end()) {
    itr->second[ichild]->getBCData(new_bc);
    bool found_it = new_bc.find(type);
    if (found_it == true) {
      return true;
    }
  }
  // Check the mat_id = "all" case
  itr = d_BCDataArray.find(-1);
  if (itr != d_BCDataArray.end()) {
    if (ichild < (int)itr->second.size()) {
      itr->second[ichild]->getBCData(new_bc_all);
      bool found_it = new_bc_all.find(type);
      if (found_it == true) {
        mat_id = -1;
        return true;
      }
    }
  }
  return false;
}

void
BCDataArray::getCellFaceIterator(int mat_id, Iterator& b_ptr, int ichild) const
{
  auto itr = d_BCDataArray.find(mat_id);
  if (itr != d_BCDataArray.end()) {
    itr->second[ichild]->getCellFaceIterator(b_ptr);
  } else {
    itr = d_BCDataArray.find(-1);
    if (itr != d_BCDataArray.end()) {
      itr->second[ichild]->getCellFaceIterator(b_ptr);
    }
  }
}

void
BCDataArray::getNodeFaceIterator(int mat_id, Iterator& b_ptr, int ichild) const
{
  auto itr = d_BCDataArray.find(mat_id);
  if (itr != d_BCDataArray.end()) {
    itr->second[ichild]->getNodeFaceIterator(b_ptr);
  } else {
    itr = d_BCDataArray.find(-1);
    if (itr != d_BCDataArray.end()) {
      itr->second[ichild]->getNodeFaceIterator(b_ptr);
    }
  }
}

int
BCDataArray::getNumberChildren(int mat_id) const
{
  auto itr = d_BCDataArray.find(mat_id);
  if (itr != d_BCDataArray.end()) {
    return itr->second.size();
  } else {
    itr = d_BCDataArray.find(-1);
    if (itr != d_BCDataArray.end()) {
      return itr->second.size();
    }
  }
  return 0;
}

BCGeomBase*
BCDataArray::getChild(int mat_id, int i) const
{
  auto itr = d_BCDataArray.find(mat_id);
  if (itr != d_BCDataArray.end()) {
    return itr->second[i];
  } else {
    itr = d_BCDataArray.find(-1);
    if (itr != d_BCDataArray.end()) {
      return itr->second[i];
    }
  }
  return nullptr;
}

void
BCDataArray::print() const
{
  bcDataArrayType::const_iterator bcda_itr;
  for (bcda_itr = d_BCDataArray.begin(); bcda_itr != d_BCDataArray.end();
       bcda_itr++) {
    BCDA_dbg << std::endl << "mat_id = " << bcda_itr->first << std::endl;
    BCDA_dbg << "Size of BCGeomBase vector = " << bcda_itr->second.size()
             << std::endl;
    for (auto i : bcda_itr->second) {
      BCDA_dbg << "BCGeometry Type = " << typeid(*i).name() << " " << i << std::endl;
      i->print();
    }
  }
}

} // namespace Uintah