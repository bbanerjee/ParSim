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

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/GeometryPiece/AbaqusMeshGeometryPiece.h>
#include <Core/GeometryPiece/BoxGeometryPiece.h>
#include <Core/GeometryPiece/ConeGeometryPiece.h>
#include <Core/GeometryPiece/CorrugEdgeGeomPiece.h>
#include <Core/GeometryPiece/CylinderGeometryPiece.h>
#include <Core/GeometryPiece/DifferenceGeometryPiece.h>
#include <Core/GeometryPiece/EllipsoidGeometryPiece.h>
#include <Core/GeometryPiece/FileGeometryPiece.h>
#include <Core/GeometryPiece/GeometryPieceFactory.h>
#include <Core/GeometryPiece/IntersectionGeometryPiece.h>
#include <Core/GeometryPiece/NaaBoxGeometryPiece.h>
#include <Core/GeometryPiece/NullGeometryPiece.h>
#include <Core/GeometryPiece/ShellGeometryFactory.h>
#include <Core/GeometryPiece/ShellGeometryPiece.h>
#include <Core/GeometryPiece/SmoothCylGeomPiece.h>
#include <Core/GeometryPiece/SmoothSphereGeomPiece.h>
#include <Core/GeometryPiece/SphereGeometryPiece.h>
#include <Core/GeometryPiece/SphereMembraneGeometryPiece.h>
#include <Core/GeometryPiece/TorusGeometryPiece.h>
#include <Core/GeometryPiece/TriGeometryPiece.h>
#include <Core/GeometryPiece/UnionGeometryPiece.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Parallel/Parallel.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Util/DOUT.hpp>
#include <Core/Util/DebugStream.h>
#include <Core/Util/RWS.h>

#include <iostream>
#include <string>

namespace {
Uintah::Dout dout_gpf("GeometryPieceFactory",
                      "GeometryPieceFactory",
                      "GeometryPieceFactory debug stream",
                      false);
}

namespace Uintah {

void
GeometryPieceFactory::create(const ProblemSpecP& ps,
                             std::vector<GeometryPieceP>& objs,
                             const std::string& input_ups_dir)
{
  for (ProblemSpecP child = ps->findBlock(); child != 0;
       child              = child->findNextBlock()) {
    std::string go_type = child->getNodeName();
    std::string go_label;

    if (!child->getAttribute("label", go_label)) {
      // "label" and "name" are both used... so check for "label"
      // first, and if it isn't found, then check for "name".
      child->getAttribute("name", go_label);
    }

    DOUTR(dout_gpf,
          "---------------------------------------------------------------: "
          "go_label: "
            << go_label);

    if (go_label != "") {
      ProblemSpecP childBlock = child->findBlock();

      // See if there is any data for this node (that is not in a sub-block)
      std::string data = child->getNodeValue();
      remove_lt_white_space(data);

      // Lookup in table to see if this piece has already be named...
      GeometryPieceP referencedPiece = s_namedPieces[go_label];

      // If it has a childBlock or data, then it is not just a reference.
      bool goHasInfo = childBlock || data != "";

      if (referencedPiece != nullptr && goHasInfo) {
        std::ostringstream out;
        out << "Error: Duplicate GeomPiece definition not allowed.\n"
            << " GeometryPiece (" << go_type << ")"
            << " labeled: '" << go_label << "' has already been specified."
            << " You cannot change its values.\n"
            << " Please just reference the original by only"
            << " using the label (no values)\n";
        throw ProblemSetupException(out.str(), __FILE__, __LINE__);
      }

      if (goHasInfo) {
        DOUTR(dout_gpf,
              "Creating new GeometryPiece: " << go_label << " (of type: "
                                             << go_type << ")");
      } else {
        if (referencedPiece != nullptr) {
          DOUTR(dout_gpf,
                "Referencing GeometryPiece: " << go_label << " (of type: "
                                              << go_type << ")");
          objs.push_back(referencedPiece);
        } else {
          std::cout << "Error... couldn't find the referenced GeomPiece: "
                    << go_label << " (" << go_type << ")\n";
          throw ProblemSetupException(
            "Referenced GeomPiece does not exist", __FILE__, __LINE__);
        }

        // Verify that the referenced piece is of the same type as
        // the originally created piece.
        if (referencedPiece->getType() != go_type) {
          std::cout << "Error... the referenced GeomPiece: " << go_label << " ("
                    << referencedPiece->getType() << "), "
                    << "is not of the same type as this new object: '"
                    << go_type << "'!\n";
          throw ProblemSetupException(
            "Referenced GeomPiece is not of the same type as original",
            __FILE__,
            __LINE__);
        }
        continue;
      }

    } else {
      DOUTR(dout_gpf,
            "Creating non-labeled GeometryPiece of type '" << go_type << "'");
    }

    GeometryPieceP newGeomPiece = nullptr;

    if (go_type == BoxGeometryPiece::TYPE_NAME) {
      newGeomPiece = std::make_shared<BoxGeometryPiece>(child);
    } else if (go_type == NaaBoxGeometryPiece::TYPE_NAME) {
      newGeomPiece = std::make_shared<NaaBoxGeometryPiece>(child);
    } else if (go_type == SphereGeometryPiece::TYPE_NAME) {
      newGeomPiece = std::make_shared<SphereGeometryPiece>(child);
    } else if (go_type == CylinderGeometryPiece::TYPE_NAME) {
      newGeomPiece = std::make_shared<CylinderGeometryPiece>(child);
    } else if (go_type == ConeGeometryPiece::TYPE_NAME) {
      newGeomPiece = std::make_shared<ConeGeometryPiece>(child);
    } else if (go_type == EllipsoidGeometryPiece::TYPE_NAME) {
      newGeomPiece = std::make_shared<EllipsoidGeometryPiece>(child);
    } else if (go_type == TorusGeometryPiece::TYPE_NAME) {
      newGeomPiece = std::make_shared<TorusGeometryPiece>(child);
    } else if (go_type == TriGeometryPiece::TYPE_NAME) {
      newGeomPiece = std::make_shared<TriGeometryPiece>(child, input_ups_dir);
    } else if (go_type == UnionGeometryPiece::TYPE_NAME) {
      newGeomPiece = std::make_shared<UnionGeometryPiece>(child);
    } else if (go_type == DifferenceGeometryPiece::TYPE_NAME) {
      newGeomPiece = std::make_shared<DifferenceGeometryPiece>(child);
    } else if (go_type == IntersectionGeometryPiece::TYPE_NAME) {
      newGeomPiece = std::make_shared<IntersectionGeometryPiece>(child);
    } else if (go_type == SmoothCylGeomPiece::TYPE_NAME) {
      newGeomPiece = std::make_shared<SmoothCylGeomPiece>(child);
    } else if (go_type == SphereMembraneGeometryPiece::TYPE_NAME) {
      newGeomPiece = std::make_shared<SphereMembraneGeometryPiece>(child);
    } else if (go_type == SmoothSphereGeomPiece::TYPE_NAME) {
      newGeomPiece = std::make_shared<SmoothSphereGeomPiece>(child);
    } else if (go_type == CorrugEdgeGeomPiece::TYPE_NAME) {
      newGeomPiece = std::make_shared<CorrugEdgeGeomPiece>(child);
    } else if (go_type == FileGeometryPiece::TYPE_NAME) {
      newGeomPiece = std::make_shared<FileGeometryPiece>(child);
    } else if (go_type == AbaqusMeshGeometryPiece::TYPE_NAME) {
      newGeomPiece = std::make_shared<AbaqusMeshGeometryPiece>(child);
    } else if (go_type == NullGeometryPiece::TYPE_NAME) {
      newGeomPiece = std::make_shared<NullGeometryPiece>(child);
    } else if (go_type == "res" || go_type == "velocity" ||
               go_type == "temperature" || go_type == "comment" ||
               go_type == "density" || go_type == "pressure" ||
               go_type == "scalar" || go_type == "color" ||
               go_type == "volumeFraction") {
      // Ignoring.
      continue; // restart loop to avoid accessing name of empty object

    } else {
      // Perhaps it is a shell piece... let's find out:
      newGeomPiece = ShellGeometryFactory::create(child);

      if (newGeomPiece == nullptr) {
        if (Parallel::getMPIRank() == 0) {
          std::cerr << "WARNING: Unknown Geometry Piece Type "
                    << "(" << go_type << ")\n";
        }
        continue; // restart loop to avoid accessing name of empty object
      }
    }
    // Look for the "name" of the object.  (Can also be referenced as "label").
    string name;
    if (child->getAttribute("name", name)) {
      newGeomPiece->setName(name);
    } else if (child->getAttribute("label", name)) {
      newGeomPiece->setName(name);
    }
    if (name != "") {
      s_namedPieces[name] = newGeomPiece;
    } else {
      s_unnamedPieces.push_back(newGeomPiece);
    }

    objs.push_back(newGeomPiece);

  } // end for( child )
  DOUTR(dout_gpf, "Done creating geometry objects");
}

void
GeometryPieceFactory::resetFactory()
{
  s_unnamedPieces.clear();
  s_namedPieces.clear();
}

bool
GeometryPieceFactory::foundInsidePoints(const std::string geomPieceName,
                                        const int patchID)
{
  typedef std::map<int, std::vector<Point>> PatchIDInsidePointsMapT;

  if (s_insidePointsMap.find(geomPieceName) != s_insidePointsMap.end()) {
    // we found this geometry, lets see if we find this patch
    PatchIDInsidePointsMapT& thisPatchIDInsidePoints =
      s_insidePointsMap[geomPieceName];

    if (thisPatchIDInsidePoints.find(patchID) !=
        thisPatchIDInsidePoints.end()) {
      // we found this patch ID
      return true;
    }
  }
  return false;
}

const std::vector<Point>&
GeometryPieceFactory::getInsidePoints(const Uintah::Patch* const patch)
{
  typedef std::map<std::string, GeometryPieceP> NameGeomPiecesMapT;
  const int patchID = patch->getID();

  if (s_allInsidePointsMap.find(patchID) != s_allInsidePointsMap.end()) {
    return s_allInsidePointsMap[patchID];
  }

  // loop over all geometry objects
  std::vector<Point> allInsidePoints;
  NameGeomPiecesMapT::iterator geomIter;

  for (geomIter = s_namedPieces.begin(); geomIter != s_namedPieces.end();
       ++geomIter) {
    GeometryPieceP geomPiece = geomIter->second;

    const string geomName = geomPiece->getName();
    const std::vector<Point>& thisGeomInsidePoints =
      getInsidePoints(geomName, patch);

    allInsidePoints.insert(allInsidePoints.end(),
                           thisGeomInsidePoints.begin(),
                           thisGeomInsidePoints.end());
  }

  s_allInsidePointsMap.insert(
    std::pair<int, std::vector<Point>>(patchID, allInsidePoints));
  return s_allInsidePointsMap[patchID];
}

const std::vector<Point>&
GeometryPieceFactory::getInsidePoints(const std::string geomPieceName,
                                      const Uintah::Patch* const patch)
{
  typedef std::map<int, std::vector<Point>> PatchIDInsidePointsMapT;

  const int patchID = patch->getID();

  if (s_insidePointsMap.find(geomPieceName) != s_insidePointsMap.end()) {
    // we found this geometry, lets see if we find this patch
    PatchIDInsidePointsMapT& patchIDInsidePoints =
      s_insidePointsMap[geomPieceName];

    if (patchIDInsidePoints.find(patchID) == patchIDInsidePoints.end()) {
      // we did not find this patch. check if there are any points in this
      // patch.
      GeometryPieceP geomPiece = s_namedPieces[geomPieceName];
      std::vector<Point> insidePoints;

      for (Uintah::CellIterator iter(patch->getCellIterator()); !iter.done();
           iter++) {
        IntVector iCell = *iter;

        Point p             = patch->getCellPosition(iCell);
        const bool isInside = geomPiece->inside(p, false);

        if (isInside) {
          insidePoints.push_back(p);
        }
      }
      patchIDInsidePoints.insert(
        std::pair<int, std::vector<Point>>(patch->getID(), insidePoints));
    }
  } else {
    // if we did not find this geometry piece
    std::vector<Point> insidePoints;
    std::map<int, std::vector<Point>> patchIDInsidePoints;
    GeometryPieceP geomPiece = s_namedPieces[geomPieceName];

    for (Uintah::CellIterator iter(patch->getCellIterator()); !iter.done();
         iter++) {
      IntVector iCell = *iter;

      Point p             = patch->getCellPosition(iCell);
      const bool isInside = geomPiece->inside(p, false);

      if (isInside) {
        insidePoints.push_back(p);
      }
    }
    patchIDInsidePoints.insert(
      std::pair<int, std::vector<Point>>(patch->getID(), insidePoints));
    s_insidePointsMap.insert(std::pair<std::string, PatchIDInsidePointsMapT>(
      geomPieceName, patchIDInsidePoints));
  }
  // at this point, we can GUARANTEE that there is a vector of points associated
  // with this geometry and patch. This vector could be empty.
  return s_insidePointsMap[geomPieceName][patchID];
}

void
GeometryPieceFactory::findInsidePoints(const Uintah::Patch* const patch)
{
  typedef std::map<std::string, GeometryPieceP> NameGeomPiecesMapT;
  typedef std::map<int, std::vector<Point>> PatchIDInsidePointsMapT;

  const int patchID = patch->getID();

  // loop over all geometry objects
  NameGeomPiecesMapT::iterator geomIter;
  for (geomIter = s_namedPieces.begin(); geomIter != s_namedPieces.end();
       ++geomIter) {

    std::vector<Point> insidePoints;
    std::map<int, std::vector<Point>> patchIDInsidePoints;

    GeometryPieceP geomPiece   = geomIter->second;
    const std::string geomName = geomPiece->getName();

    // check if we already found the inside points for this patch
    if (foundInsidePoints(geomName, patchID)) {
      continue;
    }

    for (Uintah::CellIterator iter(patch->getCellIterator()); !iter.done();
         iter++) {
      IntVector iCell = *iter;

      Point p             = patch->getCellPosition(iCell);
      const bool isInside = geomPiece->inside(p, false);

      if (isInside) {
        insidePoints.push_back(p);
      }
    }

    patchIDInsidePoints.insert(
      std::pair<int, std::vector<Point>>(patch->getID(), insidePoints));
    s_insidePointsMap.insert(std::pair<std::string, PatchIDInsidePointsMapT>(
      geomName, patchIDInsidePoints));
  }
}

const std::map<std::string, GeometryPieceP>&
GeometryPieceFactory::getNamedGeometryPieces()
{
  return s_namedPieces;
}

//  Recursively search through the problem spec for geometry pieces
//  that have already been created.  This is tricky and confusing code. Consider
//  two cases: Case A is easy but with Case B the geom_object must be
//  recursively searched.
//
//   CASE A
//         <geom_object>
//           <box label="mpm_box">
//             <min>[1, 1, 1]</min>
//             <max>[1.5, 1.5, 1.5]</max>
//           </box>
//         </geom_object>
//   CASE B
//         <geom_object>
//           <difference>                    << start recursive search
//             <box label="domain">
//               <min>[-1, -1, -1]</min>
//               <max>[4, 4, 4]</max>
//             </box>
//             <box label="mpm_box"/>         << no childBlock or goLabel
//           </difference>
//         </geom_object>
//
//   returns a negative integer if any of the geom pieces was not found.
//   returns a positive integer if all of the geom pieces were found.
//
//______________________________________________________________________

int
GeometryPieceFactory::geometryPieceExists(const ProblemSpecP& ps,
                                          const bool isTopLevel /* true */)
{

  int nFoundPieces = 0;
  for (ProblemSpecP child = ps->findBlock(); child != nullptr;
       child              = child->findNextBlock()) {

    bool hasChildBlock = false;
    if (child->findBlock()) {
      hasChildBlock = true;
    }

    string go_label;

    // search for either a label or name.
    if (!child->getAttribute("label", go_label)) {
      child->getAttribute("name", go_label);
    }

    //
    if (go_label == "") {

      if (hasChildBlock) { // This could be either a <difference> or
                           // <intersection > node, dig deeper
        nFoundPieces += geometryPieceExists(child, false);
      }
      continue;
    }

    // Is this child a geometry piece
    GeometryPieceP referencedPiece = s_namedPieces[go_label];

    if (referencedPiece.get() != nullptr) {
      nFoundPieces += 1;
      continue;
    }

    // Does the child have the spec of a geom_piece?
    // See if there is any data for this node (that is not in a sub-block)
    // If the spec exists then the geom_piece doesn't exist
    string data = child->getNodeValue();
    remove_lt_white_space(data);

    bool has_go_spec = (hasChildBlock || data != "");
    if (has_go_spec) {
      nFoundPieces -= INT_MAX;
    }

    if (isTopLevel) {
      break;
    }
  }

  return nFoundPieces;
}

//------------------------------------------------------------------

void
GeometryPieceFactory::resetGeometryPiecesOutput()
{
  DOUTR(dout_gpf,
        "resetGeometryPiecesOutput() unnamedPieces.size()"
          << s_unnamedPieces.size() << " namedPieces.end() "
          << s_namedPieces.size());

  for (unsigned int pos = 0; pos < s_unnamedPieces.size(); pos++) {
    s_unnamedPieces[pos]->resetOutput();
    DOUTR(dout_gpf, "  - Reset: " << s_unnamedPieces[pos]->getName());
  }

  std::map<std::string, GeometryPieceP>::const_iterator iter =
    s_namedPieces.begin();

  while (iter != s_namedPieces.end()) {
    DOUTR(dout_gpf, "  - Reset: " << iter->second->getName());
    iter->second->resetOutput();
    iter++;
  }
}

} // end namespace Uintah