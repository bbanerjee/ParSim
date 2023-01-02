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
#include <Core/Util/DebugStream.h>
#include <Core/Util/RWS.h>

#include <iostream>
#include <string>

namespace Uintah {

static DebugStream dbg("GeometryPieceFactory", false);

// Static class variable definition:
std::map<std::string, GeometryPieceP> GeometryPieceFactory::s_namedPieces;
std::vector<GeometryPieceP> GeometryPieceFactory::s_unnamedPieces;

void
GeometryPieceFactory::create(const ProblemSpecP& ps,
                             std::vector<GeometryPieceP>& objs) {
  for (ProblemSpecP child = ps->findBlock(); child != 0;
       child              = child->findNextBlock()) {
    std::string go_type = child->getNodeName();
    std::string go_label;

    if (!child->getAttribute("label", go_label)) {
      // "label" and "name" are both used... so check for "label"
      // first, and if it isn't found, then check for "name".
      child->getAttribute("name", go_label);
    }

    dbg << "---------------------------------------------------------------: "
           "go_label: "
        << go_label << "\n";

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
        dbg << "Creating new GeometryPiece: " << go_label
            << " (of type: " << go_type << ")\n";
      } else {
        if (referencedPiece != nullptr) {
          dbg << "Referencing GeometryPiece: " << go_label
              << " (of type: " << go_type << ")\n";
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
      dbg << "Creating non-labeled GeometryPiece of type '" << go_type << "'\n";
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
      newGeomPiece = std::make_shared<TriGeometryPiece>(child);
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
      continue;  // restart loop to avoid accessing name of empty object

    } else {
      // Perhaps it is a shell piece... let's find out:
      newGeomPiece = ShellGeometryFactory::create(child);

      if (newGeomPiece == nullptr) {
        if (Parallel::getMPIRank() == 0) {
          std::cerr << "WARNING: Unknown Geometry Piece Type "
                    << "(" << go_type << ")\n";
        }
        continue;  // restart loop to avoid accessing name of empty object
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

  }  // end for( child )
  dbg << "Done creating geometry objects\n";
}

void
GeometryPieceFactory::resetFactory() {
  s_unnamedPieces.clear();
  s_namedPieces.clear();
}

void
GeometryPieceFactory::resetGeometryPiecesOutput() {
  dbg << "resetGeometryPiecesOutput()\n";

  for (auto& piece : s_unnamedPieces) {
    piece->resetOutput();
    dbg << "  - Reset: " << piece->getName() << "\n";
  }

  for (auto& piece : s_namedPieces) {
    piece.second->resetOutput();
    dbg << "  - Reset: " << piece.second->getName() << "\n";
  }
}

}  // end namespace Uintah