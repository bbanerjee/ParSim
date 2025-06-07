/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#ifndef __ABAQUS_MESH_GEOMETRY_OBJECT_H__
#define __ABAQUS_MESH_GEOMETRY_OBJECT_H__

#include <Core/Geometry/Point.h>
#include <Core/GeometryPiece/SpecialGeomPiece.h>
#include <Core/Grid/Box.h>

#include <iostream>
#include <list>
#include <string>
#include <vector>

namespace Uintah {

/////////////////////////////////////////////////////////////////////////////
/*!
  \class AbaqusMeshGeometryPiece

  \brief Reads in a set of nodes and elements of a tetrahedral mesh and
         compute nodal masses and volume based on element connectivity.
         The input file is in Abaqus .inp ASCII format.

  \author Biswajit Banerjee

  The input tags look like this:
  \verbatim
    <abaqus_mesh>
      <file_name>file_name.inp</file_name>
      <scaling_factor> 0.5 </scaling_factor>
      <translation_vector> [-1.0, 0.0, 0.0] </translation_vector>
      <reflection_vector>  [1, -1, 1]       </reflection_vector>
      <axis_sequence>      [3, 1, 2]        </axis_sequence>
    </abaqus_mesh>
  \endverbatim
*/
/////////////////////////////////////////////////////////////////////////////

class AbaqusMeshGeometryPiece final : public SpecialGeomPiece {
 public:
  //////////////////////////////////////////////////////////////////////
  /*! \brief Constructor that takes a ProblemSpecP argument.
      It reads the xml input specification and builds a generalized box. */
  //////////////////////////////////////////////////////////////////////
  AbaqusMeshGeometryPiece(ProblemSpecP& ps);

  //////////////////////////////////////////////////////////////////////
  /*! Construct a box from a min/max point */
  //////////////////////////////////////////////////////////////////////
  AbaqusMeshGeometryPiece(const std::string& file_name);

  //////////
  // Destructor
  virtual ~AbaqusMeshGeometryPiece() = default;

  static const string TYPE_NAME;
  std::string
  getType() const override {
    return TYPE_NAME;
  }

  /// Make a clone
  GeometryPieceP
  clone() const override;

  //////////
  // Determines whether a point is inside the box.
  bool
  inside(const Point& p) const override;

  //////////
  //  Returns the bounding box surrounding the object.
  Box
  getBoundingBox() const override;

  void
  readPoints(int pid);

  unsigned int
  createPoints();

 protected:
  struct VolumeElement {
    VolumeElement(int id, std::vector<int> nodes) {
      id_     = id;
      node1_  = nodes[0];
      node2_  = nodes[1];
      node3_  = nodes[2];
      node4_  = nodes[3];
      volume_ = -1.0;
    }

    int id_;
    int node1_;
    int node2_;
    int node3_;
    int node4_;
    double volume_;
  };

  struct MeshNode {
    MeshNode(int id, double x, double y, double z) {
      id_     = id;
      x_      = x;
      y_      = y;
      z_      = z;
      volume_ = -1.0;
    }

    MeshNode(int id, double x, double y, double z, double vol) {
      id_     = id;
      x_      = x;
      y_      = y;
      z_      = z;
      volume_ = vol;
    }

    int id_;
    double x_;
    double y_;
    double z_;
    double volume_;
    std::vector<int> adjElements_;
  };

  struct SurfaceElement {
    SurfaceElement(int id, std::vector<int> nodes) {
      id_    = id;
      node1_ = nodes[0];
      node2_ = nodes[1];
      node3_ = nodes[2];
    }

    int id_;
    int node1_;
    int node2_;
    int node3_;
  };

  void
  readMeshNodesAndElements(const std::string& fileName);
  void
  readMeshNode(const std::string& inputLine, std::vector<MeshNode>& nodes);
  void
  readMeshVolumeElement(const std::string& inputLine,
                        std::vector<VolumeElement>& elements);
  void
  readMeshSurfaceElement(const std::string& inputLine,
                         std::vector<SurfaceElement>& elements);
  void
  computeElementVolumes(std::vector<MeshNode>& nodes,
                        std::vector<VolumeElement>& elements);
  void
  computeGaussPtVolumes(std::vector<MeshNode>& nodes,
                        std::vector<VolumeElement>& elements,
                        std::vector<MeshNode>& gaussPts);
  void
  computeNodalVolumes(std::vector<MeshNode>& nodes,
                      std::vector<VolumeElement>& elements,
                      std::map<int, int>& volElemMap);
  void
  findNodalAdjacentElements(std::vector<MeshNode>& nodes,
                            std::vector<VolumeElement>& elements);

 private:
  Box d_box;
  std::string d_fileName;
  double d_scalefac;
  Uintah::Vector d_translate;
  Uintah::Vector d_reflect;
  Uintah::IntVector d_axis;
  Uintah::Matrix3 d_rotate;
  bool d_use_gauss_pts;
  int d_num_gauss_pts;

  bool
  read_line(std::istream& is, Point& xmin, Point& xmax);
  void
  read_bbox(std::istream& source, Point& lowpt, Point& highpt) const;
  virtual void
  outputHelper(ProblemSpecP& ps) const;

  std::istream&
  getline_safer(std::istream& is, std::string& t);
};

}  // End namespace Uintah

#endif  // __ABAQUS_MESH_GEOMETRY_PIECE_H__
