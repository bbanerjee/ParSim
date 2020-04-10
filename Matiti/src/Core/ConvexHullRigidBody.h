/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015 Parresia Research Limited, New Zealand
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

#ifndef MATITI_CONVEX_HULL_RIGID_BODY_H
#define MATITI_CONVEX_HULL_RIGID_BODY_H

#include <Core/Geometry/Vector.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

namespace Matiti {

  class ConvexHullRigidBody
  {
  public:  

    friend std::ostream& operator<<(std::ostream& out, const Matiti::ConvexHullRigidBody& body);

  public:
   
    ConvexHullRigidBody();
    virtual ~ConvexHullRigidBody();

    void initialize(const int& id,
                    const std::vector<double>& mass,
                    const std::vector<double>& volume,
                    const std::vector<Uintah::Vector>& position,
                    const std::vector<Uintah::Vector>& velocity,
                    const double& vel_scale_fac,
                    const Uintah::Vector& bodyForce,
                    const Uintah::Vector& centerOfRotation,
                    const Uintah::Vector& angularVelocityOfRotation);

    /**
     * Get methods 
     */
    inline int id() const {return d_id;}

    const Uintah::Vector& centerOfMass() const {return d_com;}     
    inline double volume() const {return d_volume;}
    inline double mass()  const {return d_mass;}

    const Uintah::Vector& velocity() const {return d_vel;}

    const Uintah::Vector& bodyForce() const {return d_body_force;} 

    const Uintah::Vector& rotatingCoordCenter() const {return d_rot_center;} 
    const Uintah::Vector& rotatingCoordAngularVelocity() const {return d_rot_vel;}    

    const std::vector<Uintah::Vector>& getPositions() const
    {
      return d_positions;
    }
    const int getNumberOfPoints() const
    {
      return (int) d_positions.size();
    }

    /**
     * Set methods 
     */
    void setCenterOfMass(const Uintah::Vector& pos) 
    {
      d_com = pos;
    }
    void setPositions(const std::vector<Uintah::Vector>& positions) 
    {
      d_positions = positions;
    }
    void setVelocity(const Uintah::Vector& vel) {
      d_vel = vel;
    }

  private:

    int            d_id;

    double         d_mass;
    double         d_volume;

    Uintah::Vector d_com;     // Center of mass
    Uintah::Vector d_vel;     // Velocity

    Uintah::Vector d_body_force; // Gravity
    Uintah::Vector d_rot_center; // Rotating coord center
    Uintah::Vector d_rot_vel;    // Rotating coord angular velocity

    std::vector<Uintah::Vector> d_positions;  // Point positions
  };
} // end namespace

#endif
