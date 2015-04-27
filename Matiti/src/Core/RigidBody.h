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

#ifndef MATITI_RIGID_BODY_H
#define MATITI_RIGID_BODY_H

#include <Geometry/Vector3D.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

#include <string>
#include <fstream>
#include <iostream>
#include <map>

namespace Matiti {

  class RigidBody
  {
  public:  

    friend std::ostream& operator<<(std::ostream& out, const Matiti::RigidBody& body);

  public:
   
    RigidBody();
    virtual ~RigidBody();

    void initialize(Uintah::ProblemSpecP& ps);

    /**
     * Get/set methods 
     */
    inline int id() const {return d_id;}
    inline void id(const int& id) {d_id = id;}

    inline double density() const {return d_density;}
    const Vector3D& centerOfMass() const {return d_com;}     
    inline double radius() const {return d_radius};
    inline double volume() const {return d_volume};
    inline double mass()  const {return d_mass};
    const Vector3D& radiusOfGyration() const {return d_rog_sq;}  

    const Vector3D& initialVelocity const {return d_init_vel;}
    const Vector3D& initialAcceleration const {return d_init_acc;}
    const Vector3D& initialAngularVelocity const {return d_init_ang_vel;}
    const Vector3D& initialAngularAcceleration const {return d_init_ang_acc;}

    const Vector3D& externalForce() const {return d_ext_force;}
    const Vector3D& externalTorque() const {return d_ext_torque;}

    const Vector3D& bodyForce() const {return d_body_force;} 

    const Vector3D& rotatingCoordCenter() const {return d_rot_center;} 
    const Vector3D& rotatingCoordAngularVelocity() const {return d_rot_vel;}    

  private:

    int d_id;
    double d_density;
    Vector3D d_com;     // Center of mass
    double d_radius;
    double d_volume;
    double d_mass;
    Vector3D d_rog_sq;  // Radius of gyration

    Vector3D d_init_vel;
    Vector3D d_init_acc;
    Vector3D d_init_ang_vel;
    Vector3D d_init_ang_acc;

    Vector3D d_ext_force;
    Vector3D d_ext_torque;

    Vector3D d_body_force; // Gravity

    Vector3D d_rot_center; // Rotating coord center
    Vector3D d_rot_vel;    // Rotating coord angular velocity

    void computeInertialProperties();
  };
} // end namespace

#endif
