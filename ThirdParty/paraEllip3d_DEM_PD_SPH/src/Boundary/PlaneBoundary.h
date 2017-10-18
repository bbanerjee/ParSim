#ifndef PLANE_BOUNDARY_H
#define PLANE_BOUNDARY_H

#include <Boundary/Boundary.h>
#include <Boundary/BoundaryContact.h>
#include <Boundary/BoundaryTangent.h>
#include <Boundary/BoundaryContainers.h>
#include <Boundary/BoundaryConditionCurve.h>
#include <Core/Geometry/Plane.h>
#include <Core/Math/Vec.h>
#include <Core/Types/RealTypes.h>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <iostream>
#include <utility>

namespace dem {

class DEMParticle; // forward declaration, only use pointer to class DEMParticle

///////////////////////////////////////
class PlaneBoundary : public Boundary
{
public:

  enum class BCType 
  {
    NONE = 0,
    TRACTION = 1,
    DISPLACEMENT = 2
  };

  PlaneBoundary(BoundaryID id = BoundaryID::NONE, 
                Boundary::BoundaryType tp = Boundary::BoundaryType::NONE, 
                EdgeCount en = 0)
    : Boundary(id, tp, en)
    , d_area(1)
    , d_direction(0)
    , d_position(0)
    , d_previousPosition(0)
    , d_velocity(0)
    , d_previousVelocity(0)
    , d_bcType(BCType::NONE)
    , d_tractionBC()
    , d_displacementBC()
  {
  }

  PlaneBoundary(Boundary::BoundaryType type, std::istream& ifs);
  PlaneBoundary(Boundary::BoundaryType type, BoundaryID id, const XMLProblemSpec& ps);
  PlaneBoundary(Boundary::BoundaryType type, BoundaryID id, const JsonProblemSpec& ps);

  REAL getArea() const override { return d_area; }
  Vec getDirection() const { return d_direction; }
  Vec getPosition() const override { return d_position; }
  Vec getVelocity() const override { return d_velocity; }
  Vec getPreviousPosition() const override { return d_previousPosition; }
  Vec getPreviousVelocity() const override { return d_previousVelocity; }

  void setArea(REAL area) override { d_area = area; }
  void setDirection(Vec dir) { d_direction = dir; }
  void setPosition(Vec pnt) override { d_position = pnt; }
  void setVelocity(Vec vel) override { d_velocity = vel; }

  REAL distanceToBdry(Vec pos) const
  {
    return dot((pos - d_position) , normalize(d_direction));
  }
  REAL distanceToBdry(Vec pos, Plane pn) const
  {
    return dot(pos - pn.getPosition() , normalize(pn.getDirection()));
  }

  double getTraction(double time) const
  {
    double traction = 0.0;
    if (d_bcType == BCType::TRACTION) {
      traction = d_tractionBC.getBCValue(time);
    }
    return traction;
  }

  double getDisplacement(double time) const
  {
    double disp = 0.0;
    if (d_bcType == BCType::DISPLACEMENT) {
      disp = d_displacementBC.getBCValue(time);
    }
    return disp;
  }

  Vec getAppliedTraction(REAL currentTime) const override
  {
    REAL traction_val = getTraction(currentTime);
    Vec traction_vec = d_direction*traction_val;
    return traction_vec;
  }

  void updatePositionAndVelocity(double currTime, double deltaT,
                                 double area, double mass) override;

  void updateUsingTraction(double deltaT, double mass, const Vec& traction);

  void updateUsingDisplacement(double deltaT, const Vec& disp, 
                               const Vec& dispRate);

  void print(std::ostream& os) override;
  void printContactInfo(std::ostream& os) override;

  void updateIsotropic(REAL sigma, REAL areaX, REAL areaY, REAL areaZ) override;
  void updateOdometer(REAL simga, REAL areaX, REAL areaY, REAL areaZ) override;
  void updateTriaxial(REAL simga, REAL areaX, REAL areaY, REAL areaZ) override;
  void updatePlaneStrain(REAL simga, REAL areaX, REAL areaY,
                         REAL areaZ) override;
  void updateTrueTriaxial(REAL simga, REAL areaX, REAL areaY, REAL areaZ,
                          REAL sigmaX, REAL sigmaY) override;
  void findBoundaryContacts(DEMParticlePArray& ptcls) override;
  void boundaryForce(BoundaryTangentArrayMap& boundaryTangentMap) override;

  friend std::ostream& 
  operator<<(std::ostream& os, const PlaneBoundary& plane)
  {
    os << " Direction: " << plane.d_direction;
    os << " Position: " << plane.d_position;
    os << " Velocity: " << plane.d_velocity;
    os << " BCType: " << static_cast<int>(plane.d_bcType) << "\n";
    if (plane.d_bcType == PlaneBoundary::BCType::TRACTION) {
      os << " Tractions: \n\t" << plane.d_tractionBC;
    } else if (plane.d_bcType == PlaneBoundary::BCType::DISPLACEMENT) {
      os << " Displacements: \n\t" << plane.d_displacementBC;
    }
    return os;
  }

private:
  REAL d_area;
  Vec d_direction;
  Vec d_position;
  Vec d_previousPosition;
  Vec d_velocity;
  Vec d_previousVelocity;

  BCType d_bcType;
  BoundaryConditionCurve<double> d_tractionBC;
  BoundaryConditionCurve<double> d_displacementBC;

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& boost::serialization::base_object<Boundary>(*this);
    ar& d_area;
    ar& d_direction;
    ar& d_position;
    ar& d_previousPosition;
    ar& d_velocity;
    ar& d_previousVelocity;
    ar& d_bcType;
    ar& d_tractionBC;
    ar& d_displacementBC;
  }

};

} // namespace dem ends

#endif
