// Function Definitions

#include <Peridynamics/PeriBond.h>
#include <Peridynamics/PeriParticle.h>
#include <Core/Util/Utility.h>

namespace pd {

//-------------------------------------------------------------------------
// Default Constructor
PeriBond::PeriBond()
{
  isAlive = true;
  isRecv = false; // not between recvPeriParticle
  weight = 0.0;
  initLength = 0.0;
  pt1 = nullptr;
  pt2 = nullptr;
}

// Overload Constructor
PeriBond::PeriBond(REAL tmp_length, PeriParticleP tmp_pt1,
                   PeriParticleP tmp_pt2)
  : initLength(tmp_length)
  , pt1(std::move(tmp_pt1))
  , pt2(std::move(tmp_pt2))
{
  isAlive = true;
  isRecv = false;
} // PeriBond()

// Destructor
PeriBond::~PeriBond()
{
  //		pt1 = NULL;
  //		pt2 = NULL;
}

//-------------------------------------------------------------------------
// Accessor Functions
bool
PeriBond::getIsAlive() const
{
  return isAlive;
}

REAL
PeriBond::getWeight() const
{
  return weight;
}

REAL
PeriBond::getInitLength() const
{
  return initLength;
}

REAL
PeriBond::volume(bool is_pt1) const
{
  if (is_pt1) {
    return pt1->volume();
  } else {
    return pt2->volume();
  }
}

dem::Vec
PeriBond::getXi(bool is_pt1) const
{
  if (is_pt1) {
    return (pt2->getInitPosition() - pt1->getInitPosition());
  } else {
    return (pt1->getInitPosition() - pt2->getInitPosition());
  }
} // end getXi()

dem::Vec
PeriBond::getEta(bool is_pt1) const
{
  if (is_pt1) {
    return (pt2->getInitPosition() + pt2->getDisplacement() -
            pt1->getInitPosition() - pt1->getDisplacement());
  } else {
    return (pt1->getInitPosition() + pt1->getDisplacement() -
            pt2->getInitPosition() - pt2->getDisplacement());
  }
} // end getEta()

dem::Vec
PeriBond::getEtaHalf(bool is_pt1, const REAL dt) const
{
  if (is_pt1) {
    return (pt2->getInitPosition() + pt2->getDisplacement() -
            0.5 * dt * pt2->getVelocityHalf() -
            (pt1->getInitPosition() + pt1->getDisplacement() -
             0.5 * dt * pt1->getVelocityHalf()));
  } else {
    return (pt1->getInitPosition() + pt1->getDisplacement() -
            0.5 * dt * pt1->getVelocityHalf() -
            (pt2->getInitPosition() + pt2->getDisplacement() -
             0.5 * dt * pt2->getVelocityHalf()));
  }
} // end getEtaHalf()

dem::Matrix
PeriBond::getMicroK(const bool is_pt1) const
{

  dem::Vec xi;
  REAL volume;
  if (is_pt1) {
    xi = (pt2->getInitPosition() - pt1->getInitPosition());
    volume = pt1->volume();
  } else {
    xi = (pt1->getInitPosition() - pt2->getInitPosition());
    volume = pt2->volume();
  }

  return (dyadicProduct(xi, xi) * volume * weight);

} // end getMicroK()

dem::Matrix
PeriBond::getMicroN(const bool is_pt1, const bool bondIsAlive) const
{

  dem::Vec xi;
  dem::Vec eta;
  REAL volume;
  if (bondIsAlive) {
    if (is_pt1) {
      xi = (pt2->getInitPosition() - pt1->getInitPosition());
      eta = (pt2->getInitPosition() + pt2->getDisplacement() -
             pt1->getInitPosition() - pt1->getDisplacement());
      volume = pt1->volume();
    } else {
      xi = (pt1->getInitPosition() - pt2->getInitPosition());
      eta = (pt1->getInitPosition() + pt1->getDisplacement() -
             pt2->getInitPosition() - pt2->getDisplacement());
      volume = pt2->volume();
    }
  } else {
    if (is_pt1) {
      xi = (pt2->getInitPosition() - pt1->getInitPosition());
      // eta = (pt2->getInitPosition()+pt1->getDisplacement() -
      // pt1->getInitPosition()-pt1->getDisplacement());
      eta = xi;
      volume = pt1->volume();
    } else {
      xi = (pt1->getInitPosition() - pt2->getInitPosition());
      // eta = (pt1->getInitPosition()+pt2->getDisplacement() -
      // pt2->getInitPosition()-pt2->getDisplacement());
      eta = xi;
      volume = pt2->volume();
    }
  }

  /*
  std::cout << "Bond: (P1,P2) = (" << pt1->getId() << "," << pt2->getId() << "): "
            << " u1 = " << pt1->getDisplacement()
            << " u2 = " << pt2->getDisplacement()
            << " xi = " << xi << " eta = " << eta << " vol = " << volume 
            << " weight = " << weight << "\n";
  */
  return (dyadicProduct(eta, xi) * volume * weight);

} // end getMicroN()

dem::Matrix
PeriBond::getMicroNHalf(const bool is_pt1, const bool isBondAlive,
                        const REAL dt) const
{

  dem::Vec xi;
  dem::Vec eta;
  REAL volume;
  if (isBondAlive) {
    if (is_pt1) {
      xi = (pt2->getInitPosition() - pt1->getInitPosition());
      eta = (pt2->getInitPosition() + pt2->getDisplacement() -
             0.5 * dt * pt2->getVelocityHalf() -
             (pt1->getInitPosition() + pt1->getDisplacement() -
              0.5 * dt * pt1->getVelocityHalf()));
      volume = pt1->volume();
    } else {
      xi = (pt1->getInitPosition() - pt2->getInitPosition());
      eta = (pt1->getInitPosition() + pt1->getDisplacement() -
             0.5 * dt * pt1->getVelocityHalf() -
             (pt2->getInitPosition() + pt2->getDisplacement() -
              0.5 * dt * pt2->getVelocityHalf()));
      volume = pt2->volume();
    }
  } else {
    if (is_pt1) {
      xi = (pt2->getInitPosition() - pt1->getInitPosition());
      // eta =
      //(pt2->getInitPosition()+pt2->getDisplacement()-0.5*dt*pt2->getVelocityHalf()
      //	 -
      //(pt1->getInitPosition()+pt1->getDisplacement()-0.5*dt*pt1->getVelocityHalf())
      //);
      eta = xi;
      volume = pt1->volume();
    } else {
      xi = (pt1->getInitPosition() - pt2->getInitPosition());
      // eta =
      //(pt1->getInitPosition()+pt1->getDisplacement()-0.5*dt*pt1->getVelocityHalf()
      //	 -
      //(pt2->getInitPosition()+pt2->getDisplacement()-0.5*dt*pt2->getVelocityHalf())
      //);
      eta = xi;
      volume = pt2->volume();
    }
  }

  return (dyadicProduct(eta, xi) * volume * weight);

} // end getMicroNHalf()

dem::Matrix
PeriBond::getMicroNDeltaU(const bool is_pt1, const bool isBondAlive,
                          const REAL dt) const
{

  dem::Vec xi;
  dem::Vec eta;
  REAL volume;
  if (isBondAlive) {
    if (is_pt1) {
      xi = (pt2->getInitPosition() - pt1->getInitPosition());
      eta = dt * (pt2->getVelocityHalf() - pt1->getVelocityHalf());
      volume = pt1->volume();
    } else {
      xi = (pt1->getInitPosition() - pt2->getInitPosition());
      eta = dt * (pt1->getVelocityHalf() - pt2->getVelocityHalf());
      volume = pt2->volume();
    }
    return (dyadicProduct(eta, xi) * volume * weight);
  }
  // commented out, because eta = 0.0 for this case, no contribution needs to be
  // added
  // else {
  //	if(is_pt1){
  //	    xi = (pt2->getInitPosition()-pt1->getInitPosition());
  //	    // eta = dt*(pt2->getVelocityHalf() - pt1->getVelocityHalf());
  //		eta = 0.0;
  //	    volume = pt1->volume();
  //	}
  //	else{
  //	    xi = (pt1->getInitPosition()-pt2->getInitPosition());
  //	    // eta = dt*(pt1->getVelocityHalf() - pt2->getVelocityHalf());
  //		eta = 0.0;
  //	    volume = pt2->volume();
  //	}
  //}
  return dem::Matrix(3, 3);

} // end getMicroNDeltaU()

//-------------------------------------------------------------------------
// Mutator Functions
void
PeriBond::setIsAlive(bool newisAlive)
{
  isAlive = newisAlive;
}

void
PeriBond::setWeight(REAL newweight)
{
  weight = newweight;
}

void
PeriBond::setInitLength(REAL newinitLength)
{
  initLength = newinitLength;
}

//-------------------------------------------------------------------------
// Utility Functions

REAL
PeriBond::calcCurrentLength()
{
  return vnormL2(pt1->getInitPosition() + pt1->getDisplacement() -
               pt2->getInitPosition() - pt2->getDisplacement());
}

void
PeriBond::checkIfAlive()
{
  if (getIsAlive()) {
    REAL bond_length = calcCurrentLength();

    REAL init_length = getInitLength();
    REAL stretch = (bond_length - init_length) / init_length;

    if (stretch >
          util::getParam<REAL>("bondStretchLimit") ||
        stretch < -2.0)
      setAliveFalse();
  }
}

} // end pd
