#include <DiscreteElements/Boundary.h>
#include <DiscreteElements/Particle.h>
// use both pointer to and variable of class Particle

namespace dem {

Boundary::Boundary(std::size_t tp, std::ifstream &ifs) {
  type = tp;
  ifs >> extraNum;
  ifs >> id;
}

void Boundary::print(std::ostream &os) {
  os << std::endl << std::setw(OWID) << type << std::setw(OWID) << extraNum
     << std::endl << std::setw(OWID) << id;
}

void Boundary::printContactInfo(std::ostream &os) {
  os << std::setw(OWID) << id << std::endl << std::setw(OWID)
     << contactInfo.size() << std::endl << std::setw(OWID) << "pos_x"
     << std::setw(OWID) << "pos_y" << std::setw(OWID) << "pos_z"
     << std::setw(OWID) << "normal_x" << std::setw(OWID) << "normal_y"
     << std::setw(OWID) << "normal_z" << std::setw(OWID) << "tangt_x"
     << std::setw(OWID) << "tangt_y" << std::setw(OWID) << "tangt_z"
     << std::setw(OWID) << "pentr" << std::endl;

  for (std::vector<BdryContact>::iterator it = contactInfo.begin();
       it != contactInfo.end(); ++it)
    it->print(os);
}

void Boundary::clearStatForce() {
  contactNum = 0;
  normal = 0;
  tangt = 0;
  penetr = 0;
}

void Boundary::updateStatForce() {
  clearStatForce();
  contactNum = contactInfo.size();
  for (std::vector<BdryContact>::iterator it = contactInfo.begin();
       it != contactInfo.end(); ++it) {
    normal += it->normal;
    tangt += it->tangt;
    penetr += it->penetr;
  }
  if (contactNum != 0)
    penetr /= contactNum;
}

void Boundary::clearContactInfo() {
  possParticle.clear();
  contactInfo.clear();
}

planeBoundary::planeBoundary(std::size_t tp, std::ifstream &ifs)
    : Boundary(tp, ifs) {
  REAL dx, dy, dz, px, py, pz;
  ifs >> dx >> dy >> dz >> px >> py >> pz;
  direc = Vec(dx, dy, dz);
  point = Vec(px, py, pz);
  prevPoint = point;
  prevVeloc = veloc = 0;
  for (std::size_t i = 0; i < extraNum; ++i) {
    ifs >> dx >> dy >> dz >> px >> py >> pz;
    extraEdge.push_back(Plane(Vec(dx, dy, dz), Vec(px, py, pz)));
  }
}

void planeBoundary::print(std::ostream &os) {
  Boundary::print(os);
  os << std::setw(OWID) << direc.getX() << std::setw(OWID) << direc.getY()
     << std::setw(OWID) << direc.getZ() << std::setw(OWID) << point.getX()
     << std::setw(OWID) << point.getY() << std::setw(OWID) << point.getZ()
     << std::endl;

  for (std::vector<Plane>::iterator et = extraEdge.begin();
       et != extraEdge.end(); ++et)
    os << std::setw(OWID) << " " << std::setw(OWID) << et->getDirec().getX()
       << std::setw(OWID) << et->getDirec().getY() << std::setw(OWID)
       << et->getDirec().getZ() << std::setw(OWID) << et->getPoint().getX()
       << std::setw(OWID) << et->getPoint().getY() << std::setw(OWID)
       << et->getPoint().getZ() << std::endl;
}

void planeBoundary::printContactInfo(std::ostream &os) {
  Boundary::printContactInfo(os);
  os << std::setw(OWID) << " " << std::setw(OWID) << " " << std::setw(OWID)
     << " " << std::setw(OWID) << normal.getX() << std::setw(OWID)
     << normal.getY() << std::setw(OWID) << normal.getZ() << std::setw(OWID)
     << tangt.getX() << std::setw(OWID) << tangt.getY() << std::setw(OWID)
     << tangt.getZ() << std::setw(OWID) << penetr << std::endl << std::endl;
  ;
}

void planeBoundary::findBdryContact(ParticlePArray &ptcls) {
  possParticle.clear();
  contactInfo.clear();
  clearStatForce();

  for (auto it = ptcls.begin(); it != ptcls.end(); ++it) {
    if ((*it)->getType() ==
        0) { // only process free particles, excluding type 5
      REAL dist = distanceToBdry((*it)->getCurrPos());
      if (dist < 0 && fabs(dist) <= (*it)->getA()) {
        bool inside = true;
        for (std::vector<Plane>::iterator et = extraEdge.begin();
             et != extraEdge.end(); ++et) {
          REAL eDist = distanceToBdry((*it)->getCurrPos(), (*et));
          if (eDist >= 0) {
            inside = false;
            break;
          }
        }
        if (inside)
          possParticle.push_back(*it);
      }
    }
  }
}

void planeBoundary::boundaryForce(
    std::map<std::size_t, std::vector<BoundaryTgt> > &boundaryTgtMap) {
  // for each plane boundary, define a temparory variable vtmp to use,
  // better than define a member variable which needs to be cleared.
  // and vtmp is initialized as empty in each iteration.
  std::vector<BoundaryTgt> vtmp;

  // for each possible boundary particle
  for (auto it = possParticle.begin(); it != possParticle.end(); ++it)
    (*it)->planeRBForce(this, boundaryTgtMap, vtmp);

  // checkout tangential forces and displacements after each particle is
  // processed
  boundaryTgtMap[this->id] = vtmp;

  updateStatForce();
}

void planeBoundary::updateIsotropic(REAL sigma, REAL areaX, REAL areaY,
                                    REAL areaZ) {

  //REAL forceDamp = dem::Parameter::getSingleton().parameter["forceDamp"];
  //REAL massScale = dem::Parameter::getSingleton().parameter["massScale"];
  //REAL mass = dem::Parameter::getSingleton().parameter["boundaryMass"];
  REAL boundaryRate = dem::Parameter::getSingleton().parameter["boundaryRate"];
  REAL topSpeedup = dem::Parameter::getSingleton().parameter["topSpeedup"];
  REAL tol = dem::Parameter::getSingleton().parameter["tractionErrorTol"];
  //REAL atf = forceDamp * 2;

  REAL vel, pos;
  switch (id) {
  case 1:
    if (fabs(normal.getX() / areaX + sigma) / sigma > tol) {
      vel = ((normal.getX() + sigma * areaX) > 0 ? 1 : -1) * boundaryRate;
      //vel = prevVeloc.getX() * (2-atf) / (2+atf) + (normal.getX() + sigma *
      //areaX) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.getX() + vel * timeStep;
      setVeloc(Vec(vel, getVeloc().getY(), getVeloc().getZ()));
      setPoint(Vec(pos, getPoint().getY(), getPoint().getZ()));
    }
    break;
  case 2:
    if (fabs(normal.getX() / areaX - sigma) / sigma > tol) {
      vel = ((normal.getX() - sigma * areaX) > 0 ? 1 : -1) * boundaryRate;
      //vel = prevVeloc.getX() * (2-atf) / (2+atf) + (normal.getX() - sigma *
      //areaX) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.getX() + vel * timeStep;
      setVeloc(Vec(vel, getVeloc().getY(), getVeloc().getZ()));
      setPoint(Vec(pos, getPoint().getY(), getPoint().getZ()));
    }
    break;
  case 3:
    if (fabs(normal.getY() / areaY + sigma) / sigma > tol) {
      vel = ((normal.getY() + sigma * areaY) > 0 ? 1 : -1) * boundaryRate;
      //vel = prevVeloc.getY() * (2-atf) / (2+atf) + (normal.getY() + sigma *
      //areaY) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.getY() + vel * timeStep;
      setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ()));
      setPoint(Vec(getPoint().getX(), pos, getPoint().getZ()));
    }
    break;
  case 4:
    if (fabs(normal.getY() / areaY - sigma) / sigma > tol) {
      vel = ((normal.getY() - sigma * areaY) > 0 ? 1 : -1) * boundaryRate;
      //vel = prevVeloc.getY() * (2-atf) / (2+atf) + (normal.getY() - sigma *
      //areaY) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.getY() + vel * timeStep;
      setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ()));
      setPoint(Vec(getPoint().getX(), pos, getPoint().getZ()));
    }
    break;
  case 5:
    if (fabs(normal.getZ() / areaZ + sigma) / sigma > tol) {
      vel = ((normal.getZ() + sigma * areaZ) > 0 ? 1 : -1) * boundaryRate;
      //vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() + sigma *
      //areaZ) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.getZ() + vel * timeStep;
      setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel));
      setPoint(Vec(getPoint().getX(), getPoint().getY(), pos));
    }
    break;
  case 6:
    if (fabs(normal.getZ() / areaZ - sigma) / sigma > tol) {
      vel = ((normal.getZ() - sigma * areaZ) > 0 ? 1 : -1) * boundaryRate;
      if (normal.getZ() == 0)
        vel = -boundaryRate * topSpeedup;
      //vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() - sigma *
      //areaZ) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.getZ() + vel * timeStep;
      setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel));
      setPoint(Vec(getPoint().getX(), getPoint().getY(), pos));
    }
    break;
  }
  prevPoint = point;
  prevVeloc = veloc;
}

void planeBoundary::updateOdometer(REAL sigma, REAL areaX, REAL areaY,
                                   REAL areaZ) {

  //REAL forceDamp = dem::Parameter::getSingleton().parameter["forceDamp"];
  //REAL massScale = dem::Parameter::getSingleton().parameter["massScale"];
  //REAL mass = dem::Parameter::getSingleton().parameter["boundaryMass"];
  REAL boundaryRate = dem::Parameter::getSingleton().parameter["boundaryRate"];
  REAL tol = dem::Parameter::getSingleton().parameter["tractionErrorTol"];
  //REAL atf = forceDamp * 2;

  REAL vel, pos;
  switch (id) {
  case 5:
    if (fabs(normal.getZ() / areaZ + sigma) / sigma > tol) {
      vel = ((normal.getZ() + sigma * areaZ) > 0 ? 1 : -1) * boundaryRate;
      //vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() + sigma *
      //areaZ) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.getZ() + vel * timeStep;
      setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel));
      setPoint(Vec(getPoint().getX(), getPoint().getY(), pos));
    }
    break;
  case 6:
    if (fabs(normal.getZ() / areaZ - sigma) / sigma > tol) {
      vel = ((normal.getZ() - sigma * areaZ) > 0 ? 1 : -1) * boundaryRate;
      //vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() - sigma *
      //areaZ) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.getZ() + vel * timeStep;
      setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel));
      setPoint(Vec(getPoint().getX(), getPoint().getY(), pos));
    }
    break;
  }
  prevPoint = point;
  prevVeloc = veloc;
}

void planeBoundary::updateTriaxial(REAL sigma, REAL areaX, REAL areaY,
                                   REAL areaZ) {
  std::size_t triaxialType = static_cast<std::size_t>(
      dem::Parameter::getSingleton().parameter["triaxialType"]);
  std::size_t unloadStep = static_cast<std::size_t>(
      dem::Parameter::getSingleton().parameter["unloadStep"]);

  //REAL forceDamp = dem::Parameter::getSingleton().parameter["forceDamp"];
  //REAL massScale = dem::Parameter::getSingleton().parameter["massScale"];
  //REAL mass = dem::Parameter::getSingleton().parameter["boundaryMass"];
  REAL boundaryRate = dem::Parameter::getSingleton().parameter["boundaryRate"];
  REAL tol = dem::Parameter::getSingleton().parameter["tractionErrorTol"];
  //REAL atf = forceDamp * 2;

  REAL vel = 0.0, pos = 0.0;
  switch (id) {
  case 1:
    if (fabs(normal.getX() / areaX + sigma) / sigma > tol) {
      vel = ((normal.getX() + sigma * areaX) > 0 ? 1 : -1) * boundaryRate;
      //vel = prevVeloc.getX() * (2-atf) / (2+atf) + (normal.getX() + sigma *
      //areaX) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.getX() + vel * timeStep;
      setVeloc(Vec(vel, getVeloc().getY(), getVeloc().getZ()));
      setPoint(Vec(pos, getPoint().getY(), getPoint().getZ()));
    }
    break;
  case 2:
    if (fabs(normal.getX() / areaX - sigma) / sigma > tol) {
      vel = ((normal.getX() - sigma * areaX) > 0 ? 1 : -1) * boundaryRate;
      //vel = prevVeloc.getX() * (2-atf) / (2+atf) + (normal.getX() - sigma *
      //areaX) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.getX() + vel * timeStep;
      setVeloc(Vec(vel, getVeloc().getY(), getVeloc().getZ()));
      setPoint(Vec(pos, getPoint().getY(), getPoint().getZ()));
    }
    break;
  case 3:
    if (fabs(normal.getY() / areaY + sigma) / sigma > tol) {
      vel = ((normal.getY() + sigma * areaY) > 0 ? 1 : -1) * boundaryRate;
      //vel = prevVeloc.getY() * (2-atf) / (2+atf) + (normal.getY() + sigma *
      //areaY) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.getY() + vel * timeStep;
      setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ()));
      setPoint(Vec(getPoint().getX(), pos, getPoint().getZ()));
    }
    break;
  case 4:
    if (fabs(normal.getY() / areaY - sigma) / sigma > tol) {
      vel = ((normal.getY() - sigma * areaY) > 0 ? 1 : -1) * boundaryRate;
      //vel = prevVeloc.getY() * (2-atf) / (2+atf) + (normal.getY() - sigma *
      //areaY) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.getY() + vel * timeStep;
      setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ()));
      setPoint(Vec(getPoint().getX(), pos, getPoint().getZ()));
    }
    break;
  case 5:
    if (triaxialType == 1)
      vel = boundaryRate;
    else if (triaxialType == 2) {
      if (iteration <= unloadStep) // loading
        vel = boundaryRate;
      else if (iteration > unloadStep &&
               fabs(normal.getZ() / areaZ) >= 1.5 * sigma &&
               iteration <= 1.5 * unloadStep) // unloading
        vel = -boundaryRate;
      else if (iteration >
               1.5 * unloadStep) // reloading. Note there are invalid loops if
                                 // stress < and iteration <, it is OK.
        vel = boundaryRate;
    } else {
      vel = boundaryRate;
    }
    //vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() + sigma *
    //areaZ) / mass * timeStep * 2 / (2 + atf);
    pos = prevPoint.getZ() + vel * timeStep;
    setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel));
    setPoint(Vec(getPoint().getX(), getPoint().getY(), pos));
    break;
  case 6:
    if (triaxialType == 1)
      vel = -boundaryRate;
    else if (triaxialType == 2) {
      if (iteration <= unloadStep) // loading
        vel = -boundaryRate;
      else if (iteration > unloadStep &&
               fabs(normal.getZ() / areaZ) >= 1.5 * sigma &&
               iteration <= 1.5 * unloadStep) // unloading
        vel = boundaryRate;
      else if (iteration >
               1.5 * unloadStep) // reloading. Note there are invalid loops if
                                 // stress < and iteration <, it is OK.
        vel = -boundaryRate;
    } else {
      vel = boundaryRate;
    }
    //vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() - sigma *
    //areaZ) / mass * timeStep * 2 / (2 + atf);
    pos = prevPoint.getZ() + vel * timeStep;
    setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel));
    setPoint(Vec(getPoint().getX(), getPoint().getY(), pos));
    break;
  }
  prevPoint = point;
  prevVeloc = veloc;
}

void planeBoundary::updatePlaneStrain(REAL sigma, REAL areaX, REAL areaY,
                                      REAL areaZ) {
  std::size_t plnstrnType = static_cast<std::size_t>(
      dem::Parameter::getSingleton().parameter["plnstrnType"]);
  std::size_t unloadStep = static_cast<std::size_t>(
      dem::Parameter::getSingleton().parameter["unloadStep"]);

  //REAL forceDamp = dem::Parameter::getSingleton().parameter["forceDamp"];
  //REAL massScale = dem::Parameter::getSingleton().parameter["massScale"];
  //REAL mass = dem::Parameter::getSingleton().parameter["boundaryMass"];
  REAL boundaryRate = dem::Parameter::getSingleton().parameter["boundaryRate"];
  REAL sideRateRatio =
      dem::Parameter::getSingleton().parameter["sideRateRatio"];
  REAL tol = dem::Parameter::getSingleton().parameter["tractionErrorTol"];
  //REAL atf = forceDamp * 2;

  REAL vel, pos;
  switch (id) { // boundary x1(1) and boundary x2(2) do not move
  case 3:
    if (fabs(normal.getY() / areaY + sigma) / sigma > tol) {
      vel = ((normal.getY() + sigma * areaY) > 0 ? 1 : -1) * boundaryRate *
            sideRateRatio;
      //vel = prevVeloc.getY() * (2-atf) / (2+atf) + (normal.getY() + sigma *
      //areaY) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.getY() + vel * timeStep;
      setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ()));
      setPoint(Vec(getPoint().getX(), pos, getPoint().getZ()));
    }
    break;
  case 4:
    if (fabs(normal.getY() / areaY - sigma) / sigma > tol) {
      vel = ((normal.getY() - sigma * areaY) > 0 ? 1 : -1) * boundaryRate *
            sideRateRatio;
      //vel = prevVeloc.getY() * (2-atf) / (2+atf) + (normal.getY() - sigma *
      //areaY) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.getY() + vel * timeStep;
      setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ()));
      setPoint(Vec(getPoint().getX(), pos, getPoint().getZ()));
    }
    break;
  // displacement control, leading to zero volumetric strain
  case 5:
    if (plnstrnType == 1)
      vel = boundaryRate;
    else if (plnstrnType == 2) {
      if (iteration <= unloadStep) // loading
        vel = boundaryRate;
      else if (iteration > unloadStep &&
               fabs(normal.getZ() / areaZ) >= 1.5 * sigma &&
               iteration <= 1.5 * unloadStep) // unloading
        vel = -boundaryRate;
      else if (iteration >
               1.5 * unloadStep) // reloading. Note there are invalid loops if
                                 // stress < and iteration <, it is OK.
        vel = boundaryRate;
      else {
        vel = 0.0;
      }
    } else {
      vel = 0.0;
    }
    //vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() + sigma *
    //areaZ) / mass * timeStep * 2 / (2 + atf);
    pos = prevPoint.getZ() + vel * timeStep;
    setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel));
    setPoint(Vec(getPoint().getX(), getPoint().getY(), pos));
    break;
  case 6:
    if (plnstrnType == 1)
      vel = -boundaryRate;
    else if (plnstrnType == 2) {
      if (iteration <= unloadStep) // loading
        vel = -boundaryRate;
      else if (iteration > unloadStep &&
               fabs(normal.getZ() / areaZ) >= 1.5 * sigma &&
               iteration <= 1.5 * unloadStep) // unloading
        vel = boundaryRate;
      else if (iteration >
               1.5 * unloadStep) // reloading. Note there are invalid loops if
                                 // stress < and iteration <, it is OK.
        vel = -boundaryRate;
      else {
        vel = 0.0;
      }
    } else {
      vel = 0.0;
    }
    //vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() - sigma *
    //areaZ) / mass * timeStep * 2 / (2 + atf);
    pos = prevPoint.getZ() + vel * timeStep;
    setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel));
    setPoint(Vec(getPoint().getX(), getPoint().getY(), pos));
    break;
  }
  prevPoint = point;
  prevVeloc = veloc;
}

void planeBoundary::updateTrueTriaxial(REAL sigma, REAL areaX, REAL areaY,
                                       REAL areaZ, REAL sigmaX, REAL sigmaY) {
  // sigma implies sigmaZ

  //REAL forceDamp = dem::Parameter::getSingleton().parameter["forceDamp"];
  //REAL massScale = dem::Parameter::getSingleton().parameter["massScale"];
  //REAL mass = dem::Parameter::getSingleton().parameter["boundaryMass"];
  REAL boundaryRate = dem::Parameter::getSingleton().parameter["boundaryRate"];
  REAL tol = dem::Parameter::getSingleton().parameter["tractionErrorTol"];
  //REAL atf = forceDamp * 2;

  REAL vel, pos;
  switch (id) {
  case 1:
    if (fabs(normal.getX() / areaX + sigmaX) / sigmaX > tol) {
      vel = ((normal.getX() + sigmaX * areaX) > 0 ? 1 : -1) * boundaryRate;
      //vel = prevVeloc.getX() * (2-atf) / (2+atf) + (normal.getX() + sigma *
      //areaX) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.getX() + vel * timeStep;
      setVeloc(Vec(vel, getVeloc().getY(), getVeloc().getZ()));
      setPoint(Vec(pos, getPoint().getY(), getPoint().getZ()));
    }
    break;
  case 2:
    if (fabs(normal.getX() / areaX - sigmaX) / sigmaX > tol) {
      vel = ((normal.getX() - sigmaX * areaX) > 0 ? 1 : -1) * boundaryRate;
      //vel = prevVeloc.getX() * (2-atf) / (2+atf) + (normal.getX() - sigma *
      //areaX) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.getX() + vel * timeStep;
      setVeloc(Vec(vel, getVeloc().getY(), getVeloc().getZ()));
      setPoint(Vec(pos, getPoint().getY(), getPoint().getZ()));
    }
    break;
  case 3:
    if (fabs(normal.getY() / areaY + sigmaY) / sigmaY > tol) {
      vel = ((normal.getY() + sigmaY * areaY) > 0 ? 1 : -1) * boundaryRate;
      //vel = prevVeloc.getY() * (2-atf) / (2+atf) + (normal.getY() + sigma *
      //areaY) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.getY() + vel * timeStep;
      setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ()));
      setPoint(Vec(getPoint().getX(), pos, getPoint().getZ()));
    }
    break;
  case 4:
    if (fabs(normal.getY() / areaY - sigmaY) / sigmaY > tol) {
      vel = ((normal.getY() - sigmaY * areaY) > 0 ? 1 : -1) * boundaryRate;
      //vel = prevVeloc.getY() * (2-atf) / (2+atf) + (normal.getY() - sigma *
      //areaY) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.getY() + vel * timeStep;
      setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ()));
      setPoint(Vec(getPoint().getX(), pos, getPoint().getZ()));
    }
    break;
  case 5:
    if (fabs(normal.getZ() / areaZ + sigma) / sigma > tol) {
      vel = ((normal.getZ() + sigma * areaZ) > 0 ? 1 : -1) * boundaryRate;
      //vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() + sigma *
      //areaZ) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.getZ() + vel * timeStep;
      setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel));
      setPoint(Vec(getPoint().getX(), getPoint().getY(), pos));
    }
    break;
  case 6:
    if (fabs(normal.getZ() / areaZ - sigma) / sigma > tol) {
      vel = ((normal.getZ() - sigma * areaZ) > 0 ? 1 : -1) * boundaryRate;
      //vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() - sigma *
      //areaZ) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.getZ() + vel * timeStep;
      setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel));
      setPoint(Vec(getPoint().getX(), getPoint().getY(), pos));
    }
    break;
  }
  prevPoint = point;
  prevVeloc = veloc;
}

cylinderBoundary::cylinderBoundary(std::size_t tp, std::ifstream &ifs)
    : Boundary(tp, ifs) {
  REAL dx, dy, dz, px, py, pz;
  ifs >> dx >> dy >> dz >> px >> py >> pz >> radius;
  direc = Vec(dx, dy, dz);
  point = Vec(px, py, pz);
}

void cylinderBoundary::findBdryContact(ParticlePArray &ptcls) {
  possParticle.clear();
  contactInfo.clear();

  for (auto it = ptcls.begin(); it != ptcls.end(); ++it) {
    if ((*it)->getType() ==
        0) { // only process free particles, excluding type 5
      ;
    }
  }
}

void cylinderBoundary::boundaryForce(
    std::map<std::size_t, std::vector<BoundaryTgt> > &boundaryTgtMap) {
  // for each plane boundary, define a temparory variable vtmp to use,
  // better than define a member variable which needs to be cleared.
  // and vtmp is initialized as empty in each iteration.
  std::vector<BoundaryTgt> vtmp;

  // for each possible boundary particle
  for (auto it = possParticle.begin(); it != possParticle.end(); ++it)
    ; // (*it)->cylinderRBForce();

  // checkout tangential forces and displacements after each particle is
  // processed
  boundaryTgtMap[this->id] = vtmp;

  updateStatForce();
}

} // namespace dem ends
