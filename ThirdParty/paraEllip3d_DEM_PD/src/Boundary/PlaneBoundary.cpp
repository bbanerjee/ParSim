#include <Boundary/PlaneBoundary.h>
#include <DiscreteElements/Particle.h>
// use both pointer to and variable of class Particle

using namespace dem;

PlaneBoundary::PlaneBoundary(BoundaryType tp, std::ifstream& ifs)
  : Boundary()
{
  // These are defined in the base class
  b_type = tp;
  ifs >> b_extraNum;
  ifs >> b_id;

  // These variables are local to the PlaneBoundary class
  REAL dx, dy, dz, px, py, pz;
  ifs >> dx >> dy >> dz >> px >> py >> pz;
  std::cout << "Plane boundary: " << dx << " " << dy << " " << dz << " " << px
            << " " << py << " " << pz << "\n";
  direc = Vec(dx, dy, dz);
  point = Vec(px, py, pz);
  prevPoint = point;
  prevVeloc = veloc = 0;
  for (std::size_t i = 0; i < b_extraNum; ++i) {
    ifs >> dx >> dy >> dz >> px >> py >> pz;
    b_extraEdge.push_back(Plane(Vec(dx, dy, dz), Vec(px, py, pz)));
  }
}

PlaneBoundary::PlaneBoundary(BoundaryId id, BoundaryType tp,
                             const XMLProblemSpec& ps)
  : Boundary()
{
  b_id = id;
  b_type = tp;

  std::string vecStr;
  if (!ps["direction"](vecStr)) {
    std::cout
      << "**ERROR** Normal direction not found in plane boundary geometry\n";
    std::cout << "  Add the <direction> [x, y, z] </direction> tag.";
    // **TODO** Throw exception
  }
  direc = Vec::fromString(vecStr);

  if (!ps["position"](vecStr)) {
    std::cout
      << "**ERROR** Centroid poosition not found in plane boundary geometry\n";
    std::cout << "  Add the <position> [x, y, z] </position> tag.";
    // **TODO** Throw exception
  }
  point = Vec::fromString(vecStr);
  std::cout << "Plane boundary: " << direc.getX() << " " << direc.getY() << " "
            << direc.getZ() << " " << point.getX() << " " << point.getY() << " "
            << point.getZ() << "\n";
  veloc = 0;
  prevPoint = point;
  prevVeloc = veloc;

  // Read the planes defining the edges of the boundary plane
  for (auto edge_ps = ps["extraEdge"]; edge_ps; edge_ps.next()) {
    b_extraNum++;
    if (!edge_ps["direction"](vecStr)) {
      std::cout << "**ERROR** Normal direction not found in edge geometry\n";
      std::cout << "  Add the <direction> [x, y, z] </direction> tag.";
      // **TODO** Throw exception
    }
    Vec edgeDirec = Vec::fromString(vecStr);

    if (!edge_ps["position"](vecStr)) {
      std::cout << "**ERROR** Centroid poosition not found in edge geometry\n";
      std::cout << "  Add the <position> [x, y, z] </position> tag.";
      // **TODO** Throw exception
    }
    Vec edgePoint = Vec::fromString(vecStr);
    b_extraEdge.push_back(Plane(edgeDirec, edgePoint));
  }
}

PlaneBoundary::PlaneBoundary(BoundaryId id, BoundaryType tp,
                             const JsonProblemSpec& ps)
  : Boundary()
{
  b_id = id;
  b_type = tp;

  std::string vecStr;
  try {
    vecStr = ps["direction"].get<std::string>();
  } catch (std::exception) {
    std::cout
      << "**ERROR** Normal direction not found in plane boundary geometry\n";
    std::cout << "  Add the direction: [x, y, z]  key-value pair.";
    // **TODO** Throw exception
  }
  direc = Vec::fromString(vecStr);

  try {
    vecStr = ps["position"].get<std::string>();
  } catch (std::exception) {
    std::cout
      << "**ERROR** Centroid poosition not found in plane boundary geometry\n";
    std::cout << "  Add the position: [x, y, z]  key-value pair.";
    // **TODO** Throw exception
  }
  point = Vec::fromString(vecStr);
  std::cout << "Plane boundary: " << direc.getX() << " " << direc.getY() << " "
            << direc.getZ() << " " << point.getX() << " " << point.getY() << " "
            << point.getZ() << "\n";
  veloc = 0;
  prevPoint = point;
  prevVeloc = veloc;

  // Read the planes defining the edges of the boundary plane
  auto extraEdgeIter = ps.find("extraEdge");
  if (extraEdgeIter != ps.end()) {
    try {
      auto edge_ps = ps["extraEdge"];
      for (auto edge : edge_ps) {
        b_extraNum++;
        vecStr = edge["direction"].get<std::string>();
        vecStr = edge["position"].get<std::string>();
        Vec edgeDirec = Vec::fromString(vecStr);
        Vec edgePoint = Vec::fromString(vecStr);
        b_extraEdge.push_back(Plane(edgeDirec, edgePoint));
      }
    } catch (std::exception& e) {
      std::cout << "**ERROR** Normal direction/Centroid position not found in "
                   "edge geometry\n";
      std::cout << "  Add the direction: [x, y, z] tag.";
      std::cout << "  Add the position: [x, y, z]  tag.";
    }
  }
}

void
PlaneBoundary::print(std::ostream& os)
{
  Boundary::print(os);
  os << std::setw(OWID) << direc.getX() << std::setw(OWID) << direc.getY()
     << std::setw(OWID) << direc.getZ() << std::setw(OWID) << point.getX()
     << std::setw(OWID) << point.getY() << std::setw(OWID) << point.getZ()
     << std::endl;

  for (auto& et : b_extraEdge)
    os << std::setw(OWID) << " " << std::setw(OWID) << et.getDirec().getX()
       << std::setw(OWID) << et.getDirec().getY() << std::setw(OWID)
       << et.getDirec().getZ() << std::setw(OWID) << et.getPoint().getX()
       << std::setw(OWID) << et.getPoint().getY() << std::setw(OWID)
       << et.getPoint().getZ() << std::endl;
}

void
PlaneBoundary::printContactInfo(std::ostream& os)
{
  Boundary::printContactInfo(os);
  os << std::setw(OWID) << " " << std::setw(OWID) << " " << std::setw(OWID)
     << " " << std::setw(OWID) << normal.getX() << std::setw(OWID)
     << normal.getY() << std::setw(OWID) << normal.getZ() << std::setw(OWID)
     << tangt.getX() << std::setw(OWID) << tangt.getY() << std::setw(OWID)
     << tangt.getZ() << std::setw(OWID) << penetr << std::endl
     << std::endl;
  ;
}

void
PlaneBoundary::findBdryContact(ParticlePArray& ptcls)
{
  possParticle.clear();
  contactInfo.clear();
  clearStatForce();

  for (auto& ptcl : ptcls) {
    if (ptcl->getType() == 0) { // only process free particles, excluding type 5
      REAL dist = distanceToBdry(ptcl->getCurrPos());
      if (dist < 0 && fabs(dist) <= ptcl->getA()) {
        bool inside = true;
        for (auto& et : b_extraEdge) {
          REAL eDist = distanceToBdry(ptcl->getCurrPos(), et);
          if (eDist >= 0) {
            inside = false;
            break;
          }
        }
        if (inside)
          possParticle.push_back(ptcl);
      }
    }
  }
}

void
PlaneBoundary::boundaryForce(BoundaryTangentArrayMap& boundaryTgtMap)
{
  // for each plane boundary, define a temparory variable vtmp to use,
  // better than define a member variable which needs to be cleared.
  // and vtmp is initialized as empty in each iteration.
  BoundaryTangentArray vtmp;

  // for each possible boundary particle
  for (auto& it : possParticle)
    it->planeRBForce(this, boundaryTgtMap, vtmp);

  // checkout tangential forces and displacements after each particle is
  // processed
  boundaryTgtMap[this->b_id] = vtmp;

  updateStatForce();
}

void
PlaneBoundary::updateIsotropic(REAL sigma, REAL areaX, REAL areaY, REAL areaZ)
{

  // REAL forceDamp = dem::Parameter::getSingleton().parameter["forceDamp"];
  // REAL massScale = dem::Parameter::getSingleton().parameter["massScale"];
  // REAL mass = dem::Parameter::getSingleton().parameter["boundaryMass"];
  REAL boundaryRate = dem::Parameter::getSingleton().parameter["boundaryRate"];
  REAL topSpeedup = dem::Parameter::getSingleton().parameter["topSpeedup"];
  REAL tol = dem::Parameter::getSingleton().parameter["tractionErrorTol"];
  // REAL atf = forceDamp * 2;

  REAL vel, pos;
  switch (b_id) {
    case 1:
      if (fabs(normal.getX() / areaX + sigma) / sigma > tol) {
        vel = ((normal.getX() + sigma * areaX) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.getX() * (2-atf) / (2+atf) + (normal.getX() + sigma *
        // areaX) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.getX() + vel * timeStep;
        setVeloc(Vec(vel, getVeloc().getY(), getVeloc().getZ()));
        setPoint(Vec(pos, getPoint().getY(), getPoint().getZ()));
      }
      break;
    case 2:
      if (fabs(normal.getX() / areaX - sigma) / sigma > tol) {
        vel = ((normal.getX() - sigma * areaX) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.getX() * (2-atf) / (2+atf) + (normal.getX() - sigma *
        // areaX) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.getX() + vel * timeStep;
        setVeloc(Vec(vel, getVeloc().getY(), getVeloc().getZ()));
        setPoint(Vec(pos, getPoint().getY(), getPoint().getZ()));
      }
      break;
    case 3:
      if (fabs(normal.getY() / areaY + sigma) / sigma > tol) {
        vel = ((normal.getY() + sigma * areaY) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.getY() * (2-atf) / (2+atf) + (normal.getY() + sigma *
        // areaY) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.getY() + vel * timeStep;
        setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ()));
        setPoint(Vec(getPoint().getX(), pos, getPoint().getZ()));
      }
      break;
    case 4:
      if (fabs(normal.getY() / areaY - sigma) / sigma > tol) {
        vel = ((normal.getY() - sigma * areaY) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.getY() * (2-atf) / (2+atf) + (normal.getY() - sigma *
        // areaY) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.getY() + vel * timeStep;
        setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ()));
        setPoint(Vec(getPoint().getX(), pos, getPoint().getZ()));
      }
      break;
    case 5:
      if (fabs(normal.getZ() / areaZ + sigma) / sigma > tol) {
        vel = ((normal.getZ() + sigma * areaZ) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() + sigma *
        // areaZ) / mass * timeStep * 2 / (2 + atf);
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
        // vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() - sigma *
        // areaZ) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.getZ() + vel * timeStep;
        setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel));
        setPoint(Vec(getPoint().getX(), getPoint().getY(), pos));
      }
      break;
  }
  prevPoint = point;
  prevVeloc = veloc;
}

void
PlaneBoundary::updateOdometer(REAL sigma, REAL areaX, REAL areaY, REAL areaZ)
{

  // REAL forceDamp = dem::Parameter::getSingleton().parameter["forceDamp"];
  // REAL massScale = dem::Parameter::getSingleton().parameter["massScale"];
  // REAL mass = dem::Parameter::getSingleton().parameter["boundaryMass"];
  REAL boundaryRate = dem::Parameter::getSingleton().parameter["boundaryRate"];
  REAL tol = dem::Parameter::getSingleton().parameter["tractionErrorTol"];
  // REAL atf = forceDamp * 2;

  REAL vel, pos;
  switch (b_id) {
    case 5:
      if (fabs(normal.getZ() / areaZ + sigma) / sigma > tol) {
        vel = ((normal.getZ() + sigma * areaZ) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() + sigma *
        // areaZ) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.getZ() + vel * timeStep;
        setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel));
        setPoint(Vec(getPoint().getX(), getPoint().getY(), pos));
      }
      break;
    case 6:
      if (fabs(normal.getZ() / areaZ - sigma) / sigma > tol) {
        vel = ((normal.getZ() - sigma * areaZ) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() - sigma *
        // areaZ) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.getZ() + vel * timeStep;
        setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel));
        setPoint(Vec(getPoint().getX(), getPoint().getY(), pos));
      }
      break;
  }
  prevPoint = point;
  prevVeloc = veloc;
}

void
PlaneBoundary::updateTriaxial(REAL sigma, REAL areaX, REAL areaY, REAL areaZ)
{
  std::size_t triaxialType = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["triaxialType"]);
  std::size_t unloadStep = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["unloadStep"]);

  // REAL forceDamp = dem::Parameter::getSingleton().parameter["forceDamp"];
  // REAL massScale = dem::Parameter::getSingleton().parameter["massScale"];
  // REAL mass = dem::Parameter::getSingleton().parameter["boundaryMass"];
  REAL boundaryRate = dem::Parameter::getSingleton().parameter["boundaryRate"];
  REAL tol = dem::Parameter::getSingleton().parameter["tractionErrorTol"];
  // REAL atf = forceDamp * 2;

  REAL vel = 0.0, pos = 0.0;
  switch (b_id) {
    case 1:
      if (fabs(normal.getX() / areaX + sigma) / sigma > tol) {
        vel = ((normal.getX() + sigma * areaX) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.getX() * (2-atf) / (2+atf) + (normal.getX() + sigma *
        // areaX) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.getX() + vel * timeStep;
        setVeloc(Vec(vel, getVeloc().getY(), getVeloc().getZ()));
        setPoint(Vec(pos, getPoint().getY(), getPoint().getZ()));
      }
      break;
    case 2:
      if (fabs(normal.getX() / areaX - sigma) / sigma > tol) {
        vel = ((normal.getX() - sigma * areaX) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.getX() * (2-atf) / (2+atf) + (normal.getX() - sigma *
        // areaX) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.getX() + vel * timeStep;
        setVeloc(Vec(vel, getVeloc().getY(), getVeloc().getZ()));
        setPoint(Vec(pos, getPoint().getY(), getPoint().getZ()));
      }
      break;
    case 3:
      if (fabs(normal.getY() / areaY + sigma) / sigma > tol) {
        vel = ((normal.getY() + sigma * areaY) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.getY() * (2-atf) / (2+atf) + (normal.getY() + sigma *
        // areaY) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.getY() + vel * timeStep;
        setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ()));
        setPoint(Vec(getPoint().getX(), pos, getPoint().getZ()));
      }
      break;
    case 4:
      if (fabs(normal.getY() / areaY - sigma) / sigma > tol) {
        vel = ((normal.getY() - sigma * areaY) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.getY() * (2-atf) / (2+atf) + (normal.getY() - sigma *
        // areaY) / mass * timeStep * 2 / (2 + atf);
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
      // vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() + sigma *
      // areaZ) / mass * timeStep * 2 / (2 + atf);
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
      // vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() - sigma *
      // areaZ) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.getZ() + vel * timeStep;
      setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel));
      setPoint(Vec(getPoint().getX(), getPoint().getY(), pos));
      break;
  }
  prevPoint = point;
  prevVeloc = veloc;
}

void
PlaneBoundary::updatePlaneStrain(REAL sigma, REAL areaX, REAL areaY, REAL areaZ)
{
  std::size_t plnstrnType = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["plnstrnType"]);
  std::size_t unloadStep = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["unloadStep"]);

  // REAL forceDamp = dem::Parameter::getSingleton().parameter["forceDamp"];
  // REAL massScale = dem::Parameter::getSingleton().parameter["massScale"];
  // REAL mass = dem::Parameter::getSingleton().parameter["boundaryMass"];
  REAL boundaryRate = dem::Parameter::getSingleton().parameter["boundaryRate"];
  REAL sideRateRatio =
    dem::Parameter::getSingleton().parameter["sideRateRatio"];
  REAL tol = dem::Parameter::getSingleton().parameter["tractionErrorTol"];
  // REAL atf = forceDamp * 2;

  REAL vel, pos;
  switch (b_id) { // boundary x1(1) and boundary x2(2) do not move
    case 3:
      if (fabs(normal.getY() / areaY + sigma) / sigma > tol) {
        vel = ((normal.getY() + sigma * areaY) > 0 ? 1 : -1) * boundaryRate *
              sideRateRatio;
        // vel = prevVeloc.getY() * (2-atf) / (2+atf) + (normal.getY() + sigma *
        // areaY) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.getY() + vel * timeStep;
        setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ()));
        setPoint(Vec(getPoint().getX(), pos, getPoint().getZ()));
      }
      break;
    case 4:
      if (fabs(normal.getY() / areaY - sigma) / sigma > tol) {
        vel = ((normal.getY() - sigma * areaY) > 0 ? 1 : -1) * boundaryRate *
              sideRateRatio;
        // vel = prevVeloc.getY() * (2-atf) / (2+atf) + (normal.getY() - sigma *
        // areaY) / mass * timeStep * 2 / (2 + atf);
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
      // vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() + sigma *
      // areaZ) / mass * timeStep * 2 / (2 + atf);
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
      // vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() - sigma *
      // areaZ) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.getZ() + vel * timeStep;
      setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel));
      setPoint(Vec(getPoint().getX(), getPoint().getY(), pos));
      break;
  }
  prevPoint = point;
  prevVeloc = veloc;
}

void
PlaneBoundary::updateTrueTriaxial(REAL sigma, REAL areaX, REAL areaY,
                                  REAL areaZ, REAL sigmaX, REAL sigmaY)
{
  // sigma implies sigmaZ

  // REAL forceDamp = dem::Parameter::getSingleton().parameter["forceDamp"];
  // REAL massScale = dem::Parameter::getSingleton().parameter["massScale"];
  // REAL mass = dem::Parameter::getSingleton().parameter["boundaryMass"];
  REAL boundaryRate = dem::Parameter::getSingleton().parameter["boundaryRate"];
  REAL tol = dem::Parameter::getSingleton().parameter["tractionErrorTol"];
  // REAL atf = forceDamp * 2;

  REAL vel, pos;
  switch (b_id) {
    case 1:
      if (fabs(normal.getX() / areaX + sigmaX) / sigmaX > tol) {
        vel = ((normal.getX() + sigmaX * areaX) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.getX() * (2-atf) / (2+atf) + (normal.getX() + sigma *
        // areaX) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.getX() + vel * timeStep;
        setVeloc(Vec(vel, getVeloc().getY(), getVeloc().getZ()));
        setPoint(Vec(pos, getPoint().getY(), getPoint().getZ()));
      }
      break;
    case 2:
      if (fabs(normal.getX() / areaX - sigmaX) / sigmaX > tol) {
        vel = ((normal.getX() - sigmaX * areaX) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.getX() * (2-atf) / (2+atf) + (normal.getX() - sigma *
        // areaX) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.getX() + vel * timeStep;
        setVeloc(Vec(vel, getVeloc().getY(), getVeloc().getZ()));
        setPoint(Vec(pos, getPoint().getY(), getPoint().getZ()));
      }
      break;
    case 3:
      if (fabs(normal.getY() / areaY + sigmaY) / sigmaY > tol) {
        vel = ((normal.getY() + sigmaY * areaY) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.getY() * (2-atf) / (2+atf) + (normal.getY() + sigma *
        // areaY) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.getY() + vel * timeStep;
        setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ()));
        setPoint(Vec(getPoint().getX(), pos, getPoint().getZ()));
      }
      break;
    case 4:
      if (fabs(normal.getY() / areaY - sigmaY) / sigmaY > tol) {
        vel = ((normal.getY() - sigmaY * areaY) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.getY() * (2-atf) / (2+atf) + (normal.getY() - sigma *
        // areaY) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.getY() + vel * timeStep;
        setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ()));
        setPoint(Vec(getPoint().getX(), pos, getPoint().getZ()));
      }
      break;
    case 5:
      if (fabs(normal.getZ() / areaZ + sigma) / sigma > tol) {
        vel = ((normal.getZ() + sigma * areaZ) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() + sigma *
        // areaZ) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.getZ() + vel * timeStep;
        setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel));
        setPoint(Vec(getPoint().getX(), getPoint().getY(), pos));
      }
      break;
    case 6:
      if (fabs(normal.getZ() / areaZ - sigma) / sigma > tol) {
        vel = ((normal.getZ() - sigma * areaZ) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() - sigma *
        // areaZ) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.getZ() + vel * timeStep;
        setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel));
        setPoint(Vec(getPoint().getX(), getPoint().getY(), pos));
      }
      break;
  }
  prevPoint = point;
  prevVeloc = veloc;
}
