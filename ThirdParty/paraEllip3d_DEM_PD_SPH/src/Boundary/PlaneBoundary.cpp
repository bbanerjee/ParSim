#include <Boundary/PlaneBoundary.h>
#include <Core/Util/Utility.h>
#include <DiscreteElements/DEMParticle.h>

// use both pointer to and variable of class DEMParticle

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
  //std::cout << "Plane boundary: " << dx << " " << dy << " " << dz << " " << px
  //          << " " << py << " " << pz << "\n";
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
    std::cerr
      << "**ERROR** Normal direction not found in plane boundary geometry\n";
    std::cerr << "  Add the <direction> [x, y, z] </direction> tag.";
    // **TODO** Throw exception
  }
  direc = Vec::fromString(vecStr);

  if (!ps["position"](vecStr)) {
    std::cerr
      << "**ERROR** Centroid poosition not found in plane boundary geometry\n";
    std::cerr << "  Add the <position> [x, y, z] </position> tag.";
    // **TODO** Throw exception
  }
  point = Vec::fromString(vecStr);
  //std::cout << "Plane boundary: " << direc.x() << " " << direc.y() << " "
  //          << direc.z() << " " << point.x() << " " << point.y() << " "
  //          << point.z() << "\n";
  veloc = 0;
  prevPoint = point;
  prevVeloc = veloc;

  // Read the planes defining the edges of the boundary plane
  b_extraNum = 0;
  for (auto edge_ps = ps["extraEdge"]; edge_ps; edge_ps.next()) {
    b_extraNum++;
    if (!edge_ps["direction"](vecStr)) {
      std::cerr << "**ERROR** Normal direction not found in edge geometry\n";
      std::cerr << "  Add the <direction> [x, y, z] </direction> tag.";
      // **TODO** Throw exception
    }
    Vec edgeDirec = Vec::fromString(vecStr);

    if (!edge_ps["position"](vecStr)) {
      std::cerr << "**ERROR** Centroid poosition not found in edge geometry\n";
      std::cerr << "  Add the <position> [x, y, z] </position> tag.";
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
    std::cerr
      << "**ERROR** Normal direction not found in plane boundary geometry\n";
    std::cerr << "  Add the direction: [x, y, z]  key-value pair.";
    // **TODO** Throw exception
  }
  direc = Vec::fromString(vecStr);

  try {
    vecStr = ps["position"].get<std::string>();
  } catch (std::exception) {
    std::cerr
      << "**ERROR** Centroid poosition not found in plane boundary geometry\n";
    std::cerr << "  Add the position: [x, y, z]  key-value pair.";
    // **TODO** Throw exception
  }
  point = Vec::fromString(vecStr);
  //std::cout << "Plane boundary: " << direc.x() << " " << direc.y() << " "
  //          << direc.z() << " " << point.x() << " " << point.y() << " "
  //          << point.z() << "\n";
  veloc = 0;
  prevPoint = point;
  prevVeloc = veloc;

  // Read the planes defining the edges of the boundary plane
  b_extraNum = 0;
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
      std::cerr << "**ERROR** Normal direction/Centroid position not found in "
                   "edge geometry\n";
      std::cerr << "  Add the direction: [x, y, z] tag.";
      std::cerr << "  Add the position: [x, y, z]  tag.";
    }
  }
}

void
PlaneBoundary::print(std::ostream& os)
{
  Boundary::print(os);
  os << std::setw(OWID) << direc.x() << std::setw(OWID) << direc.y()
     << std::setw(OWID) << direc.z() << std::setw(OWID) << point.x()
     << std::setw(OWID) << point.y() << std::setw(OWID) << point.z()
     << std::endl;

  for (auto& et : b_extraEdge)
    os << std::setw(OWID) << " " << std::setw(OWID) << et.getDirec().x()
       << std::setw(OWID) << et.getDirec().y() << std::setw(OWID)
       << et.getDirec().z() << std::setw(OWID) << et.getPoint().x()
       << std::setw(OWID) << et.getPoint().y() << std::setw(OWID)
       << et.getPoint().z() << std::endl;
}

void
PlaneBoundary::printContactInfo(std::ostream& os)
{
  Boundary::printContactInfo(os);
  os << std::setw(OWID) << " " << std::setw(OWID) << " " << std::setw(OWID)
     << " " << std::setw(OWID) << normal.x() << std::setw(OWID) << normal.y()
     << std::setw(OWID) << normal.z() << std::setw(OWID) << tangt.x()
     << std::setw(OWID) << tangt.y() << std::setw(OWID) << tangt.z()
     << std::setw(OWID) << penetr << std::endl
     << std::endl;
  ;
}

void
PlaneBoundary::findBdryContact(DEMParticlePArray& ptcls)
{
  possParticle.clear();
  contactInfo.clear();
  clearStatForce();

  for (auto& ptcl : ptcls) {
    if (ptcl->getType() == 0) { // only process free particles, excluding type 5
      REAL dist = distanceToBdry(ptcl->currentPosition());
      /*
      if (ptcl->getId() == 2 || ptcl->getId() == 95) {
        //std::cout << "Boundary distance: DEMParticle " << ptcl->getId() << ":"
        //          << std::setprecision(16) << dist << "\n";
      }
      */
      if (dist < 0 && fabs(dist) <= ptcl->getA()) {
        bool inside = true;
        for (auto& et : b_extraEdge) {
          REAL eDist = distanceToBdry(ptcl->currentPosition(), et);
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

  // REAL forceDamp = "forceDamp"];
  // REAL massScale = "massScale"];
  // REAL mass = "boundaryMass"];
  REAL boundaryRate = util::getParam<REAL>("boundaryRate");
  REAL topSpeedup = util::getParam<REAL>("topSpeedup");
  REAL tol = util::getParam<REAL>("tractionErrorTol");
  // REAL atf = forceDamp * 2;

  REAL vel, pos;
  switch (b_id) {
    case 1:
      if (fabs(normal.x() / areaX + sigma) / sigma > tol) {
        vel = ((normal.x() + sigma * areaX) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.x() * (2-atf) / (2+atf) + (normal.x() + sigma *
        // areaX) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.x() + vel * timeStep;
        setVeloc(Vec(vel, getVeloc().y(), getVeloc().z()));
        setPoint(Vec(pos, getPoint().y(), getPoint().z()));
      }
      break;
    case 2:
      if (fabs(normal.x() / areaX - sigma) / sigma > tol) {
        vel = ((normal.x() - sigma * areaX) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.x() * (2-atf) / (2+atf) + (normal.x() - sigma *
        // areaX) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.x() + vel * timeStep;
        setVeloc(Vec(vel, getVeloc().y(), getVeloc().z()));
        setPoint(Vec(pos, getPoint().y(), getPoint().z()));
      }
      break;
    case 3:
      if (fabs(normal.y() / areaY + sigma) / sigma > tol) {
        vel = ((normal.y() + sigma * areaY) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.y() * (2-atf) / (2+atf) + (normal.y() + sigma *
        // areaY) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.y() + vel * timeStep;
        setVeloc(Vec(getVeloc().x(), vel, getVeloc().z()));
        setPoint(Vec(getPoint().x(), pos, getPoint().z()));
      }
      break;
    case 4:
      if (fabs(normal.y() / areaY - sigma) / sigma > tol) {
        vel = ((normal.y() - sigma * areaY) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.y() * (2-atf) / (2+atf) + (normal.y() - sigma *
        // areaY) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.y() + vel * timeStep;
        setVeloc(Vec(getVeloc().x(), vel, getVeloc().z()));
        setPoint(Vec(getPoint().x(), pos, getPoint().z()));
      }
      break;
    case 5:
      if (fabs(normal.z() / areaZ + sigma) / sigma > tol) {
        vel = ((normal.z() + sigma * areaZ) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.z() * (2-atf) / (2+atf) + (normal.z() + sigma *
        // areaZ) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.z() + vel * timeStep;
        setVeloc(Vec(getVeloc().x(), getVeloc().y(), vel));
        setPoint(Vec(getPoint().x(), getPoint().y(), pos));
      }
      break;
    case 6:
      if (fabs(normal.z() / areaZ - sigma) / sigma > tol) {
        vel = ((normal.z() - sigma * areaZ) > 0 ? 1 : -1) * boundaryRate;
        if (normal.z() == 0)
          vel = -boundaryRate * topSpeedup;
        // vel = prevVeloc.z() * (2-atf) / (2+atf) + (normal.z() - sigma *
        // areaZ) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.z() + vel * timeStep;
        setVeloc(Vec(getVeloc().x(), getVeloc().y(), vel));
        setPoint(Vec(getPoint().x(), getPoint().y(), pos));
      }
      break;
  }
  prevPoint = point;
  prevVeloc = veloc;
}

void
PlaneBoundary::updateOdometer(REAL sigma, REAL areaX, REAL areaY, REAL areaZ)
{

  // REAL forceDamp = "forceDamp"];
  // REAL massScale = "massScale"];
  // REAL mass = "boundaryMass"];
  REAL boundaryRate = util::getParam<REAL>("boundaryRate");
  REAL tol = util::getParam<REAL>("tractionErrorTol");
  // REAL atf = forceDamp * 2;

  REAL vel, pos;
  switch (b_id) {
    case 5:
      if (fabs(normal.z() / areaZ + sigma) / sigma > tol) {
        vel = ((normal.z() + sigma * areaZ) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.z() * (2-atf) / (2+atf) + (normal.z() + sigma *
        // areaZ) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.z() + vel * timeStep;
        setVeloc(Vec(getVeloc().x(), getVeloc().y(), vel));
        setPoint(Vec(getPoint().x(), getPoint().y(), pos));
      }
      break;
    case 6:
      if (fabs(normal.z() / areaZ - sigma) / sigma > tol) {
        vel = ((normal.z() - sigma * areaZ) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.z() * (2-atf) / (2+atf) + (normal.z() - sigma *
        // areaZ) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.z() + vel * timeStep;
        setVeloc(Vec(getVeloc().x(), getVeloc().y(), vel));
        setPoint(Vec(getPoint().x(), getPoint().y(), pos));
      }
      break;
  }
  prevPoint = point;
  prevVeloc = veloc;
}

void
PlaneBoundary::updateTriaxial(REAL sigma, REAL areaX, REAL areaY, REAL areaZ)
{
  auto triaxialType = util::getParam<std::size_t>("triaxialType");
  auto unloadStep = util::getParam<std::size_t>("unloadStep");

  // REAL forceDamp = util::getParam<REAL>("forceDamp");
  // REAL massScale = util::getParam<REAL>("massScale");
  // REAL mass = util::getParam<REAL>("boundaryMass");
  REAL boundaryRate = util::getParam<REAL>("boundaryRate");
  REAL tol = util::getParam<REAL>("tractionErrorTol");
  // REAL atf = forceDamp * 2;

  REAL vel = 0.0, pos = 0.0;
  switch (b_id) {
    case 1:
      if (fabs(normal.x() / areaX + sigma) / sigma > tol) {
        vel = ((normal.x() + sigma * areaX) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.x() * (2-atf) / (2+atf) + (normal.x() + sigma *
        // areaX) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.x() + vel * timeStep;
        setVeloc(Vec(vel, getVeloc().y(), getVeloc().z()));
        setPoint(Vec(pos, getPoint().y(), getPoint().z()));
      }
      break;
    case 2:
      if (fabs(normal.x() / areaX - sigma) / sigma > tol) {
        vel = ((normal.x() - sigma * areaX) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.x() * (2-atf) / (2+atf) + (normal.x() - sigma *
        // areaX) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.x() + vel * timeStep;
        setVeloc(Vec(vel, getVeloc().y(), getVeloc().z()));
        setPoint(Vec(pos, getPoint().y(), getPoint().z()));
      }
      break;
    case 3:
      if (fabs(normal.y() / areaY + sigma) / sigma > tol) {
        vel = ((normal.y() + sigma * areaY) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.y() * (2-atf) / (2+atf) + (normal.y() + sigma *
        // areaY) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.y() + vel * timeStep;
        setVeloc(Vec(getVeloc().x(), vel, getVeloc().z()));
        setPoint(Vec(getPoint().x(), pos, getPoint().z()));
      }
      break;
    case 4:
      if (fabs(normal.y() / areaY - sigma) / sigma > tol) {
        vel = ((normal.y() - sigma * areaY) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.y() * (2-atf) / (2+atf) + (normal.y() - sigma *
        // areaY) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.y() + vel * timeStep;
        setVeloc(Vec(getVeloc().x(), vel, getVeloc().z()));
        setPoint(Vec(getPoint().x(), pos, getPoint().z()));
      }
      break;
    case 5:
      if (triaxialType == 1)
        vel = boundaryRate;
      else if (triaxialType == 2) {
        if (iteration <= unloadStep) // loading
          vel = boundaryRate;
        else if (iteration > unloadStep &&
                 fabs(normal.z() / areaZ) >= 1.5 * sigma &&
                 iteration <= 1.5 * unloadStep) // unloading
          vel = -boundaryRate;
        else if (iteration >
                 1.5 * unloadStep) // reloading. Note there are invalid loops if
                                   // stress < and iteration <, it is OK.
          vel = boundaryRate;
      } else {
        vel = boundaryRate;
      }
      // vel = prevVeloc.z() * (2-atf) / (2+atf) + (normal.z() + sigma *
      // areaZ) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.z() + vel * timeStep;
      setVeloc(Vec(getVeloc().x(), getVeloc().y(), vel));
      setPoint(Vec(getPoint().x(), getPoint().y(), pos));
      break;
    case 6:
      if (triaxialType == 1)
        vel = -boundaryRate;
      else if (triaxialType == 2) {
        if (iteration <= unloadStep) // loading
          vel = -boundaryRate;
        else if (iteration > unloadStep &&
                 fabs(normal.z() / areaZ) >= 1.5 * sigma &&
                 iteration <= 1.5 * unloadStep) // unloading
          vel = boundaryRate;
        else if (iteration >
                 1.5 * unloadStep) // reloading. Note there are invalid loops if
                                   // stress < and iteration <, it is OK.
          vel = -boundaryRate;
      } else {
        vel = boundaryRate;
      }
      // vel = prevVeloc.z() * (2-atf) / (2+atf) + (normal.z() - sigma *
      // areaZ) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.z() + vel * timeStep;
      setVeloc(Vec(getVeloc().x(), getVeloc().y(), vel));
      setPoint(Vec(getPoint().x(), getPoint().y(), pos));
      break;
  }
  prevPoint = point;
  prevVeloc = veloc;
}

void
PlaneBoundary::updatePlaneStrain(REAL sigma, REAL areaX, REAL areaY, REAL areaZ)
{
  auto plnstrnType = util::getParam<std::size_t>("plnstrnType");
  auto unloadStep = util::getParam<std::size_t>("unloadStep");

  // REAL forceDamp = util::getParam<REAL>("forceDamp");
  // REAL massScale = util::getParam<REAL>("massScale");
  // REAL mass = util::getParam<REAL>("boundaryMass");
  REAL boundaryRate = util::getParam<REAL>("boundaryRate");
  REAL sideRateRatio = util::getParam<REAL>("sideRateRatio");
  REAL tol = util::getParam<REAL>("tractionErrorTol");
  // REAL atf = forceDamp * 2;

  REAL vel, pos;
  switch (b_id) { // boundary x1(1) and boundary x2(2) do not move
    case 3:
      if (fabs(normal.y() / areaY + sigma) / sigma > tol) {
        vel = ((normal.y() + sigma * areaY) > 0 ? 1 : -1) * boundaryRate *
              sideRateRatio;
        // vel = prevVeloc.y() * (2-atf) / (2+atf) + (normal.y() + sigma *
        // areaY) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.y() + vel * timeStep;
        setVeloc(Vec(getVeloc().x(), vel, getVeloc().z()));
        setPoint(Vec(getPoint().x(), pos, getPoint().z()));
      }
      break;
    case 4:
      if (fabs(normal.y() / areaY - sigma) / sigma > tol) {
        vel = ((normal.y() - sigma * areaY) > 0 ? 1 : -1) * boundaryRate *
              sideRateRatio;
        // vel = prevVeloc.y() * (2-atf) / (2+atf) + (normal.y() - sigma *
        // areaY) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.y() + vel * timeStep;
        setVeloc(Vec(getVeloc().x(), vel, getVeloc().z()));
        setPoint(Vec(getPoint().x(), pos, getPoint().z()));
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
                 fabs(normal.z() / areaZ) >= 1.5 * sigma &&
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
      // vel = prevVeloc.z() * (2-atf) / (2+atf) + (normal.z() + sigma *
      // areaZ) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.z() + vel * timeStep;
      setVeloc(Vec(getVeloc().x(), getVeloc().y(), vel));
      setPoint(Vec(getPoint().x(), getPoint().y(), pos));
      break;
    case 6:
      if (plnstrnType == 1)
        vel = -boundaryRate;
      else if (plnstrnType == 2) {
        if (iteration <= unloadStep) // loading
          vel = -boundaryRate;
        else if (iteration > unloadStep &&
                 fabs(normal.z() / areaZ) >= 1.5 * sigma &&
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
      // vel = prevVeloc.z() * (2-atf) / (2+atf) + (normal.z() - sigma *
      // areaZ) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.z() + vel * timeStep;
      setVeloc(Vec(getVeloc().x(), getVeloc().y(), vel));
      setPoint(Vec(getPoint().x(), getPoint().y(), pos));
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

  // REAL forceDamp = "forceDamp"];
  // REAL massScale = "massScale"];
  // REAL mass = "boundaryMass"];
  REAL boundaryRate = util::getParam<REAL>("boundaryRate");
  REAL tol = util::getParam<REAL>("tractionErrorTol");
  // REAL atf = forceDamp * 2;

  REAL vel, pos;
  switch (b_id) {
    case 1:
      if (fabs(normal.x() / areaX + sigmaX) / sigmaX > tol) {
        vel = ((normal.x() + sigmaX * areaX) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.x() * (2-atf) / (2+atf) + (normal.x() + sigma *
        // areaX) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.x() + vel * timeStep;
        setVeloc(Vec(vel, getVeloc().y(), getVeloc().z()));
        setPoint(Vec(pos, getPoint().y(), getPoint().z()));
      }
      break;
    case 2:
      if (fabs(normal.x() / areaX - sigmaX) / sigmaX > tol) {
        vel = ((normal.x() - sigmaX * areaX) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.x() * (2-atf) / (2+atf) + (normal.x() - sigma *
        // areaX) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.x() + vel * timeStep;
        setVeloc(Vec(vel, getVeloc().y(), getVeloc().z()));
        setPoint(Vec(pos, getPoint().y(), getPoint().z()));
      }
      break;
    case 3:
      if (fabs(normal.y() / areaY + sigmaY) / sigmaY > tol) {
        vel = ((normal.y() + sigmaY * areaY) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.y() * (2-atf) / (2+atf) + (normal.y() + sigma *
        // areaY) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.y() + vel * timeStep;
        setVeloc(Vec(getVeloc().x(), vel, getVeloc().z()));
        setPoint(Vec(getPoint().x(), pos, getPoint().z()));
      }
      break;
    case 4:
      if (fabs(normal.y() / areaY - sigmaY) / sigmaY > tol) {
        vel = ((normal.y() - sigmaY * areaY) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.y() * (2-atf) / (2+atf) + (normal.y() - sigma *
        // areaY) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.y() + vel * timeStep;
        setVeloc(Vec(getVeloc().x(), vel, getVeloc().z()));
        setPoint(Vec(getPoint().x(), pos, getPoint().z()));
      }
      break;
    case 5:
      if (fabs(normal.z() / areaZ + sigma) / sigma > tol) {
        vel = ((normal.z() + sigma * areaZ) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.z() * (2-atf) / (2+atf) + (normal.z() + sigma *
        // areaZ) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.z() + vel * timeStep;
        setVeloc(Vec(getVeloc().x(), getVeloc().y(), vel));
        setPoint(Vec(getPoint().x(), getPoint().y(), pos));
      }
      break;
    case 6:
      if (fabs(normal.z() / areaZ - sigma) / sigma > tol) {
        vel = ((normal.z() - sigma * areaZ) > 0 ? 1 : -1) * boundaryRate;
        // vel = prevVeloc.z() * (2-atf) / (2+atf) + (normal.z() - sigma *
        // areaZ) / mass * timeStep * 2 / (2 + atf);
        pos = prevPoint.z() + vel * timeStep;
        setVeloc(Vec(getVeloc().x(), getVeloc().y(), vel));
        setPoint(Vec(getPoint().x(), getPoint().y(), pos));
      }
      break;
  }
  prevPoint = point;
  prevVeloc = veloc;
}
