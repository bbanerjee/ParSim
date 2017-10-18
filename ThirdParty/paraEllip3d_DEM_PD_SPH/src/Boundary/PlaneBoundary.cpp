#include <Boundary/PlaneBoundary.h>
#include <Core/Util/Utility.h>
#include <DiscreteElements/DEMParticle.h>

// use both pointer to and variable of class DEMParticle

using namespace dem;

PlaneBoundary::PlaneBoundary(Boundary::BoundaryType tp, std::istream& ifs)
  : Boundary()
{
  // These are defined in the base class
  b_type = tp;
  ifs >> b_extraNum;
  int id;
  ifs >> id;
  b_id = Boundary::getBoundaryID(id);

  // These variables are local to the PlaneBoundary class
  REAL dx, dy, dz, px, py, pz;
  ifs >> dx >> dy >> dz >> px >> py >> pz;
  //std::cout << "Plane boundary: " << dx << " " << dy << " " << dz << " " << px
  //          << " " << py << " " << pz << "\n";
  d_direction = Vec(dx, dy, dz);
  d_position = Vec(px, py, pz);
  d_previousPosition = d_position;
  d_previousVelocity = d_velocity = 0;
  for (std::size_t i = 0; i < b_extraNum; ++i) {
    ifs >> dx >> dy >> dz >> px >> py >> pz;
    b_extraEdge.push_back(Plane(Vec(dx, dy, dz), Vec(px, py, pz)));
  }

  // Boundary conditions
  int bcType = 0;
  ifs >> bcType;
  d_bcType = static_cast<BCType>(bcType);
  if (d_bcType == BCType::TRACTION) {
    d_tractionBC.read(ifs);
  } else if (d_bcType == BCType::DISPLACEMENT) {
    d_displacementBC.read(ifs);
  } else {
    d_bcType = BCType::NONE; 
  }

}

PlaneBoundary::PlaneBoundary(Boundary::BoundaryType tp, BoundaryID id, 
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
  d_direction = Vec::fromString(vecStr);

  if (!ps["position"](vecStr)) {
    std::cerr
      << "**ERROR** Centroid position not found in plane boundary geometry\n";
    std::cerr << "  Add the <position> [x, y, z] </position> tag.";
    // **TODO** Throw exception
  }
  d_position = Vec::fromString(vecStr);
  d_previousPosition = d_position;

  d_velocity = 0;
  if (ps["initial_velocity"](vecStr)) {
    d_velocity = Vec::fromString(vecStr);
  } else {
    std::cerr
      << "**WARNING** Initial velocity not found for plane boundary at "
      << d_position << " with normal " << d_direction << "\n"
      << "\t\t Using initial velocity = 0. Add the <initial_velocity> [x, y, z] </initial_velocity> tag.\n";
  }
  d_previousVelocity = d_velocity;
  //std::cout << "Plane boundary: " 
  //  << d_direction.x() << " " << d_direction.y() << " " << d_direction.z() 
  //  << " " 
  //  << d_position.x() << " " << d_position.y() << " " << d_position.z() 
  //  << " " 
  //  << d_velocity.x() << " " << d_velocity.y() << " " << d_velocity.z() 
  // << "\n";

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

  // Read the boundary conditions (if any)
  d_bcType = BCType::NONE; 
  auto traction_bc_ps = ps["tractionBC"];
  auto disp_bc_ps = ps["displacementBC"];
  if (traction_bc_ps) {
    d_bcType = BCType::TRACTION; 
    d_tractionBC.read(traction_bc_ps);
  } else if (disp_bc_ps) {
    d_bcType = BCType::DISPLACEMENT; 
    d_displacementBC.read(disp_bc_ps);
  }
}

PlaneBoundary::PlaneBoundary(Boundary::BoundaryType tp, BoundaryID id, 
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
  d_direction = Vec::fromString(vecStr);

  try {
    vecStr = ps["position"].get<std::string>();
  } catch (std::exception) {
    std::cerr
      << "**ERROR** Centroid poosition not found in plane boundary geometry\n";
    std::cerr << "  Add the position: [x, y, z]  key-value pair.";
    // **TODO** Throw exception
  }
  d_position = Vec::fromString(vecStr);

  d_velocity = 0;
  try {
    vecStr = ps["initial_velocity"].get<std::string>();
    d_velocity = Vec::fromString(vecStr);
  } catch (std::exception) {
    std::cerr
      << "**WARNING** Centroid poosition not found in plane boundary geometry\n"
      << "  Using 0 initial velocity. Add the initial_velocity: [x, y, z]  key-value pair.";
    // **TODO** Throw exception
  }
  //std::cout << "Plane boundary: " << d_direction.x() << " " << d_direction.y() << " "
  //          << d_direction.z() << " " << d_position.x() << " " << d_position.y() << " "
  //          << d_position.z() << "\n";
  d_previousPosition = d_position;
  d_previousVelocity = d_velocity;

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

  // Read the boundary conditions (if any)
  d_bcType = BCType::NONE; 
  auto traction_bc_iter = ps.find("tractionBC");
  auto disp_bc_iter = ps.find("displacementBC");
  if (traction_bc_iter != ps.end()) {
    auto traction_bc_ps = ps["tractionBC"];
    d_bcType = BCType::TRACTION; 
    d_tractionBC.read(traction_bc_ps);
  } else if (disp_bc_iter != ps.end()) {
    auto disp_bc_ps = ps["displacementBC"];
    d_bcType = BCType::DISPLACEMENT; 
    d_displacementBC.read(disp_bc_ps);
  }
}

void
PlaneBoundary::print(std::ostream& os)
{
  Boundary::print(os);
  os << std::setw(OWID) << d_direction.x() << std::setw(OWID) << d_direction.y()
     << std::setw(OWID) << d_direction.z() << std::setw(OWID) << d_position.x()
     << std::setw(OWID) << d_position.y() << std::setw(OWID) << d_position.z()
     << std::endl;

  for (auto& et : b_extraEdge)
    os << std::setw(OWID) << " " << std::setw(OWID) << et.getDirection().x()
       << std::setw(OWID) << et.getDirection().y() << std::setw(OWID)
       << et.getDirection().z() << std::setw(OWID) << et.getPosition().x()
       << std::setw(OWID) << et.getPosition().y() << std::setw(OWID)
       << et.getPosition().z() << std::endl;
}

void
PlaneBoundary::printContactInfo(std::ostream& os)
{
  Boundary::printContactInfo(os);
  os << std::setw(OWID) << " " << std::setw(OWID) << " " << std::setw(OWID)
     << " " << std::setw(OWID) << b_normalForce.x() << std::setw(OWID) << b_normalForce.y()
     << std::setw(OWID) << b_normalForce.z() << std::setw(OWID) << b_tangentForce.x()
     << std::setw(OWID) << b_tangentForce.y() << std::setw(OWID) << b_tangentForce.z()
     << std::setw(OWID) << b_penetration << std::endl
     << std::endl;
  ;
}

void
PlaneBoundary::findBoundaryContacts(DEMParticlePArray& particles)
{
  b_probableBoundaryParticles.clear();
  b_contacts.clear();
  clearStatForce();

  for (auto& particle : particles) {
    // only process free particles, excluding type 5
    if (particle->getType() == DEMParticle::DEMParticleType::FREE) { 
      REAL dist = distanceToBdry(particle->currentPosition());
      /*
      if (particle->getId() == 2 || particle->getId() == 95) {
        //std::cout << "Boundary distance: DEMParticle " << particle->getId() << ":"
        //          << std::setprecision(16) << dist << "\n";
      }
      */
      if (dist < 0 && fabs(dist) <= particle->getA()) {
        bool inside = true;
        for (auto& et : b_extraEdge) {
          REAL eDist = distanceToBdry(particle->currentPosition(), et);
          if (eDist >= 0) {
            inside = false;
            break;
          }
        }
        if (inside)
          b_probableBoundaryParticles.push_back(particle);
      }
    }
  }
}

void
PlaneBoundary::boundaryForce(BoundaryTangentArrayMap& boundaryTangentMap)
{
  // for each plane boundary, define a temparory variable vtmp to use,
  // better than define a member variable which needs to be cleared.
  // and vtmp is initialized as empty in each iteration.
  BoundaryTangentArray vtmp;

  // for each possible boundary particle
  for (auto& it : b_probableBoundaryParticles)
    it->planeRBForce(this, boundaryTangentMap, vtmp);

  // checkout tangential forces and displacements after each particle is
  // processed
  boundaryTangentMap[static_cast<size_t>(this->b_id)] = vtmp;

  updateStatForce();
}

void 
PlaneBoundary::updatePositionAndVelocity(double currTime, double deltaT,
                                         double area, double mass)
{
  if (d_bcType == BCType::TRACTION) {

    double traction_val = d_tractionBC.getBCValue(currTime);
    Vec tractionForce = d_direction * traction_val;
    tractionForce.setX(tractionForce.x()*area);
    tractionForce.setY(tractionForce.y()*area);
    tractionForce.setZ(tractionForce.z()*area);
    updateUsingTraction(deltaT, mass, tractionForce);

  } else if (d_bcType == BCType::DISPLACEMENT) {

    double disp_val = d_displacementBC.getBCValue(currTime);
    double disp_val_prev = d_displacementBC.getBCValue(currTime-deltaT/2);
    double disp_val_next = d_displacementBC.getBCValue(currTime+deltaT/2);
    double disp_rate = (disp_val_next - disp_val_prev)/deltaT;
    Vec disp = d_direction * disp_val;
    Vec rate = d_direction * disp_rate;
    updateUsingDisplacement(deltaT, disp, rate);

  } else {

    d_previousPosition = d_position;
    d_position += (d_velocity*deltaT);
    d_previousVelocity = d_velocity;

  }
}

void 
PlaneBoundary::updateUsingTraction(double deltaT, double mass,
                                   const Vec& tractionForce) 
{
  Vec resultantForce = tractionForce - b_normalForce;

  d_previousVelocity = d_velocity;
  d_velocity += (resultantForce/mass)*deltaT;

  d_previousPosition = d_position;
  d_position += d_velocity*deltaT;
}

void 
PlaneBoundary::updateUsingDisplacement(double deltaT, 
                                       const Vec& disp, const Vec& dispRate)
{
  d_previousVelocity = d_velocity;
  d_velocity = dispRate;

  d_previousPosition = d_position;
  d_position += disp;
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
    case BoundaryID::XMINUS:
      if (fabs(b_normalForce.x() / areaX + sigma) / sigma > tol) {
        vel = ((b_normalForce.x() + sigma * areaX) > 0 ? 1 : -1) * boundaryRate;
        // vel = d_previousVelocity.x() * (2-atf) / (2+atf) + (b_normalForce.x() + sigma *
        // areaX) / mass * timeStep * 2 / (2 + atf);
        pos = d_previousPosition.x() + vel * timeStep;
        setVelocity(Vec(vel, getVelocity().y(), getVelocity().z()));
        setPosition(Vec(pos, getPosition().y(), getPosition().z()));
      }
      break;
    case BoundaryID::XPLUS:
      if (fabs(b_normalForce.x() / areaX - sigma) / sigma > tol) {
        vel = ((b_normalForce.x() - sigma * areaX) > 0 ? 1 : -1) * boundaryRate;
        // vel = d_previousVelocity.x() * (2-atf) / (2+atf) + (b_normalForce.x() - sigma *
        // areaX) / mass * timeStep * 2 / (2 + atf);
        pos = d_previousPosition.x() + vel * timeStep;
        setVelocity(Vec(vel, getVelocity().y(), getVelocity().z()));
        setPosition(Vec(pos, getPosition().y(), getPosition().z()));
      }
      break;
    case BoundaryID::YMINUS:
      if (fabs(b_normalForce.y() / areaY + sigma) / sigma > tol) {
        vel = ((b_normalForce.y() + sigma * areaY) > 0 ? 1 : -1) * boundaryRate;
        // vel = d_previousVelocity.y() * (2-atf) / (2+atf) + (b_normalForce.y() + sigma *
        // areaY) / mass * timeStep * 2 / (2 + atf);
        pos = d_previousPosition.y() + vel * timeStep;
        setVelocity(Vec(getVelocity().x(), vel, getVelocity().z()));
        setPosition(Vec(getPosition().x(), pos, getPosition().z()));
      }
      break;
    case BoundaryID::YPLUS:
      if (fabs(b_normalForce.y() / areaY - sigma) / sigma > tol) {
        vel = ((b_normalForce.y() - sigma * areaY) > 0 ? 1 : -1) * boundaryRate;
        // vel = d_previousVelocity.y() * (2-atf) / (2+atf) + (b_normalForce.y() - sigma *
        // areaY) / mass * timeStep * 2 / (2 + atf);
        pos = d_previousPosition.y() + vel * timeStep;
        setVelocity(Vec(getVelocity().x(), vel, getVelocity().z()));
        setPosition(Vec(getPosition().x(), pos, getPosition().z()));
      }
      break;
    case BoundaryID::ZMINUS:
      if (fabs(b_normalForce.z() / areaZ + sigma) / sigma > tol) {
        vel = ((b_normalForce.z() + sigma * areaZ) > 0 ? 1 : -1) * boundaryRate;
        // vel = d_previousVelocity.z() * (2-atf) / (2+atf) + (b_normalForce.z() + sigma *
        // areaZ) / mass * timeStep * 2 / (2 + atf);
        pos = d_previousPosition.z() + vel * timeStep;
        setVelocity(Vec(getVelocity().x(), getVelocity().y(), vel));
        setPosition(Vec(getPosition().x(), getPosition().y(), pos));
      }
      break;
    case BoundaryID::ZPLUS:
      if (fabs(b_normalForce.z() / areaZ - sigma) / sigma > tol) {
        vel = ((b_normalForce.z() - sigma * areaZ) > 0 ? 1 : -1) * boundaryRate;
        if (b_normalForce.z() == 0)
          vel = -boundaryRate * topSpeedup;
        // vel = d_previousVelocity.z() * (2-atf) / (2+atf) + (b_normalForce.z() - sigma *
        // areaZ) / mass * timeStep * 2 / (2 + atf);
        pos = d_previousPosition.z() + vel * timeStep;
        setVelocity(Vec(getVelocity().x(), getVelocity().y(), vel));
        setPosition(Vec(getPosition().x(), getPosition().y(), pos));
      }
      break;
    default:
      break;
  }
  d_previousPosition = d_position;
  d_previousVelocity = d_velocity;
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
    case BoundaryID::ZMINUS:
      if (fabs(b_normalForce.z() / areaZ + sigma) / sigma > tol) {
        vel = ((b_normalForce.z() + sigma * areaZ) > 0 ? 1 : -1) * boundaryRate;
        // vel = d_previousVelocity.z() * (2-atf) / (2+atf) + (b_normalForce.z() + sigma *
        // areaZ) / mass * timeStep * 2 / (2 + atf);
        pos = d_previousPosition.z() + vel * timeStep;
        setVelocity(Vec(getVelocity().x(), getVelocity().y(), vel));
        setPosition(Vec(getPosition().x(), getPosition().y(), pos));
      }
      break;
    case BoundaryID::ZPLUS:
      if (fabs(b_normalForce.z() / areaZ - sigma) / sigma > tol) {
        vel = ((b_normalForce.z() - sigma * areaZ) > 0 ? 1 : -1) * boundaryRate;
        // vel = d_previousVelocity.z() * (2-atf) / (2+atf) + (b_normalForce.z() - sigma *
        // areaZ) / mass * timeStep * 2 / (2 + atf);
        pos = d_previousPosition.z() + vel * timeStep;
        setVelocity(Vec(getVelocity().x(), getVelocity().y(), vel));
        setPosition(Vec(getPosition().x(), getPosition().y(), pos));
      }
      break;
    default:
      break;
  }
  d_previousPosition = d_position;
  d_previousVelocity = d_velocity;
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
    case BoundaryID::XMINUS:
      if (fabs(b_normalForce.x() / areaX + sigma) / sigma > tol) {
        vel = ((b_normalForce.x() + sigma * areaX) > 0 ? 1 : -1) * boundaryRate;
        // vel = d_previousVelocity.x() * (2-atf) / (2+atf) + (b_normalForce.x() + sigma *
        // areaX) / mass * timeStep * 2 / (2 + atf);
        pos = d_previousPosition.x() + vel * timeStep;
        setVelocity(Vec(vel, getVelocity().y(), getVelocity().z()));
        setPosition(Vec(pos, getPosition().y(), getPosition().z()));
      }
      break;
    case BoundaryID::XPLUS:
      if (fabs(b_normalForce.x() / areaX - sigma) / sigma > tol) {
        vel = ((b_normalForce.x() - sigma * areaX) > 0 ? 1 : -1) * boundaryRate;
        // vel = d_previousVelocity.x() * (2-atf) / (2+atf) + (b_normalForce.x() - sigma *
        // areaX) / mass * timeStep * 2 / (2 + atf);
        pos = d_previousPosition.x() + vel * timeStep;
        setVelocity(Vec(vel, getVelocity().y(), getVelocity().z()));
        setPosition(Vec(pos, getPosition().y(), getPosition().z()));
      }
      break;
    case BoundaryID::YMINUS:
      if (fabs(b_normalForce.y() / areaY + sigma) / sigma > tol) {
        vel = ((b_normalForce.y() + sigma * areaY) > 0 ? 1 : -1) * boundaryRate;
        // vel = d_previousVelocity.y() * (2-atf) / (2+atf) + (b_normalForce.y() + sigma *
        // areaY) / mass * timeStep * 2 / (2 + atf);
        pos = d_previousPosition.y() + vel * timeStep;
        setVelocity(Vec(getVelocity().x(), vel, getVelocity().z()));
        setPosition(Vec(getPosition().x(), pos, getPosition().z()));
      }
      break;
    case BoundaryID::YPLUS:
      if (fabs(b_normalForce.y() / areaY - sigma) / sigma > tol) {
        vel = ((b_normalForce.y() - sigma * areaY) > 0 ? 1 : -1) * boundaryRate;
        // vel = d_previousVelocity.y() * (2-atf) / (2+atf) + (b_normalForce.y() - sigma *
        // areaY) / mass * timeStep * 2 / (2 + atf);
        pos = d_previousPosition.y() + vel * timeStep;
        setVelocity(Vec(getVelocity().x(), vel, getVelocity().z()));
        setPosition(Vec(getPosition().x(), pos, getPosition().z()));
      }
      break;
    case BoundaryID::ZMINUS:
      if (triaxialType == 1)
        vel = boundaryRate;
      else if (triaxialType == 2) {
        if (iteration <= unloadStep) // loading
          vel = boundaryRate;
        else if (iteration > unloadStep &&
                 fabs(b_normalForce.z() / areaZ) >= 1.5 * sigma &&
                 iteration <= 1.5 * unloadStep) // unloading
          vel = -boundaryRate;
        else if (iteration >
                 1.5 * unloadStep) // reloading. Note there are invalid loops if
                                   // stress < and iteration <, it is OK.
          vel = boundaryRate;
      } else {
        vel = boundaryRate;
      }
      // vel = d_previousVelocity.z() * (2-atf) / (2+atf) + (b_normalForce.z() + sigma *
      // areaZ) / mass * timeStep * 2 / (2 + atf);
      pos = d_previousPosition.z() + vel * timeStep;
      setVelocity(Vec(getVelocity().x(), getVelocity().y(), vel));
      setPosition(Vec(getPosition().x(), getPosition().y(), pos));
      break;
    case BoundaryID::ZPLUS:
      if (triaxialType == 1)
        vel = -boundaryRate;
      else if (triaxialType == 2) {
        if (iteration <= unloadStep) // loading
          vel = -boundaryRate;
        else if (iteration > unloadStep &&
                 fabs(b_normalForce.z() / areaZ) >= 1.5 * sigma &&
                 iteration <= 1.5 * unloadStep) // unloading
          vel = boundaryRate;
        else if (iteration >
                 1.5 * unloadStep) // reloading. Note there are invalid loops if
                                   // stress < and iteration <, it is OK.
          vel = -boundaryRate;
      } else {
        vel = boundaryRate;
      }
      // vel = d_previousVelocity.z() * (2-atf) / (2+atf) + (b_normalForce.z() - sigma *
      // areaZ) / mass * timeStep * 2 / (2 + atf);
      pos = d_previousPosition.z() + vel * timeStep;
      setVelocity(Vec(getVelocity().x(), getVelocity().y(), vel));
      setPosition(Vec(getPosition().x(), getPosition().y(), pos));
      break;
    default:
      break;
  }
  d_previousPosition = d_position;
  d_previousVelocity = d_velocity;
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
    case BoundaryID::YMINUS:
      if (fabs(b_normalForce.y() / areaY + sigma) / sigma > tol) {
        vel = ((b_normalForce.y() + sigma * areaY) > 0 ? 1 : -1) * boundaryRate *
              sideRateRatio;
        // vel = d_previousVelocity.y() * (2-atf) / (2+atf) + (b_normalForce.y() + sigma *
        // areaY) / mass * timeStep * 2 / (2 + atf);
        pos = d_previousPosition.y() + vel * timeStep;
        setVelocity(Vec(getVelocity().x(), vel, getVelocity().z()));
        setPosition(Vec(getPosition().x(), pos, getPosition().z()));
      }
      break;
    case BoundaryID::YPLUS:
      if (fabs(b_normalForce.y() / areaY - sigma) / sigma > tol) {
        vel = ((b_normalForce.y() - sigma * areaY) > 0 ? 1 : -1) * boundaryRate *
              sideRateRatio;
        // vel = d_previousVelocity.y() * (2-atf) / (2+atf) + (b_normalForce.y() - sigma *
        // areaY) / mass * timeStep * 2 / (2 + atf);
        pos = d_previousPosition.y() + vel * timeStep;
        setVelocity(Vec(getVelocity().x(), vel, getVelocity().z()));
        setPosition(Vec(getPosition().x(), pos, getPosition().z()));
      }
      break;
    // displacement control, leading to zero volumetric strain
    case BoundaryID::ZMINUS:
      if (plnstrnType == 1)
        vel = boundaryRate;
      else if (plnstrnType == 2) {
        if (iteration <= unloadStep) // loading
          vel = boundaryRate;
        else if (iteration > unloadStep &&
                 fabs(b_normalForce.z() / areaZ) >= 1.5 * sigma &&
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
      // vel = d_previousVelocity.z() * (2-atf) / (2+atf) + (b_normalForce.z() + sigma *
      // areaZ) / mass * timeStep * 2 / (2 + atf);
      pos = d_previousPosition.z() + vel * timeStep;
      setVelocity(Vec(getVelocity().x(), getVelocity().y(), vel));
      setPosition(Vec(getPosition().x(), getPosition().y(), pos));
      break;
    case BoundaryID::ZPLUS:
      if (plnstrnType == 1)
        vel = -boundaryRate;
      else if (plnstrnType == 2) {
        if (iteration <= unloadStep) // loading
          vel = -boundaryRate;
        else if (iteration > unloadStep &&
                 fabs(b_normalForce.z() / areaZ) >= 1.5 * sigma &&
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
      // vel = d_previousVelocity.z() * (2-atf) / (2+atf) + (b_normalForce.z() - sigma *
      // areaZ) / mass * timeStep * 2 / (2 + atf);
      pos = d_previousPosition.z() + vel * timeStep;
      setVelocity(Vec(getVelocity().x(), getVelocity().y(), vel));
      setPosition(Vec(getPosition().x(), getPosition().y(), pos));
      break;
    default:
      break;
  }
  d_previousPosition = d_position;
  d_previousVelocity = d_velocity;
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
    case BoundaryID::XMINUS:
      if (fabs(b_normalForce.x() / areaX + sigmaX) / sigmaX > tol) {
        vel = ((b_normalForce.x() + sigmaX * areaX) > 0 ? 1 : -1) * boundaryRate;
        // vel = d_previousVelocity.x() * (2-atf) / (2+atf) + (b_normalForce.x() + sigma *
        // areaX) / mass * timeStep * 2 / (2 + atf);
        pos = d_previousPosition.x() + vel * timeStep;
        setVelocity(Vec(vel, getVelocity().y(), getVelocity().z()));
        setPosition(Vec(pos, getPosition().y(), getPosition().z()));
      }
      break;
    case BoundaryID::XPLUS:
      if (fabs(b_normalForce.x() / areaX - sigmaX) / sigmaX > tol) {
        vel = ((b_normalForce.x() - sigmaX * areaX) > 0 ? 1 : -1) * boundaryRate;
        // vel = d_previousVelocity.x() * (2-atf) / (2+atf) + (b_normalForce.x() - sigma *
        // areaX) / mass * timeStep * 2 / (2 + atf);
        pos = d_previousPosition.x() + vel * timeStep;
        setVelocity(Vec(vel, getVelocity().y(), getVelocity().z()));
        setPosition(Vec(pos, getPosition().y(), getPosition().z()));
      }
      break;
    case BoundaryID::YMINUS:
      if (fabs(b_normalForce.y() / areaY + sigmaY) / sigmaY > tol) {
        vel = ((b_normalForce.y() + sigmaY * areaY) > 0 ? 1 : -1) * boundaryRate;
        // vel = d_previousVelocity.y() * (2-atf) / (2+atf) + (b_normalForce.y() + sigma *
        // areaY) / mass * timeStep * 2 / (2 + atf);
        pos = d_previousPosition.y() + vel * timeStep;
        setVelocity(Vec(getVelocity().x(), vel, getVelocity().z()));
        setPosition(Vec(getPosition().x(), pos, getPosition().z()));
      }
      break;
    case BoundaryID::YPLUS:
      if (fabs(b_normalForce.y() / areaY - sigmaY) / sigmaY > tol) {
        vel = ((b_normalForce.y() - sigmaY * areaY) > 0 ? 1 : -1) * boundaryRate;
        // vel = d_previousVelocity.y() * (2-atf) / (2+atf) + (b_normalForce.y() - sigma *
        // areaY) / mass * timeStep * 2 / (2 + atf);
        pos = d_previousPosition.y() + vel * timeStep;
        setVelocity(Vec(getVelocity().x(), vel, getVelocity().z()));
        setPosition(Vec(getPosition().x(), pos, getPosition().z()));
      }
      break;
    case BoundaryID::ZMINUS:
      if (fabs(b_normalForce.z() / areaZ + sigma) / sigma > tol) {
        vel = ((b_normalForce.z() + sigma * areaZ) > 0 ? 1 : -1) * boundaryRate;
        // vel = d_previousVelocity.z() * (2-atf) / (2+atf) + (b_normalForce.z() + sigma *
        // areaZ) / mass * timeStep * 2 / (2 + atf);
        pos = d_previousPosition.z() + vel * timeStep;
        setVelocity(Vec(getVelocity().x(), getVelocity().y(), vel));
        setPosition(Vec(getPosition().x(), getPosition().y(), pos));
      }
      break;
    case BoundaryID::ZPLUS:
      if (fabs(b_normalForce.z() / areaZ - sigma) / sigma > tol) {
        vel = ((b_normalForce.z() - sigma * areaZ) > 0 ? 1 : -1) * boundaryRate;
        // vel = d_previousVelocity.z() * (2-atf) / (2+atf) + (b_normalForce.z() - sigma *
        // areaZ) / mass * timeStep * 2 / (2 + atf);
        pos = d_previousPosition.z() + vel * timeStep;
        setVelocity(Vec(getVelocity().x(), getVelocity().y(), vel));
        setPosition(Vec(getPosition().x(), getPosition().y(), pos));
      }
      break;
    default:
      break;
  }
  d_previousPosition = d_position;
  d_previousVelocity = d_velocity;
}
