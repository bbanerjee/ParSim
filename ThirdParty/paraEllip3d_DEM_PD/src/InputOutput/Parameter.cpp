#include <Core/Math/IntVec.h>
#include <InputOutput/Parameter.h>
#include <InputOutput/zenxml/xml.h>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>

using namespace dem;

bool
Parameter::readInXML(const std::string& inputFileName)
{
  // Read the input file
  zen::XmlDoc doc;
  try {
    std::cout << "Input file name= " << inputFileName << "\n";
    doc = zen::load(inputFileName);
  } catch (const zen::XmlFileError& err) {
    std::cout << "*ERROR** Could not read input file " << inputFileName << "\n";
    std::cout << "    Error # = " << err.lastError << "\n";
    return false;
  } catch (const zen::XmlParsingError& err) {
    std::cout << "*ERROR** Could not read input file " << inputFileName << "\n";
    std::cout << "    Parse Error in line: " << err.row + 1
              << " col: " << err.col << "\n";
    return false;
  }

  // Check whether this is the right type of input file
  if (doc.root().getNameAs<std::string>() != "Ellip3D_input") {
    std::cout << "*ERROR** Could not find tag <Ellip3D_input> in input file "
              << inputFileName << "\n";
    return false;
  }

  // Load the document into input proxy for easier element access
  zen::XmlIn ps(doc);

  // Set up a lambda function to trim strings
  auto trim = [](std::string& str) {
    str.erase(0, str.find_first_not_of(" \n\r\t"));
    str.erase(str.find_last_not_of(" \n\r\t") + 1);
    return str;
  };

  // Read the title
  std::string title;
  if (!ps["Meta"]["title"](title)) {
    std::cout << "*ERROR** Could not find simulation title in input file "
              << inputFileName << "\n";
    std::cout << "  Add the <title> tag inside a <Meta> tag\n";
    return false;
  }
  std::cout << "title = " << trim(title) << "\n";

  // Read the simulation type
  int simType = 0;
  if (!ps["SimulationType"](simType)) {
    std::cout << "*ERROR** Could not find simulation type in input file "
              << inputFileName << "\n";
    std::cout << "  Add the <SimulationType> tag.\n";
    return false;
  }
  std::cout << "simulationType = " << simType << "\n";
  parameter["simuType"] = simType;

  // Read the parallel setup
  std::string mpiProcStr;
  if (!ps["Parallel"]["mpiProc"](mpiProcStr)) {
    std::cout << "*ERROR** Could not find mpi proc info in input file "
              << inputFileName << "\n";
    std::cout << "  Add the <mpiProc> tag inside the <Parallel> tag.\n";
    return false;
  }
  IntVec mpiProc = IntVec::fromString(mpiProcStr);
  std::cout << "mpiProcX = " << mpiProc.getX()
            << " mpiProcY = " << mpiProc.getY()
            << " mpiProcZ = " << mpiProc.getZ() << "\n";
  parameter["mpiProcX"] = mpiProc.getX();
  parameter["mpiProcY"] = mpiProc.getY();
  parameter["mpiProcZ"] = mpiProc.getZ();

  int ompThreads = 1;
  if (!ps["Parallel"]["ompThreads"](ompThreads)) {
    std::cout << "*ERROR** Could not find omp thread info in input file "
              << inputFileName << "\n";
    std::cout << "  Add the <ompThreads> tag inside the <Parallel> tag.\n";
    return false;
  }
  std::cout << "ompThreads = " << ompThreads << "\n";
  parameter["ompThreads"] = ompThreads;

  // Read time stepping info
  int startStep = 1;
  int endStep = 100;
  double timeAccrued = 0.0;
  double timeStep = 0.1;
  ps["Time"]["startStep"](startStep);
  ps["Time"]["endStep"](endStep);
  ps["Time"]["timeAccrued"](timeAccrued);
  ps["Time"]["timeStep"](timeStep);

  parameter["startStep"] = startStep;
  parameter["endStep"] = endStep;
  parameter["timeAccrued"] = timeAccrued;
  parameter["timeStep"] = timeStep;
  std::cout << "startStep = " << startStep << "\n"
            << "endStep = " << endStep << "\n"
            << "timeAccrued = " << timeAccrued << "\n"
            << "timeStep = " << timeStep << "\n";

  // Read the output info
  int startSnapshot = 1;
  int endSnapshot = 100;
  ps["Output"]["startSnapshot"](startSnapshot);
  ps["Output"]["endSnapshot"](endSnapshot);

  parameter["startSnap"] = startSnapshot;
  parameter["endSnap"] = endSnapshot;
  std::cout << "startSnapshot = " << startSnapshot << "\n"
            << "endSnapshot = " << endSnapshot << "\n";

  // Read the physical constants
  double gravity = 9.8;
  double gravityScale = 0.0;
  ps["PhysicalConstants"]["gravityAcceleration"](gravity);
  ps["PhysicalConstants"]["gravityScaleFactor"](gravityScale);

  parameter["gravAccel"] = gravity;
  parameter["gravScale"] = gravityScale;
  std::cout << "gravityAcc = " << gravity << "\n"
            << "gravityScale = " << gravityScale << "\n";

  // Read the boundary information
  std::string boundaryFile;
  double boundaryFriction = 0.0;
  ps["Boundary"]["boundaryFile"](boundaryFile);
  ps["Boundary"]["boundaryFriction"](boundaryFriction);

  datafile["boundaryFile"] = trim(boundaryFile);
  parameter["boundaryFric"] = boundaryFriction;
  std::cout << "boundaryFile = " << trim(boundaryFile) << "\n"
            << "boundaryFriction = " << boundaryFriction << "\n";

  // Read the DEM base information
  std::string particleFile;
  double massScaleFactor = 1.0;
  double momentScaleFactor = 1.0;
  ps["DEM"]["particleFile"](particleFile);
  ps["DEM"]["massScaleFactor"](massScaleFactor);
  ps["DEM"]["momentScaleFactor"](momentScaleFactor);

  datafile["particleFile"] = trim(particleFile);
  parameter["massScale"] = massScaleFactor;
  parameter["mntScale"] = momentScaleFactor;

  std::cout << "particleFile = " << trim(particleFile) << "\n"
            << "massScaleFactor = " << massScaleFactor << "\n"
            << "momentScaleFactor = " << momentScaleFactor << "\n";

  // Read the DEM material information
  double youngModulus = 1.0e10;
  double poissonRatio = 0.3;
  double specificGravity = 1.0;
  double membraneYoungModulus = 1.0e6;
  double forceDamping = 0.0;
  double momentDamping = 0.0;
  ps["DEM"]["Material"]["youngModulus"](youngModulus);
  ps["DEM"]["Material"]["poissonRatio"](poissonRatio);
  ps["DEM"]["Material"]["specificGravity"](specificGravity);
  ps["DEM"]["Material"]["membraneYoungModulus"](membraneYoungModulus);
  ps["DEM"]["Material"]["forceDamping"](forceDamping);
  ps["DEM"]["Material"]["momentDamping"](momentDamping);

  parameter["young"] = youngModulus;
  parameter["poisson"] = poissonRatio;
  parameter["specificG"] = specificGravity;
  parameter["memYoung"] = membraneYoungModulus;
  parameter["forceDamp"] = forceDamping;
  parameter["momentDamp"] = momentDamping;

  std::cout << "youngModulus = " << youngModulus << "\n"
            << " poissonRatio = " << poissonRatio << "\n"
            << " specificGravity = " << specificGravity << "\n"
            << " membraneYoungModulus = " << membraneYoungModulus << "\n"
            << " forceDamping = " << forceDamping << "\n"
            << " momentDamping = " << momentDamping << "\n";

  // Read the DEM contact information
  ps["DEM"]["Contact"]["contactDamping"](parameter["contactDamp"]);
  ps["DEM"]["Contact"]["contactFriction"](parameter["contactFric"]);
  ps["DEM"]["Contact"]["contactCohesion"](parameter["contactCohesion"]);
  ps["DEM"]["Contact"]["minRelativeOverlap"](parameter["minRelaOverlap"]);
  ps["DEM"]["Contact"]["maxRelativeOverlap"](parameter["maxRelaOverlap"]);
  ps["DEM"]["Contact"]["measurableOverlap"](parameter["measureOverlap"]);

  std::cout << "contactDamping = " << parameter["contactDamp"] << "\n"
            << "contactFriction = " << parameter["contactFric"] << "\n"
            << "contactCohesion = " << parameter["contactCohesion"] << "\n"
            << "minRelativeOverlap = " << parameter["minRelaOverlap"] << "\n"
            << "maxRelativeOverlap = " << parameter["maxRelaOverlap"] << "\n"
            << "measurableOverlap = " << parameter["measureOverlap"] << "\n";

  // Check if a peridynamics section exists
  auto peri_ps = ps["Peridynamics"];
  if (peri_ps) {

    std::string periFile;
    peri_ps["periFile"](periFile);
    datafile["periFile"] = trim(periFile);
    std::cout << "periFile = " << trim(periFile) << "\n";

    int initializeFromFile = 0;
    peri_ps["initializeFromFile"](initializeFromFile);
    std::cout << "initializeFromFile = " << initializeFromFile << "\n";
    parameter["toInitParticle"] = initializeFromFile;

    std::string intvecStr;
    if (!peri_ps["minPeriDomain"](intvecStr)) {
      std::cout
        << "*ERROR** Could not find min peridynamic domain info in input file "
        << inputFileName << "\n";
      std::cout
        << "  Add the <minPeriDomain> tag inside the <Peridynamics> tag.\n";
      return false;
    }
    IntVec minPeriDomain = IntVec::fromString(intvecStr);
    std::cout << "minPeriX = " << minPeriDomain.getX()
              << " minPeriY = " << minPeriDomain.getY()
              << " minPeriZ = " << minPeriDomain.getZ() << "\n";
    parameter["Xmin"] = minPeriDomain.getX();
    parameter["Ymin"] = minPeriDomain.getY();
    parameter["Zmin"] = minPeriDomain.getZ();

    if (!peri_ps["maxPeriDomain"](intvecStr)) {
      std::cout
        << "*ERROR** Could not find max peridynamic domain info in input file "
        << inputFileName << "\n";
      std::cout
        << "  Add the <maxPeriDomain> tag inside the <Peridynamics> tag.\n";
      return false;
    }
    IntVec maxPeriDomain = IntVec::fromString(intvecStr);
    std::cout << "maxPeriX = " << maxPeriDomain.getX()
              << " maxPeriY = " << maxPeriDomain.getY()
              << " maxPeriZ = " << maxPeriDomain.getZ() << "\n";
    parameter["Xmax"] = maxPeriDomain.getX();
    parameter["Ymax"] = maxPeriDomain.getY();
    parameter["Zmax"] = maxPeriDomain.getZ();

    // Peridynamics material properties
    auto peri_mat_ps = peri_ps["Material"];
    if (!peri_mat_ps) {
      std::cout
        << "*ERROR** No peridynamics material properties found in input file "
        << inputFileName << "\n";
      std::cout << "  Add the <Material> tag inside the <Peridynamics> tag.\n";
      return false;
    }

    peri_mat_ps["periDensity"](parameter["periDensity"]);
    peri_mat_ps["bodyDensity"](parameter["bodyDensity"]);
    peri_mat_ps["hchi"](parameter["hchi"]);
    peri_mat_ps["chi"](parameter["chi"]);
    peri_mat_ps["c"](parameter["c"]);
    peri_mat_ps["phi"](parameter["phi"]);
    peri_mat_ps["psi"](parameter["psi"]);
    peri_mat_ps["kappa"](parameter["kappa"]);
    peri_mat_ps["rEllip"](parameter["rEllip"]);
    peri_mat_ps["beta"](parameter["beta"]);
    peri_mat_ps["bondStretchLimit"](parameter["bondStretchLimit"]);

    std::cout << "periDensity " << parameter["periDensity"] << "\n"
              << " bodyDensity " << parameter["bodyDensity"] << "\n"
              << " hchi " << parameter["hchi"] << "\n"
              << " chi " << parameter["chi"] << "\n"
              << " c " << parameter["c"] << "\n"
              << " phi " << parameter["phi"] << "\n"
              << " psi " << parameter["psi"] << "\n"
              << " kappa " << parameter["kappa"] << "\n"
              << " rEllip " << parameter["rEllip"] << "\n"
              << " beta " << parameter["beta"] << "\n"
              << " bondStretchLimit " << parameter["bondStretchLimit"] << "\n";

    // Peridynamics constitutive model
    auto peri_cm_ps = peri_mat_ps["constitutive_model"];
    if (!peri_cm_ps) {
      std::cout << "*ERROR** No peridynamics material constitutive model found "
                   "in input file "
                << inputFileName << "\n";
      std::cout
        << "  Add the <constitutive_model> tag inside the <Material> tag "
        << " inside the <Peridynamics> tag.\n";
      return false;
    }

    // Get the material model type
    std::string model_type;
    if (!peri_cm_ps.attribute("type", model_type)) {
      std::cout
        << "**ERROR** Peridynamics constitutive model type not provided."
        << " Specify <constitutive_model type=\"xxxx\" >"
        << "\n";
      return false;
    }

    if (model_type == "linear_elastic") {
      peri_cm_ps["poissonRatio"](parameter["periPoisson"]);
      peri_cm_ps["youngModulus"](parameter["periYoung"]);
    } else {
      std::cout << "**ERROR** Only linear_elastic models are allowed\n";
      return false;
    }
    std::cout << "periPoisson = " << parameter["periPoisson"] << "\n"
              << "periYoung = " << parameter["periYoung"] << "\n";

    parameter["lambda"] = parameter["periPoisson"] * parameter["periYoung"] /
                          ((1.0 + parameter["periPoisson"]) *
                           (1.0 - 2.0 * parameter["periPoisson"]));
    parameter["mu"] =
      parameter["periYoung"] / (2.0 * (1.0 + parameter["periPoisson"]));
    parameter["kBulk"] =
      parameter["periYoung"] / (3.0 * (1.0 - 2.0 * parameter["periPoisson"]));
    parameter["tangentModulus11"] = parameter["lambda"] + 2.0 * parameter["mu"];
    parameter["tangentModulus12"] = parameter["lambda"];
    parameter["tangentModulus13"] = parameter["lambda"];
    parameter["tangentModulus21"] = parameter["lambda"];
    parameter["tangentModulus22"] = parameter["lambda"] + 2.0 * parameter["mu"];
    parameter["tangentModulus23"] = parameter["lambda"];
    parameter["tangentModulus31"] = parameter["lambda"];
    parameter["tangentModulus32"] = parameter["lambda"];
    parameter["tangentModulus33"] = parameter["lambda"] + 2.0 * parameter["mu"];
    parameter["tangentModulus44"] = parameter["mu"];
    parameter["tangentModulus55"] = parameter["mu"];
    parameter["tangentModulus66"] = parameter["mu"];
    parameter["Aphi"] = 2 * sqrt(6.0) * cos(parameter["phi"]) /
                        (3.0 + parameter["beta"] * sin(parameter["phi"]));
    parameter["Bphi"] = 2 * sqrt(6.0) * sin(parameter["phi"]) /
                        (3.0 + parameter["beta"] * sin(parameter["phi"]));
    parameter["Apsi"] = 2 * sqrt(6.0) * cos(parameter["psi"]) /
                        (3.0 + parameter["beta"] * sin(parameter["psi"]));
    parameter["Bpsi"] = 2 * sqrt(6.0) * sin(parameter["psi"]) /
                        (3.0 + parameter["beta"] * sin(parameter["psi"]));
  }

  return true;
}

void
Parameter::writeOutXML()
{
}

void
Parameter::readIn(const char* input)
{
  if (readInXML(input))
    return;

  std::cout << "**WARNING** Failed to read XML input file " << input << "\n";
  std::cout << "            Trying to read the file as ordinary text\n";

  std::ifstream ifs;
  ifs.open(input);
  if (!ifs) {
    std::cout << "stream error: Parameter.cpp" << std::endl;
    exit(-1);
  }
  std::string line;
  std::istringstream ssline;
  std::string str, str2;
  REAL val;

  // 28 generic parameters
  for (std::size_t i = 0; i < 28; ++i) {
    while (getline(ifs, line))
      if (line[0] != '#' && line.compare("") != 0)
        break;
    ssline.clear();
    ssline.str(line);
    ssline >> str >> val;
    parameter[str] = val;
  }

  // for different types of simulation
  std::size_t simuType = static_cast<std::size_t>(parameter["simuType"]);
  switch (simuType) {
    case 001: // proceed from preset state
      for (std::size_t i = 0; i < 2; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> str2;
        datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 1; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> val;
        parameter[str] = val;
      }
      break;

    case 002: // tuneMassPercentage
      for (std::size_t i = 0; i < 12; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> val;
        parameter[str] = val;
      }
      for (std::size_t i = 0;
           i < static_cast<std::size_t>(parameter["sieveNum"]); ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        REAL percent, size;
        ssline >> percent >> size;
        gradation.push_back(std::make_pair(percent, size));
      }
      break;

    case 003: // trimOnly
      for (std::size_t i = 0; i < 2; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> str2;
        datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 1; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> val;
        parameter[str] = val;
      }
      break;

    case 101: // depositIntoContainer
      for (std::size_t i = 0; i < 12; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> val;
        parameter[str] = val;
      }
      for (std::size_t i = 0;
           i < static_cast<std::size_t>(parameter["sieveNum"]); ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        REAL percent, size;
        ssline >> percent >> size;
        gradation.push_back(std::make_pair(percent, size));
      }
      break;

    case 102: // resumeDepositIntoContainer
      for (std::size_t i = 0; i < 2; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> str2;
        datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 2; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> val;
        parameter[str] = val;
      }
      break;

    case 201: // isotropic 1
      for (std::size_t i = 0; i < 2; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> str2;
        datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 6; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> val;
        parameter[str] = val;
      }
      break;

    case 202: // isotropic 2
      for (std::size_t i = 0; i < 2; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> str2;
        datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 7; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> val;
        parameter[str] = val;
      }
      break;

    case 203: // isotropic 3
      for (std::size_t i = 0; i < 2; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> str2;
        datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 3; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> val;
        parameter[str] = val;
      }
      for (std::size_t i = 0;
           i < static_cast<std::size_t>(parameter["sigmaPoints"]); ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        REAL sigma;
        ssline >> sigma;
        sigmaPath.push_back(sigma);
      }
      for (std::size_t i = 0; i < 3; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> val;
        parameter[str] = val;
      }
      break;

    case 301: // odometer 1
      for (std::size_t i = 0; i < 2; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> str2;
        datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 7; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> val;
        parameter[str] = val;
      }
      break;

    case 302: // odometer 2
      for (std::size_t i = 0; i < 2; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> str2;
        datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 3; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> val;
        parameter[str] = val;
      }
      for (std::size_t i = 0;
           i < static_cast<std::size_t>(parameter["sigmaPoints"]); ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        REAL sigma;
        ssline >> sigma;
        sigmaPath.push_back(sigma);
      }
      for (std::size_t i = 0; i < 3; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> val;
        parameter[str] = val;
      }
      break;

    case 401: // triaxial 1
      for (std::size_t i = 0; i < 2; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> str2;
        datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 5; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> val;
        parameter[str] = val;
      }
      break;

    case 402: // triaxial 2
      for (std::size_t i = 0; i < 2; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> str2;
        datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 6; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> val;
        parameter[str] = val;
      }
      break;

    case 411: // plane strain 1
      for (std::size_t i = 0; i < 2; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> str2;
        datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 6; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> val;
        parameter[str] = val;
      }
      break;

    case 412: // plain strain 2
      for (std::size_t i = 0; i < 2; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> str2;
        datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 7; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> val;
        parameter[str] = val;
      }
      break;

    case 501: // true triaxial 1
      for (std::size_t i = 0; i < 2; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> str2;
        datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 9; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> val;
        parameter[str] = val;
      }
      break;

    case 502: // true triaxial 2
      for (std::size_t i = 0; i < 2; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> str2;
        datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 10; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> val;
        parameter[str] = val;
      }
      break;

    case 601: // expandCavityParticle
      for (std::size_t i = 0; i < 2; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> str2;
        datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 7; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> val;
        parameter[str] = val;
      }
      break;

    case 602: // resumeExpandCavityParticle
      for (std::size_t i = 0; i < 2; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> str2;
        datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 1; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> val;
        parameter[str] = val;
      }
      break;

    case 701: // couple with sonic fluid flow
      for (std::size_t i = 0; i < 2; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> str2;
        datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 19; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> val;
        parameter[str] = val;
      }
      break;

    case 3001: // couple with peridynamics clay
      for (std::size_t i = 0; i < 2; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> str2;
        datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 26; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> val;
        parameter[str] = val;
      }
      parameter["lambda"] = parameter["periPoisson"] * parameter["periYoung"] /
                            ((1.0 + parameter["periPoisson"]) *
                             (1.0 - 2.0 * parameter["periPoisson"]));
      parameter["mu"] =
        parameter["periYoung"] / (2.0 * (1.0 + parameter["periPoisson"]));
      parameter["kBulk"] =
        parameter["periYoung"] / (3.0 * (1.0 - 2.0 * parameter["periPoisson"]));
      parameter["tangentModulus11"] =
        parameter["lambda"] + 2.0 * parameter["mu"];
      parameter["tangentModulus12"] = parameter["lambda"];
      parameter["tangentModulus13"] = parameter["lambda"];
      parameter["tangentModulus21"] = parameter["lambda"];
      parameter["tangentModulus22"] =
        parameter["lambda"] + 2.0 * parameter["mu"];
      parameter["tangentModulus23"] = parameter["lambda"];
      parameter["tangentModulus31"] = parameter["lambda"];
      parameter["tangentModulus32"] = parameter["lambda"];
      parameter["tangentModulus33"] =
        parameter["lambda"] + 2.0 * parameter["mu"];
      parameter["tangentModulus44"] = parameter["mu"];
      parameter["tangentModulus55"] = parameter["mu"];
      parameter["tangentModulus66"] = parameter["mu"];
      parameter["Aphi"] = 2 * sqrt(6.0) * cos(parameter["phi"]) /
                          (3.0 + parameter["beta"] * sin(parameter["phi"]));
      parameter["Bphi"] = 2 * sqrt(6.0) * sin(parameter["phi"]) /
                          (3.0 + parameter["beta"] * sin(parameter["phi"]));
      parameter["Apsi"] = 2 * sqrt(6.0) * cos(parameter["psi"]) /
                          (3.0 + parameter["beta"] * sin(parameter["psi"]));
      parameter["Bpsi"] = 2 * sqrt(6.0) * sin(parameter["psi"]) /
                          (3.0 + parameter["beta"] * sin(parameter["psi"]));

      break;

    case 3002: // couple with peridynamics clay, pull out DEM particles in
               // peri-domain
      for (std::size_t i = 0; i < 2; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> str2;
        datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 21; ++i) {
        while (getline(ifs, line))
          if (line[0] != '#' && line.compare("") != 0)
            break;
        ssline.clear();
        ssline.str(line);
        ssline >> str >> val;
        parameter[str] = val;
      }
      parameter["lambda"] = parameter["periPoisson"] * parameter["periYoung"] /
                            ((1.0 + parameter["periPoisson"]) *
                             (1.0 - 2.0 * parameter["periPoisson"]));
      parameter["mu"] =
        parameter["periYoung"] / (2.0 * (1.0 + parameter["periPoisson"]));
      parameter["kBulk"] =
        parameter["periYoung"] / (3.0 * (1.0 - 2.0 * parameter["periPoisson"]));
      parameter["tangentModulus11"] =
        parameter["lambda"] + 2.0 * parameter["mu"];
      parameter["tangentModulus12"] = parameter["lambda"];
      parameter["tangentModulus13"] = parameter["lambda"];
      parameter["tangentModulus21"] = parameter["lambda"];
      parameter["tangentModulus22"] =
        parameter["lambda"] + 2.0 * parameter["mu"];
      parameter["tangentModulus23"] = parameter["lambda"];
      parameter["tangentModulus31"] = parameter["lambda"];
      parameter["tangentModulus32"] = parameter["lambda"];
      parameter["tangentModulus33"] =
        parameter["lambda"] + 2.0 * parameter["mu"];
      parameter["tangentModulus44"] = parameter["mu"];
      parameter["tangentModulus55"] = parameter["mu"];
      parameter["tangentModulus66"] = parameter["mu"];
      parameter["Aphi"] = 2 * sqrt(6.0) * cos(parameter["phi"]) /
                          (3.0 + parameter["beta"] * sin(parameter["phi"]));
      parameter["Bphi"] = 2 * sqrt(6.0) * sin(parameter["phi"]) /
                          (3.0 + parameter["beta"] * sin(parameter["phi"]));
      parameter["Apsi"] = 2 * sqrt(6.0) * cos(parameter["psi"]) /
                          (3.0 + parameter["beta"] * sin(parameter["psi"]));
      parameter["Bpsi"] = 2 * sqrt(6.0) * sin(parameter["psi"]) /
                          (3.0 + parameter["beta"] * sin(parameter["psi"]));

      break;
  }

  ifs.close();
}

void
Parameter::writeOut()
{
  std::map<std::string, REAL>& param = dem::Parameter::getSingleton().parameter;
  std::vector<std::pair<REAL, REAL>>& grada =
    dem::Parameter::getSingleton().gradation;
  std::map<std::string, std::string>& file =
    dem::Parameter::getSingleton().datafile;
  std::vector<REAL>& sigma = dem::Parameter::getSingleton().sigmaPath;

  for (std::map<std::string, REAL>::const_iterator it = param.begin();
       it != param.end(); ++it)
    std::cout << it->first << "  " << it->second << std::endl;

  for (auto& i : grada)
    std::cout << i.first << "  " << i.second << std::endl;

  for (std::map<std::string, std::string>::const_iterator it = file.begin();
       it != file.end(); ++it)
    std::cout << it->first << "  " << it->second << std::endl;

  for (std::vector<REAL>::const_iterator it = sigma.begin(); it != sigma.end();
       ++it)
    std::cout << (*it) << std::endl;
}
