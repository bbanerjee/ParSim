#include <Core/Math/IntVec.h>
#include <Core/Math/Vec.h>
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
    //std::cout << "Input file name= " << inputFileName << "\n";
    doc = zen::load(inputFileName);
  } catch (const zen::XmlFileError& err) {
    std::cerr << "*ERROR** Could not read input file " << inputFileName << "\n";
    std::cerr << "    Error # = " << err.lastError << "\n";
    return false;
  } catch (const zen::XmlParsingError& err) {
    std::cerr << "*ERROR** Could not read input file " << inputFileName << "\n";
    std::cerr << "    Parse Error in line: " << err.row + 1
              << " col: " << err.col << "\n";
    return false;
  }

  // Check whether this is the right type of input file
  if (doc.root().getNameAs<std::string>() != "Ellip3D_input") {
    std::cerr << "*ERROR** Could not find tag <Ellip3D_input> in input file "
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
    std::cerr << "*ERROR** Could not find simulation title in input file "
              << inputFileName << "\n";
    std::cerr << "  Add the <title> tag inside a <Meta> tag\n";
    return false;
  }
  //std::cout << "title = " << trim(title) << "\n";

  // Read the simulation type
  int simType = 0;
  if (!ps["SimulationType"](simType)) {
    std::cerr << "*ERROR** Could not find simulation type in input file "
              << inputFileName << "\n";
    std::cerr << "  Add the <SimulationType> tag.\n";
    return false;
  }
  //std::cout << "simulationType = " << simType << "\n";
  param["simuType"] = simType;

  // Read the parallel setup
  std::string mpiProcStr;
  if (!ps["Parallel"]["mpiProc"](mpiProcStr)) {
    std::cerr << "*ERROR** Could not find mpi proc info in input file "
              << inputFileName << "\n";
    std::cerr << "  Add the <mpiProc> tag inside the <Parallel> tag.\n";
    return false;
  }
  IntVec mpiProc = IntVec::fromString(mpiProcStr);
  //std::cout << "mpiProcX = " << mpiProc.x() << " mpiProcY = " << mpiProc.y()
  //          << " mpiProcZ = " << mpiProc.z() << "\n";
  param["mpiProcX"] = mpiProc.x();
  param["mpiProcY"] = mpiProc.y();
  param["mpiProcZ"] = mpiProc.z();

  int ompThreads = 1;
  if (!ps["Parallel"]["ompThreads"](ompThreads)) {
    std::cerr << "*ERROR** Could not find omp thread info in input file "
              << inputFileName << "\n";
    std::cerr << "  Add the <ompThreads> tag inside the <Parallel> tag.\n";
    return false;
  }
  //std::cout << "ompThreads = " << ompThreads << "\n";
  param["ompThreads"] = ompThreads;

  // Read time stepping info
  int startStep = 1;
  int endStep = 100;
  double timeAccrued = 0.0;
  double timeStep = 0.1;
  ps["Time"]["startStep"](startStep);
  ps["Time"]["endStep"](endStep);
  ps["Time"]["timeAccrued"](timeAccrued);
  ps["Time"]["timeStep"](timeStep);

  param["startStep"] = startStep;
  param["endStep"] = endStep;
  param["timeAccrued"] = timeAccrued;
  param["timeStep"] = timeStep;
  //std::cout << "startStep = " << startStep << "\n"
  //          << "endStep = " << endStep << "\n"
  //          << "timeAccrued = " << timeAccrued << "\n"
  //          << "timeStep = " << timeStep << "\n";

  // Read the output info
  int startSnapshot = 1;
  int endSnapshot = 100;
  std::string outputFolderName = "deposit";
  ps["Output"]["outputFolder"](outputFolderName);
  ps["Output"]["startSnapshot"](startSnapshot);
  ps["Output"]["endSnapshot"](endSnapshot);

  datafile["outputFolder"] = trim(outputFolderName);
  param["startSnap"] = startSnapshot;
  param["endSnap"] = endSnapshot;
  //std::cout << "startSnapshot = " << startSnapshot << "\n"
  //          << "endSnapshot = " << endSnapshot << "\n";

  // Read the physical constants
  double gravity = 9.8;
  double gravityScale = 0.0;
  ps["PhysicalConstants"]["gravityAcceleration"](gravity);
  ps["PhysicalConstants"]["gravityScaleFactor"](gravityScale);

  param["gravAccel"] = gravity;
  param["gravScale"] = gravityScale;
  //std::cout << "gravityAcc = " << gravity << "\n"
  //          << "gravityScale = " << gravityScale << "\n";

  // Read the boundary information
  std::string boundaryFile;
  double boundaryFriction = 0.0;
  ps["Boundary"]["boundaryFile"](boundaryFile);
  ps["Boundary"]["boundaryFriction"](boundaryFriction);

  datafile["boundaryFile"] = trim(boundaryFile);
  param["boundaryFric"] = boundaryFriction;
  //std::cout << "boundaryFile = " << trim(boundaryFile) << "\n"
  //          << "boundaryFriction = " << boundaryFriction << "\n";

  // Read the DEM base information
  std::string particleFile;
  double massScaleFactor = 1.0;
  double momentScaleFactor = 1.0;
  double pileRate = 0.0;
  ps["DEM"]["particleFile"](particleFile);
  ps["DEM"]["massScaleFactor"](massScaleFactor);
  ps["DEM"]["momentScaleFactor"](momentScaleFactor);
  ps["DEM"]["pileRate"](pileRate);

  datafile["particleFile"] = trim(particleFile);
  param["massScale"] = massScaleFactor;
  param["mntScale"] = momentScaleFactor;
  param["pileRate"] = pileRate;

  //std::cout << "particleFile = " << trim(particleFile) << "\n"
  //          << "massScaleFactor = " << massScaleFactor << "\n"
  //          << "momentScaleFactor = " << momentScaleFactor << "\n";

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

  param["young"] = youngModulus;
  param["poisson"] = poissonRatio;
  param["specificG"] = specificGravity;
  param["memYoung"] = membraneYoungModulus;
  param["forceDamp"] = forceDamping;
  param["momentDamp"] = momentDamping;

  //std::cout << "youngModulus = " << youngModulus << "\n"
  //          << " poissonRatio = " << poissonRatio << "\n"
  //          << " specificGravity = " << specificGravity << "\n"
  //          << " membraneYoungModulus = " << membraneYoungModulus << "\n"
  //          << " forceDamping = " << forceDamping << "\n"
  //          << " momentDamping = " << momentDamping << "\n";

  // Read the DEM contact information
  ps["DEM"]["Contact"]["contactDamping"](param["contactDamp"]);
  ps["DEM"]["Contact"]["contactFriction"](param["contactFric"]);
  ps["DEM"]["Contact"]["contactCohesion"](param["contactCohesion"]);
  ps["DEM"]["Contact"]["minRelativeOverlap"](param["minRelaOverlap"]);
  ps["DEM"]["Contact"]["maxRelativeOverlap"](param["maxRelaOverlap"]);
  ps["DEM"]["Contact"]["measurableOverlap"](param["measureOverlap"]);

  //std::cout << "contactDamping = " << param["contactDamp"] << "\n"
  //          << "contactFriction = " << param["contactFric"] << "\n"
  //          << "contactCohesion = " << param["contactCohesion"] << "\n"
  //          << "minRelativeOverlap = " << param["minRelaOverlap"] << "\n"
  //          << "maxRelativeOverlap = " << param["maxRelaOverlap"] << "\n"
  //          << "measurableOverlap = " << param["measureOverlap"] << "\n";

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
    param["toInitParticle"] = initializeFromFile;

    std::string vecStr;
    if (!peri_ps["minPeriDomain"](vecStr)) {
      std::cerr
        << "*ERROR** Could not find min peridynamic domain info in input file "
        << inputFileName << "\n";
      std::cerr
        << "  Add the <minPeriDomain> tag inside the <Peridynamics> tag.\n";
      return false;
    }
    Vec minPeriDomain = Vec::fromString(vecStr);
    std::cout << "minPeriX = " << minPeriDomain.x()
              << " minPeriY = " << minPeriDomain.y()
              << " minPeriZ = " << minPeriDomain.z() << "\n";
    param["Xmin"] = minPeriDomain.x();
    param["Ymin"] = minPeriDomain.y();
    param["Zmin"] = minPeriDomain.z();

    if (!peri_ps["maxPeriDomain"](vecStr)) {
      std::cerr
        << "*ERROR** Could not find max peridynamic domain info in input file "
        << inputFileName << "\n";
      std::cerr
        << "  Add the <maxPeriDomain> tag inside the <Peridynamics> tag.\n";
      return false;
    }
    Vec maxPeriDomain = Vec::fromString(vecStr);
    std::cout << "maxPeriX = " << maxPeriDomain.x()
              << " maxPeriY = " << maxPeriDomain.y()
              << " maxPeriZ = " << maxPeriDomain.z() << "\n";
    param["Xmax"] = maxPeriDomain.x();
    param["Ymax"] = maxPeriDomain.y();
    param["Zmax"] = maxPeriDomain.z();

    // Peridynamics material properties
    auto peri_mat_ps = peri_ps["Material"];
    if (!peri_mat_ps) {
      std::cerr
        << "*ERROR** No peridynamics material properties found in input file "
        << inputFileName << "\n";
      std::cerr << "  Add the <Material> tag inside the <Peridynamics> tag.\n";
      return false;
    }

    peri_mat_ps["typeConstitutive"](param["typeConstitutive"]);

    peri_mat_ps["periDensity"](param["periDensity"]);
    peri_mat_ps["bodyDensity"](param["bodyDensity"]);
    peri_mat_ps["hchi"](param["hchi"]);
    peri_mat_ps["chi"](param["chi"]);
    peri_mat_ps["c"](param["c"]);
    peri_mat_ps["phi"](param["phi"]);
    peri_mat_ps["psi"](param["psi"]);
    peri_mat_ps["kappa"](param["kappa"]);
    peri_mat_ps["rEllip"](param["rEllip"]);
    peri_mat_ps["beta"](param["beta"]);
    peri_mat_ps["bondStretchLimit"](param["bondStretchLimit"]);

    std::cout << "periDensity " << param["periDensity"] << "\n"
              << " bodyDensity " << param["bodyDensity"] << "\n"
              << " hchi " << param["hchi"] << "\n"
              << " chi " << param["chi"] << "\n"
              << " c " << param["c"] << "\n"
              << " phi " << param["phi"] << "\n"
              << " psi " << param["psi"] << "\n"
              << " kappa " << param["kappa"] << "\n"
              << " rEllip " << param["rEllip"] << "\n"
              << " beta " << param["beta"] << "\n"
              << " bondStretchLimit " << param["bondStretchLimit"] << "\n";

    // Peridynamics constitutive model
    auto peri_cm_ps = peri_mat_ps["constitutive_model"];
    if (!peri_cm_ps) {
      std::cerr << "*ERROR** No peridynamics material constitutive model found "
                   "in input file "
                << inputFileName << "\n";
      std::cerr
        << "  Add the <constitutive_model> tag inside the <Material> tag "
        << " inside the <Peridynamics> tag.\n";
      return false;
    }

    // Get the material model type
    std::string model_type;
    if (!peri_cm_ps.attribute("type", model_type)) {
      std::cerr
        << "**ERROR** Peridynamics constitutive model type not provided."
        << " Specify <constitutive_model type=\"xxxx\" >"
        << "\n";
      return false;
    }

    if (model_type == "linear_elastic") {
      peri_cm_ps["poissonRatio"](param["periPoisson"]);
      peri_cm_ps["youngModulus"](param["periYoung"]);
    } else {
      std::cerr << "**ERROR** Only linear_elastic models are allowed\n";
      return false;
    }
    std::cout << "periPoisson = " << param["periPoisson"] << "\n"
              << "periYoung = " << param["periYoung"] << "\n";

    param["lambda"] =
      param["periPoisson"] * param["periYoung"] /
      ((1.0 + param["periPoisson"]) * (1.0 - 2.0 * param["periPoisson"]));
    param["mu"] = param["periYoung"] / (2.0 * (1.0 + param["periPoisson"]));
    param["kBulk"] =
      param["periYoung"] / (3.0 * (1.0 - 2.0 * param["periPoisson"]));
    param["tangentModulus11"] = param["lambda"] + 2.0 * param["mu"];
    param["tangentModulus12"] = param["lambda"];
    param["tangentModulus13"] = param["lambda"];
    param["tangentModulus21"] = param["lambda"];
    param["tangentModulus22"] = param["lambda"] + 2.0 * param["mu"];
    param["tangentModulus23"] = param["lambda"];
    param["tangentModulus31"] = param["lambda"];
    param["tangentModulus32"] = param["lambda"];
    param["tangentModulus33"] = param["lambda"] + 2.0 * param["mu"];
    param["tangentModulus44"] = param["mu"];
    param["tangentModulus55"] = param["mu"];
    param["tangentModulus66"] = param["mu"];
    param["Aphi"] = 2 * sqrt(6.0) * cos(param["phi"]) /
                    (3.0 + param["beta"] * sin(param["phi"]));
    param["Bphi"] = 2 * sqrt(6.0) * sin(param["phi"]) /
                    (3.0 + param["beta"] * sin(param["phi"]));
    param["Apsi"] = 2 * sqrt(6.0) * cos(param["psi"]) /
                    (3.0 + param["beta"] * sin(param["psi"]));
    param["Bpsi"] = 2 * sqrt(6.0) * sin(param["psi"]) /
                    (3.0 + param["beta"] * sin(param["psi"]));

    // Peridynamics BCs
    auto peri_bc_ps = peri_ps["BoundaryConditions"];
    if (!peri_bc_ps) {
      std::cerr
        << "*ERROR** No peridynamics displacement/load BC found in input file "
        << inputFileName << "\n";
      std::cerr << "  Add the <BoundaryConditions> tag inside the <Peridynamics> tag.\n";
      return false;
    }
    peri_bc_ps["fixRadius"](param["fixRadius"]);
    peri_bc_ps["periFixCentroidX"](param["periFixCentroidX"]);
    peri_bc_ps["periFixCentroidY"](param["periFixCentroidY"]);
    peri_bc_ps["periFixCentroidZ"](param["periFixCentroidZ"]);
    peri_bc_ps["periForce"](param["periForce"]);
    peri_bc_ps["rampStep"](param["rampStep"]);
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

  std::cerr << "**WARNING** Failed to read XML input file " << input << "\n";
  std::cerr << "            Trying to read the file as ordinary text\n";

  std::ifstream ifs;
  ifs.open(input);
  if (!ifs) {
    std::cerr << "stream error: Parameter.cpp" << std::endl;
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
    param[str] = val;
  }

  // for different types of simulation
  std::size_t simuType = static_cast<std::size_t>(param["simuType"]);
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
        param[str] = val;
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
        param[str] = val;
      }
      for (std::size_t i = 0; i < static_cast<std::size_t>(param["sieveNum"]);
           ++i) {
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
        param[str] = val;
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
        param[str] = val;
      }
      for (std::size_t i = 0; i < static_cast<std::size_t>(param["sieveNum"]);
           ++i) {
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
        param[str] = val;
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
        param[str] = val;
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
        param[str] = val;
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
        param[str] = val;
      }
      for (std::size_t i = 0;
           i < static_cast<std::size_t>(param["sigmaPoints"]); ++i) {
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
        param[str] = val;
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
        param[str] = val;
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
        param[str] = val;
      }
      for (std::size_t i = 0;
           i < static_cast<std::size_t>(param["sigmaPoints"]); ++i) {
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
        param[str] = val;
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
        param[str] = val;
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
        param[str] = val;
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
        param[str] = val;
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
        param[str] = val;
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
        param[str] = val;
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
        param[str] = val;
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
        param[str] = val;
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
        param[str] = val;
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
        param[str] = val;
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
        param[str] = val;
      }
      param["lambda"] =
        param["periPoisson"] * param["periYoung"] /
        ((1.0 + param["periPoisson"]) * (1.0 - 2.0 * param["periPoisson"]));
      param["mu"] = param["periYoung"] / (2.0 * (1.0 + param["periPoisson"]));
      param["kBulk"] =
        param["periYoung"] / (3.0 * (1.0 - 2.0 * param["periPoisson"]));
      param["tangentModulus11"] = param["lambda"] + 2.0 * param["mu"];
      param["tangentModulus12"] = param["lambda"];
      param["tangentModulus13"] = param["lambda"];
      param["tangentModulus21"] = param["lambda"];
      param["tangentModulus22"] = param["lambda"] + 2.0 * param["mu"];
      param["tangentModulus23"] = param["lambda"];
      param["tangentModulus31"] = param["lambda"];
      param["tangentModulus32"] = param["lambda"];
      param["tangentModulus33"] = param["lambda"] + 2.0 * param["mu"];
      param["tangentModulus44"] = param["mu"];
      param["tangentModulus55"] = param["mu"];
      param["tangentModulus66"] = param["mu"];
      param["Aphi"] = 2 * sqrt(6.0) * cos(param["phi"]) /
                      (3.0 + param["beta"] * sin(param["phi"]));
      param["Bphi"] = 2 * sqrt(6.0) * sin(param["phi"]) /
                      (3.0 + param["beta"] * sin(param["phi"]));
      param["Apsi"] = 2 * sqrt(6.0) * cos(param["psi"]) /
                      (3.0 + param["beta"] * sin(param["psi"]));
      param["Bpsi"] = 2 * sqrt(6.0) * sin(param["psi"]) /
                      (3.0 + param["beta"] * sin(param["psi"]));

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
        param[str] = val;
      }
      param["lambda"] =
        param["periPoisson"] * param["periYoung"] /
        ((1.0 + param["periPoisson"]) * (1.0 - 2.0 * param["periPoisson"]));
      param["mu"] = param["periYoung"] / (2.0 * (1.0 + param["periPoisson"]));
      param["kBulk"] =
        param["periYoung"] / (3.0 * (1.0 - 2.0 * param["periPoisson"]));
      param["tangentModulus11"] = param["lambda"] + 2.0 * param["mu"];
      param["tangentModulus12"] = param["lambda"];
      param["tangentModulus13"] = param["lambda"];
      param["tangentModulus21"] = param["lambda"];
      param["tangentModulus22"] = param["lambda"] + 2.0 * param["mu"];
      param["tangentModulus23"] = param["lambda"];
      param["tangentModulus31"] = param["lambda"];
      param["tangentModulus32"] = param["lambda"];
      param["tangentModulus33"] = param["lambda"] + 2.0 * param["mu"];
      param["tangentModulus44"] = param["mu"];
      param["tangentModulus55"] = param["mu"];
      param["tangentModulus66"] = param["mu"];
      param["Aphi"] = 2 * sqrt(6.0) * cos(param["phi"]) /
                      (3.0 + param["beta"] * sin(param["phi"]));
      param["Bphi"] = 2 * sqrt(6.0) * sin(param["phi"]) /
                      (3.0 + param["beta"] * sin(param["phi"]));
      param["Apsi"] = 2 * sqrt(6.0) * cos(param["psi"]) /
                      (3.0 + param["beta"] * sin(param["psi"]));
      param["Bpsi"] = 2 * sqrt(6.0) * sin(param["psi"]) /
                      (3.0 + param["beta"] * sin(param["psi"]));

      break;
  }

  ifs.close();
}

void
Parameter::writeOut()
{
  std::map<std::string, REAL>& param = Parameter::get().param;
  std::vector<std::pair<REAL, REAL>>& grada = Parameter::get().gradation;
  std::map<std::string, std::string>& file = Parameter::get().datafile;
  std::vector<REAL>& sigma = Parameter::get().sigmaPath;

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
