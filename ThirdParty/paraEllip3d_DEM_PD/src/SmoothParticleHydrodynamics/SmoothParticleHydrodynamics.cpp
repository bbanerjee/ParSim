#include <SmoothParticleHydrodynamics/SmoothParticleHydrodynamics.h>

#include <InputOutput/PeriParticleFileReader.h>
#include <InputOutput/OutputTecplot.h>
#include <InputOutput/OutputVTK.h>
#include <Core/Const/const.h>
#include <Core/Math/IntVec.h>
#include <chrono>

using namespace sph;

using Timer = std::chrono::steady_clock;
using Seconds = std::chrono::seconds;
using IntVec = dem::IntVec;
using Vec = dem::Vec;
using Matrix = dem::Matrix;
using Box = dem::Box;
using ParticlePArray = dem::ParticlePArray;
using OutputVTK = dem::OutputVTK<SPHParticlePArray>;
using OutputTecplot = dem::OutputTecplot<SPHParticlePArray>;

SmoothParticleHydrodynamics::SmoothParticleHydrodynamics() 
{
}

SmoothParticleHydrodynamics::~SmoothParticleHydrodynamics()
{
  allSPHParticleVec.clear();
  SPHParticleVec.clear();
} 

// here the free particles, ghost particles and boundary particles 
// will all stored in the same vector
// unlike the implementation in serial version.
void 
SmoothParticleHydrodynamics::generateSPHParticle2D(const dem::Box& allContainer)
{
  if (getMPIRank() != 0) return;  // make only primary cpu generate particles

  // sph parameters
  auto gravAccel   = util::getParam<REAL>("gravAccel");
  auto gravScale   = util::getParam<REAL>("gravScale");
  auto waterLength = util::getParam<REAL>("waterLength");
  auto numSPHPoint = util::getParam<REAL>("nSPHPoint");
  auto numLayers   = util::getParam<int>("numLayers");
  auto gamma       = util::getParam<REAL>("gamma");
  auto P0          = util::getParam<REAL>("P0");
  auto SPHInitialDensity = util::getParam<REAL>("SPHInitialDensity");

  REAL spaceInterval = waterLength/(numSPHPoint-1);
  REAL L_over_N = waterLength/numSPHPoint;
  REAL SPHmass = SPHInitialDensity*L_over_N*L_over_N;
  REAL small_value = 0.01*spaceInterval;
  REAL smoothLength = 1.5*spaceInterval;
  REAL kernelSize = 3*smoothLength;

  // get the dimensions of the sph domain
  Vec vmin = allContainer.getMinCorner();
  Vec vmax = allContainer.getMaxCorner();
  REAL xmin = vmin.getX();
  REAL ymin = vmin.getY();
  REAL zmin = vmin.getZ();
  REAL xmax = vmax.getX();
  REAL ymax = vmax.getY();
  REAL zmax = vmax.getZ();

  // Create the domain buffer length
  REAL bufferLength = spaceInterval*numLayers;

  // Modify the domain size using the buffer length
  REAL xminBuffered = xmin - bufferLength;
  REAL xmaxBuffered = xmax + bufferLength;
  REAL zminBuffered = zmin - bufferLength;
  REAL zmaxBuffered = zmax + bufferLength;

  // Create an linearly spaced array of zcoords from zminBuffered to zmaxBuffered

  // create and store SPHParticle objects into sphParticleVec
  int isGhost = 0;  // is Ghost particle
  REAL radius_a, radius_b, radius_c;
  REAL inner_a,  inner_b,  inner_c;
  dem::Vec pt_position;
  dem::Vec local_x;  // local position of ghost point in dem particle
  dem::Vec tmp_xyz, tmp_local;

  for(REAL tmp_z=zmin-spaceInterval*(numLayers); tmp_z<zmax+spaceInterval*numLayers; tmp_z=tmp_z+spaceInterval){
    for(REAL tmp_x=xmin-spaceInterval*(numLayers); tmp_x<xmax+spaceInterval*numLayers; tmp_x=tmp_x+spaceInterval){
    isGhost = 0;
    tmp_xyz = dem::Vec(tmp_x, 0, tmp_z);
      for(std::vector<Particle*>::iterator pt=allParticleVec.begin(); pt!=allParticleVec.end(); pt++){
      radius_a = (*pt)->getA(); radius_b = (*pt)->getB(); radius_c = (*pt)->getC();
      inner_a = radius_a-kernelSize; inner_b = radius_b-kernelSize; inner_c = radius_c-kernelSize;
      pt_position = (*pt)->getCurrPos();
      pt_position.setY(0);
      tmp_local = (*pt)->globalToLocal(tmp_xyz-pt_position);
      if( tmp_local.getX()*tmp_local.getX()/(radius_a*radius_a)+tmp_local.getY()*tmp_local.getY()/(radius_b*radius_b)+tmp_local.getZ()*tmp_local.getZ()/(radius_c*radius_c) <= 1 ){
        isGhost = 1;  // is Ghost particle
      if(tmp_local.getX()*tmp_local.getX()/(inner_a*inner_a)+tmp_local.getY()*tmp_local.getY()/(inner_b*inner_b)+tmp_local.getZ()*tmp_local.getZ()/(inner_c*inner_c) > 1 || inner_a<=0 || inner_b<=0 || inner_c<=0){
        local_x = (*pt)->globalToLocal( dem::Vec(tmp_x, 0, tmp_z)-pt_position );
        sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity, tmp_x, 0, tmp_z, local_x, 2);  // 2 is ghost SPH particle
        (*pt)->SPHGhostParticleVec.push_back(tmp_pt);   // at current, the scatterSPHParticle is not valide for ghost particles. July 15, 2015

  //          break;  // if this sph point is ghost for dem particle 1, then it cannot be ghost for any others,
        // should avoid the initial overlap of the different dem particles
        }
      }
      } // end dem particle

    if(isGhost==0){  // free/boundary particles
      if(tmp_x<=xmin-spaceInterval+small_value || tmp_x>=xmax-small_value || tmp_z<=zmin-spaceInterval+small_value || tmp_z>=zmax-small_value){  // boundary sph particles
      if(tmp_z>=L+spaceInterval){
          sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity, tmp_x, 0, tmp_z, 3);  // 3 is boundary SPH particle
        allSPHParticleVec.push_back(tmp_pt);
      }
      else{
          sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity*pow(1+SPHInitialDensity*gravAccel*gravScale*(L-tmp_z)/P0, 1.0/gamma), tmp_x, 0, tmp_z, 3);  // 3 is boundary SPH particle
        allSPHParticleVec.push_back(tmp_pt);
      }
      } 
      else if(tmp_x<=L && tmp_z<=L){  // free sph particles
        sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity*pow(1+SPHInitialDensity*gravAccel*gravScale*(L-tmp_z)/P0, 1.0/gamma), tmp_x, 0, tmp_z, 1);  // 1 is free SPH particle
  //      sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity, tmp_x, 0, tmp_z, 1);  // 1 is free SPH particle
          allSPHParticleVec.push_back(tmp_pt);
      }
    }
    }
  }

  for(std::vector<sph::SPHParticle*>::iterator pt=allSPHParticleVec.begin(); pt!=allSPHParticleVec.end(); pt++){
    (*pt)->initial();
  }

} // generateSPHParticle2D

  // here the free particles, ghost particles and boundary particles will all stored in the same vector
  // unlike the implementation in serial version.
  void SmoothParticleHydrodynamics::generateSPHParticle3D(){
  
        if(mpiRank!=0) return;  // make only primary cpu generate particles

  int numLayers = util::getParam<>("numLayers"];
   REAL L = util::getParam<>("waterLength"];
  REAL nSPHPoint = util::getParam<>("nSPHPoint"];
  REAL gamma = util::getParam<>("gamma"];
  REAL P0 = util::getParam<>("P0"];
  REAL gravAccel = util::getParam<>("gravAccel"];
  REAL gravScale = util::getParam<>("gravScale"];
  REAL SPHInitialDensity = util::getParam<>("SPHInitialDensity"];

    spaceInterval = L/(nSPHPoint-1);
  REAL SPHmass = SPHInitialDensity*L*L*L/(nSPHPoint*nSPHPoint*nSPHPoint);

  // get the dimensions of the sph domain
      Vec vmin = allContainer.getMinCorner();
      Vec vmax = allContainer.getMaxCorner();
      REAL xmin = vmin.getX();
      REAL ymin = vmin.getY();
      REAL zmin = vmin.getZ();
      REAL xmax = vmax.getX();
      REAL ymax = vmax.getY();
      REAL zmax = vmax.getZ();

  // create and store PeriParticle objects into periParticleVec
  int isGhost = 0;  // is Ghost particle
  REAL radius_a, radius_b, radius_c;
  REAL inner_a,  inner_b,  inner_c;
    dem::Vec pt_position;
  dem::Vec local_x;  // local position of ghost point in dem particle
  dem::Vec tmp_xyz, tmp_local;
  REAL small_value = 0.01*spaceInterval;
      smoothLength = 1.5*spaceInterval;
  kernelSize = 3*smoothLength;
    for(REAL tmp_y=ymin-spaceInterval*numLayers; tmp_y<ymax+spaceInterval*numLayers; tmp_y=tmp_y+spaceInterval){
  for(REAL tmp_x=xmin-spaceInterval*numLayers; tmp_x<xmax+spaceInterval*numLayers; tmp_x=tmp_x+spaceInterval){
      for(REAL tmp_z=zmin-spaceInterval*numLayers; tmp_z<zmax+spaceInterval*numLayers; tmp_z=tmp_z+spaceInterval){
    isGhost = 0;
    tmp_xyz = dem::Vec(tmp_x, tmp_y, tmp_z);
        for(std::vector<Particle*>::iterator pt=allParticleVec.begin(); pt!=allParticleVec.end(); pt++){
        radius_a = (*pt)->getA(); radius_b = (*pt)->getB(); radius_c = (*pt)->getC();
        inner_a = radius_a-kernelSize; inner_b = radius_b-kernelSize; inner_c = radius_c-kernelSize;
        pt_position = (*pt)->getCurrPos();
        tmp_local = (*pt)->globalToLocal(tmp_xyz-pt_position);
        if( tmp_local.getX()*tmp_local.getX()/(radius_a*radius_a)+tmp_local.getY()*tmp_local.getY()/(radius_b*radius_b)+tmp_local.getZ()*tmp_local.getZ()/(radius_c*radius_c) <= 1 ){
          isGhost = 1;  // is Ghost particle
      if(tmp_local.getX()*tmp_local.getX()/(inner_a*inner_a)+tmp_local.getY()*tmp_local.getY()/(inner_b*inner_b)+tmp_local.getZ()*tmp_local.getZ()/(inner_c*inner_c) > 1 || inner_a<=0 || inner_b<=0 || inner_c<=0){
          local_x = (*pt)->globalToLocal( dem::Vec(tmp_x, tmp_y, tmp_z)-pt_position );
          sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity, tmp_x, tmp_y, tmp_z, local_x, 2);  // 2 is ghost SPH particle
          (*pt)->SPHGhostParticleVec.push_back(tmp_pt);   // at current, the scatterSPHParticle is not valide for ghost particles. July 15, 2015

//          break;  // if this sph point is ghost for dem particle 1, then it cannot be ghost for any others,
        // should avoid the initial overlap of the different dem particles
          }
        }
        } // end dem particle

    if(isGhost==0){  // free/boundary particles
        if(tmp_x<=xmin-spaceInterval+small_value || tmp_x>=xmax-small_value 
        || tmp_y<=ymin-spaceInterval+small_value || tmp_y>=ymax+spaceInterval-small_value
        || tmp_z<=zmin-spaceInterval+small_value || tmp_z>=zmax-small_value){  // boundary sph particles
      if(tmp_z>=L+spaceInterval){
              sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity, tmp_x, tmp_y, tmp_z, 3);  // 3 is boundary SPH particle
          allSPHParticleVec.push_back(tmp_pt);
      }
      else{
              sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity*pow(1+SPHInitialDensity*gravAccel*gravScale*(L-tmp_z)/P0, 1.0/gamma), tmp_x, tmp_y, tmp_z, 3);  // 3 is boundary SPH particle
          allSPHParticleVec.push_back(tmp_pt);
      }
        } 
        else if(tmp_x<=L && tmp_z<=L){  // free sph particles
          sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity*pow(1+SPHInitialDensity*gravAccel*gravScale*(L-tmp_z)/P0, 1.0/gamma), tmp_x, tmp_y, tmp_z, 1);  // 1 is free SPH particle
//      sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity, tmp_x, 0, tmp_z, 1);  // 1 is free SPH particle
              allSPHParticleVec.push_back(tmp_pt);
        }
    }
      }
  }
    }

  for(std::vector<sph::SPHParticle*>::iterator pt=allSPHParticleVec.begin(); pt!=allSPHParticleVec.end(); pt++){
      (*pt)->initial();
  }

  } // generateSPHParticle3D


  // here the free particles and boundary particles will all stored in the same vector while sph ghost particles will be stored in dem particles
  void SmoothParticleHydrodynamics::generateSPHParticleNoBottom3D(){  // this is for drainage problem
  
        if(mpiRank!=0) return;  // make only primary cpu generate particles

  int numLayers = util::getParam<>("numLayers"];
   REAL L = allContainer.getMaxCorner().getX()-allContainer.getMinCorner().getX();  // based on x direction
  REAL gamma = util::getParam<>("gamma"];
  REAL P0 = util::getParam<>("P0"];
  REAL gravAccel = util::getParam<>("gravAccel"];
  REAL gravScale = util::getParam<>("gravScale"];
  REAL SPHInitialDensity = util::getParam<>("SPHInitialDensity"];

    spaceInterval = util::getParam<>("spaceInterval"];
  REAL SPHmass = SPHInitialDensity*spaceInterval*spaceInterval*spaceInterval;

  // get the dimensions of the sph domain
      Vec vmin = allContainer.getMinCorner();
      Vec vmax = allContainer.getMaxCorner();
      REAL xmin = vmin.getX();
      REAL ymin = vmin.getY();
      REAL zmin = vmin.getZ();
      REAL xmax = vmax.getX();
      REAL ymax = vmax.getY();
      REAL zmax = vmax.getZ();

  // create and store PeriParticle objects into periParticleVec
  int isGhost = 0;  // is Ghost particle
  REAL radius_a, radius_b, radius_c;
    dem::Vec pt_position;
  dem::Vec local_x;  // local position of ghost point in dem particle
  dem::Vec tmp_xyz, tmp_local;
  REAL small_value = 0.01*spaceInterval;
    for(REAL tmp_y=ymin-spaceInterval*numLayers; tmp_y<ymax+spaceInterval*numLayers; tmp_y=tmp_y+spaceInterval){
  for(REAL tmp_x=xmin-spaceInterval*numLayers; tmp_x<xmax+spaceInterval*numLayers; tmp_x=tmp_x+spaceInterval){
      for(REAL tmp_z=zmin-spaceInterval*numLayers; tmp_z<zmax+spaceInterval*numLayers; tmp_z=tmp_z+spaceInterval){
    isGhost = 0;
    tmp_xyz = dem::Vec(tmp_x, tmp_y, tmp_z);
        for(std::vector<Particle*>::iterator pt=allParticleVec.begin(); pt!=allParticleVec.end(); pt++){
        radius_a = (*pt)->getA(); radius_b = (*pt)->getB(); radius_c = (*pt)->getC();
        pt_position = (*pt)->getCurrPos();
        tmp_local = (*pt)->globalToLocal(tmp_xyz-pt_position);
        if( tmp_local.getX()*tmp_local.getX()/(radius_a*radius_a)+tmp_local.getY()*tmp_local.getY()/(radius_b*radius_b)+tmp_local.getZ()*tmp_local.getZ()/(radius_c*radius_c) <= 1 ){
          isGhost = 1;  // is Ghost particle
          local_x = (*pt)->globalToLocal( dem::Vec(tmp_x, tmp_y, tmp_z)-pt_position );
          sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity, tmp_x, tmp_y, tmp_z, local_x, 2);  // 2 is ghost SPH particle
          (*pt)->SPHGhostParticleVec.push_back(tmp_pt);   // at current, the scatterSPHParticle is not valide for ghost particles. July 15, 2015

          break;  // if this sph point is ghost for dem particle 1, then it cannot be ghost for any others,
        // should avoid the initial overlap of the different dem particles
        }
        } // end dem particle

    if(isGhost==0){  // free/boundary particles
        if(tmp_x<=xmin-spaceInterval+small_value || tmp_x>=xmax+spaceInterval-small_value 
        || tmp_y<=ymin-spaceInterval+small_value || tmp_y>=ymax+spaceInterval-small_value
    /*|| tmp_z<=zmin-spaceInterval+small_value*/ || tmp_z>=zmax+spaceInterval-small_value){  // boundary sph particles
      if(tmp_z>=zmax+spaceInterval-small_value){
              sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity, tmp_x, tmp_y, tmp_z, 3);  // 3 is boundary SPH particle
          allSPHParticleVec.push_back(tmp_pt);
      }
      else{
              sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity*pow(1+SPHInitialDensity*gravAccel*gravScale*(zmax-tmp_z)/P0, 1.0/gamma), tmp_x, tmp_y, tmp_z, 3);  // 3 is boundary SPH particle
          allSPHParticleVec.push_back(tmp_pt);
      }
        } 
        else if(tmp_z>=zmin){  // free sph particles
          sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity*pow(1+SPHInitialDensity*gravAccel*gravScale*(zmax-tmp_z)/P0, 1.0/gamma), tmp_x, tmp_y, tmp_z, 1);  // 1 is free SPH particle
//      sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity, tmp_x, 0, tmp_z, 1);  // 1 is free SPH particle
              allSPHParticleVec.push_back(tmp_pt);
        }
    }
      }
  }
    }

  for(std::vector<sph::SPHParticle*>::iterator pt=allSPHParticleVec.begin(); pt!=allSPHParticleVec.end(); pt++){
      (*pt)->initial();
  }

  } // generateSPHParticleNoBottom3D


  // here the free particles and boundary particles will all stored in the same vector while sph ghost particles will be stored in dem particles
  void SmoothParticleHydrodynamics::generateSPHParticleMiddleLayers3D(){  // this is for drainage middleLayers
  
        if(mpiRank!=0) return;  // make only primary cpu generate particles

  int numLayers = util::getParam<>("numLayers"];
   REAL L = allContainer.getMaxCorner().getX()-allContainer.getMinCorner().getX();  // based on x direction
  REAL gamma = util::getParam<>("gamma"];
  REAL P0 = util::getParam<>("P0"];
  REAL gravAccel = util::getParam<>("gravAccel"];
  REAL gravScale = util::getParam<>("gravScale"];
  REAL SPHInitialDensity = util::getParam<>("SPHInitialDensity"];

    spaceInterval = util::getParam<>("spaceInterval"];
  REAL SPHmass = SPHInitialDensity*spaceInterval*spaceInterval*spaceInterval;

  // get the dimensions of the sph domain
      Vec vmin = allContainer.getMinCorner();
      Vec vmax = allContainer.getMaxCorner();
      REAL xmin = vmin.getX();
      REAL ymin = vmin.getY();
      REAL zmin = vmin.getZ();
      REAL xmax = vmax.getX();
      REAL ymax = vmax.getY();
      REAL zmax = vmax.getZ();

  REAL waterZmin = util::getParam<>("waterZmin"];
  REAL waterZmax = util::getParam<>("waterZmax"];

  // create and store PeriParticle objects into periParticleVec
  int isGhost = 0;  // is Ghost particle
  REAL radius_a, radius_b, radius_c;
    dem::Vec pt_position;
  dem::Vec local_x;  // local position of ghost point in dem particle
  dem::Vec tmp_xyz, tmp_local;
  REAL small_value = 0.01*spaceInterval;
    for(REAL tmp_y=ymin-spaceInterval*numLayers; tmp_y<ymax+spaceInterval*numLayers; tmp_y=tmp_y+spaceInterval){
  for(REAL tmp_x=xmin-spaceInterval*numLayers; tmp_x<xmax+spaceInterval*numLayers; tmp_x=tmp_x+spaceInterval){
      for(REAL tmp_z=zmin-spaceInterval*numLayers; tmp_z<zmax+spaceInterval*numLayers; tmp_z=tmp_z+spaceInterval){
    isGhost = 0;
    tmp_xyz = dem::Vec(tmp_x, tmp_y, tmp_z);
        for(std::vector<Particle*>::iterator pt=allParticleVec.begin(); pt!=allParticleVec.end(); pt++){
        radius_a = (*pt)->getA(); radius_b = (*pt)->getB(); radius_c = (*pt)->getC();
        pt_position = (*pt)->getCurrPos();
        tmp_local = (*pt)->globalToLocal(tmp_xyz-pt_position);
        if( tmp_local.getX()*tmp_local.getX()/(radius_a*radius_a)+tmp_local.getY()*tmp_local.getY()/(radius_b*radius_b)+tmp_local.getZ()*tmp_local.getZ()/(radius_c*radius_c) <= 1 ){
          isGhost = 1;  // is Ghost particle
          local_x = (*pt)->globalToLocal( dem::Vec(tmp_x, tmp_y, tmp_z)-pt_position );
          sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity, tmp_x, tmp_y, tmp_z, local_x, 2);  // 2 is ghost SPH particle
          (*pt)->SPHGhostParticleVec.push_back(tmp_pt);   // at current, the scatterSPHParticle is not valide for ghost particles. July 15, 2015

          break;  // if this sph point is ghost for dem particle 1, then it cannot be ghost for any others,
        // should avoid the initial overlap of the different dem particles
        }
        } // end dem particle

    if(isGhost==0){  // free/boundary particles
        if(tmp_x<=xmin-spaceInterval+small_value || tmp_x>=xmax+spaceInterval-small_value 
        || tmp_y<=ymin-spaceInterval+small_value || tmp_y>=ymax+spaceInterval-small_value
        || tmp_z<=zmin-spaceInterval+small_value /*|| tmp_z>=zmax+spaceInterval-small_value*/){  // boundary sph particles, no top boundary
          sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity, tmp_x, tmp_y, tmp_z, 3);  // 3 is boundary SPH particle
      allSPHParticleVec.push_back(tmp_pt);
        } 
        else if(tmp_z>=waterZmin && tmp_z<=waterZmax){  // free sph particles
          sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity*pow(1+SPHInitialDensity*gravAccel*gravScale*(waterZmax-tmp_z)/P0, 1.0/gamma), tmp_x, tmp_y, tmp_z, 1);  // 1 is free SPH particle
//      sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity, tmp_x, 0, tmp_z, 1);  // 1 is free SPH particle
              allSPHParticleVec.push_back(tmp_pt);
        }
    }
      }
  }
    }

  for(std::vector<sph::SPHParticle*>::iterator pt=allSPHParticleVec.begin(); pt!=allSPHParticleVec.end(); pt++){
      (*pt)->initial();
  }

  } // generateSPHParticleNoBottom3D

  void SmoothParticleHydrodynamics::findSPHParticleInRectangle(const Rectangle &container,
        const std::vector<sph::SPHParticle*> &inputParticle,
        std::vector<sph::SPHParticle*> &foundParticle) {
    foundParticle.reserve(inputParticle.size());
    Vec  v1 = container.getMinCorner();
    Vec  v2 = container.getMaxCorner();
    REAL x1 = v1.getX();
    REAL y1 = v1.getY();
    REAL z1 = v1.getZ();
    REAL x2 = v2.getX();
    REAL y2 = v2.getY();
    REAL z2 = v2.getZ();
    for (std::size_t pt = 0; pt < inputParticle.size(); ++pt) {
      Vec center = inputParticle[pt]->getCurrPosition();
      // it is critical to use EPS
      if (center.getX() - x1 >= -EPS && center.getX() - x2 < -EPS &&
    center.getY() - y1 >= -EPS && center.getY() - y2 < -EPS &&
    center.getZ() - z1 >= -EPS && center.getZ() - z2 < -EPS)
  foundParticle.push_back(inputParticle[pt]);
    }
    std::vector<sph::SPHParticle*>(foundParticle).swap(foundParticle);
  }

  void SmoothParticleHydrodynamics::removeParticleOutRectangle() {
    Vec  v1 = container.getMinCorner();
    Vec  v2 = container.getMaxCorner();
    REAL x1 = v1.getX();
    REAL y1 = v1.getY();
    REAL z1 = v1.getZ();
    REAL x2 = v2.getX();
    REAL y2 = v2.getY();
    REAL z2 = v2.getZ();

    std::vector<Particle*>::iterator itr;
    Vec center;
    //std::size_t flag = 0;

    for (itr = particleVec.begin(); itr != particleVec.end(); ) {
      center=(*itr)->getCurrPos();
      // it is critical to use EPS
      if ( !(center.getX() - x1 >= -EPS && center.getX() - x2 < -EPS &&
       center.getY() - y1 >= -EPS && center.getY() - y2 < -EPS &&
       center.getZ() - z1 >= -EPS && center.getZ() - z2 < -EPS) )
  {
    /*
      debugInf << "iter=" << std::setw(8) << iteration << " rank=" << std::setw(2) << mpiRank
      << " removed=" << std::setw(3) << (*itr)->getId();  
      flag = 1;
    */
    // this is important to free the memory of these sph ghost particles
    for(std::vector<sph::SPHParticle*>::iterator st=(*itr)->SPHGhostParticleVec.begin(); st!=(*itr)->SPHGhostParticleVec.end(); st++){
       delete (*st);  
    }
    (*itr)->SPHGhostParticleVec.clear();

    delete (*itr); // release memory
    itr = particleVec.erase(itr); 
  }
      else
  ++itr;
    }
    /*
      if (flag == 1) {
      debugInf << " now " << particleVec.size() << ": ";
      for (std::vector<Particle*>::const_iterator it = particleVec.begin(); it != particleVec.end(); ++it)
      debugInf << std::setw(3) << (*it)->getId();
      debugInf << std::endl;
      }
    */

  }

  void SmoothParticleHydrodynamics::removeSPHParticleOutRectangle() {
    Vec  v1 = container.getMinCorner();
    Vec  v2 = container.getMaxCorner();
    REAL x1 = v1.getX();
    REAL y1 = v1.getY();
    REAL z1 = v1.getZ();
    REAL x2 = v2.getX();
    REAL y2 = v2.getY();
    REAL z2 = v2.getZ();

    std::vector<sph::SPHParticle*>::iterator itr;
    Vec center;
    //std::size_t flag = 0;

    // this for loop may be replaced by for_each, remove, erase and swap (shrink to fit)
    // see item 14 in "effective STL"
    for (itr = SPHParticleVec.begin(); itr != SPHParticleVec.end(); ) {
      center=(*itr)->getCurrPosition();
      // it is critical to use EPS
      if ( !(center.getX() - x1 >= -EPS && center.getX() - x2 < -EPS &&
       center.getY() - y1 >= -EPS && center.getY() - y2 < -EPS &&
       center.getZ() - z1 >= -EPS && center.getZ() - z2 < -EPS) )
  {
    /*
      debugInf << "iter=" << std::setw(8) << iteration << " rank=" << std::setw(2) << mpiRank
      << " removed=" << std::setw(3) << (*itr)->getId();  
      flag = 1;
    */
    delete (*itr); // release memory
    itr = SPHParticleVec.erase(itr); 
  }
      else
  ++itr;
    }
    /*
      if (flag == 1) {
      debugInf << " now " << particleVec.size() << ": ";
      for (std::vector<Particle*>::const_iterator it = particleVec.begin(); it != particleVec.end(); ++it)
      debugInf << std::setw(3) << (*it)->getId();
      debugInf << std::endl;
      }
    */

  }

  // this is to scatter the dem and sph particle
  // two point: (1) the sph ghost particles will be partitioned with dem particles together. There is no explicite partition for sph ghost particles.
  //        After the partition of the dem particles, the sph ghost particles will be paritioned as a class member of dem particle.
  //        Before partition, SPHParticle.demParticle points to NULL; after receive, SPHParticle.demParticle can point to their dem particle. July 27, 2015
  //            (2) block partitions for dem domain and sph domain have to be the same, 
  //        otherwise some dem particles and their nearby sph particles will belong to different block
  //        (grid is expanded to include the boundary sph particles)
  //    (3) this partition method here is not very suitable for the free surface flow problem, such as bursting dam problem
  //        since the total domain is divided, while there are lots of voids in the domain. (partition only free sph particles will be better)
  //        But our goal is to simulate the porous media in triaxial, insotropic or SHPB simulations, particles are filled in the container.
  void SmoothParticleHydrodynamics::scatterDEMSPHParticle() {
    // partition particles and send to each process
    if (mpiRank == 0) { // process 0
      int numLayers = util::getParam<>("numLayers"];
      setGrid(Rectangle(allContainer.getMinCorner().getX() - spaceInterval*numLayers,
      allContainer.getMinCorner().getY() - spaceInterval*numLayers,
      allContainer.getMinCorner().getZ() - spaceInterval*numLayers,
      allContainer.getMaxCorner().getX() + spaceInterval*numLayers,
      allContainer.getMaxCorner().getY() + spaceInterval*numLayers,
      allContainer.getMaxCorner().getZ() + spaceInterval*numLayers ));
    
      Vec v1 = grid.getMinCorner();
      Vec v2 = grid.getMaxCorner();
      Vec vspan = v2 - v1;

      boost::mpi::request *reqs = new boost::mpi::request [mpiSize - 1];
      std::vector<Particle*> tmpParticleVec;
      for (int iRank = mpiSize - 1; iRank >= 0; --iRank) {
  tmpParticleVec.clear(); // do not release memory!
  int ndim = 3;
  int coords[3];
  MPI_Cart_coords(cartComm, iRank, ndim, coords);
  Rectangle container(v1.getX() + vspan.getX() / mpiProcX * coords[0],
          v1.getY() + vspan.getY() / mpiProcY * coords[1],
          v1.getZ() + vspan.getZ() / mpiProcZ * coords[2],
          v1.getX() + vspan.getX() / mpiProcX * (coords[0] + 1),
          v1.getY() + vspan.getY() / mpiProcY * (coords[1] + 1),
          v1.getZ() + vspan.getZ() / mpiProcZ * (coords[2] + 1));
  findParticleInRectangle(container, allParticleVec, tmpParticleVec);
  if (iRank != 0)
    reqs[iRank - 1] = boostWorld.isend(iRank, mpiTag, tmpParticleVec); // non-blocking send
    // before send, the SPHParticle.demParticle == NULL, since NULL is assigned when SPHParticle is created
  if (iRank == 0) {
    particleVec.resize(tmpParticleVec.size());
    for (int i = 0; i < particleVec.size(); ++i){
      // default synthesized copy constructor
      particleVec[i] = new Particle(*tmpParticleVec[i]);  // at this point, particleVec and allParticleVec are pointing to the same SPHGhoastParticle
      particleVec[i]->SPHGhostParticleVec.clear();  // at this point, particleVec is pointing to nothing
      for(std::vector<sph::SPHParticle*>::iterator st=tmpParticleVec[i]->SPHGhostParticleVec.begin(); st!=tmpParticleVec[i]->SPHGhostParticleVec.end(); st++){
      sph::SPHParticle* tmp_sph = new sph::SPHParticle(**st);  // create a new SPHGhost particle, which is the same as the one in allParticleVec
      particleVec[i]->SPHGhostParticleVec.push_back(tmp_sph);  // now particleVec points to the new SPHGhostParticle
      }
    }
  } // now particleVec do not share memeory with allParticleVec
      }
      boost::mpi::wait_all(reqs, reqs + mpiSize - 1); // for non-blocking send
      delete [] reqs;

    } else { // other processes except 0
      boostWorld.recv(0, mpiTag, particleVec);
    }

    // this is important to pointer the demParticle in SPHParticle to the correct dem particle
    for(std::vector<Particle*>::iterator it=particleVec.begin(); it!=particleVec.end(); it++){
  (*it)->setDemParticleInSPHParticle();  // set SPHParticle.demParticle to (*it)
    }

    // content of allParticleVec may need to be printed, so do not clear it. 
    //if (mpiRank == 0) releaseGatheredParticle();


    ///////////////////////////////////////////////////////////////////////////////
    // partition SPH particles (free and boundary) and send to each process
    if (mpiRank == 0) { // process 0
      int numLayers = util::getParam<>("numLayers"];
      setGrid(Rectangle(allContainer.getMinCorner().getX() - spaceInterval*numLayers,
      allContainer.getMinCorner().getY() - spaceInterval*numLayers,
      allContainer.getMinCorner().getZ() - spaceInterval*numLayers,
      allContainer.getMaxCorner().getX() + spaceInterval*numLayers,
      allContainer.getMaxCorner().getY() + spaceInterval*numLayers,
      allContainer.getMaxCorner().getZ() + spaceInterval*numLayers ));
    
      Vec v1 = grid.getMinCorner();
      Vec v2 = grid.getMaxCorner();
      Vec vspan = v2 - v1;

      boost::mpi::request *reqs = new boost::mpi::request [mpiSize - 1];
      std::vector<sph::SPHParticle*> tmpSPHParticleVec;
      for (int iRank = mpiSize - 1; iRank >= 0; --iRank) {
  tmpSPHParticleVec.clear(); // do not release memory!
  int ndim = 3;
  int coords[3];
  MPI_Cart_coords(cartComm, iRank, ndim, coords);
  Rectangle container(v1.getX() + vspan.getX() / mpiProcX * coords[0],
          v1.getY() + vspan.getY() / mpiProcY * coords[1],
          v1.getZ() + vspan.getZ() / mpiProcZ * coords[2],
          v1.getX() + vspan.getX() / mpiProcX * (coords[0] + 1),
          v1.getY() + vspan.getY() / mpiProcY * (coords[1] + 1),
          v1.getZ() + vspan.getZ() / mpiProcZ * (coords[2] + 1));
  findSPHParticleInRectangle(container, allSPHParticleVec, tmpSPHParticleVec);
  if (iRank != 0)
    reqs[iRank - 1] = boostWorld.isend(iRank, mpiTag, tmpSPHParticleVec); // non-blocking send
  if (iRank == 0) {
    SPHParticleVec.resize(tmpSPHParticleVec.size());
    for (int i = 0; i < SPHParticleVec.size(); ++i)
      SPHParticleVec[i] = new sph::SPHParticle(*tmpSPHParticleVec[i]); // default synthesized copy constructor
  } // now particleVec do not share memeory with allParticleVec
      }
      boost::mpi::wait_all(reqs, reqs + mpiSize - 1); // for non-blocking send
      delete [] reqs;

    } else { // other processes except 0
      boostWorld.recv(0, mpiTag, SPHParticleVec);
    }


    // broadcast necessary info
    broadcast(boostWorld, gradation, 0);
    broadcast(boostWorld, boundaryVec, 0);
    broadcast(boostWorld, allContainer, 0);
    broadcast(boostWorld, grid, 0);
  } // scatterDEMSPHParticle


  void SmoothParticleHydrodynamics::scatterDEMSPHParticleCopyDEM() {
    // partition particles and send to each process
    if (mpiRank == 0) { // process 0
      int numLayers = util::getParam<>("numLayers"];
      setGrid(Rectangle(allContainer.getMinCorner().getX() - spaceInterval*numLayers,
      allContainer.getMinCorner().getY() - spaceInterval*numLayers,
      allContainer.getMinCorner().getZ() - spaceInterval*numLayers,
      allContainer.getMaxCorner().getX() + spaceInterval*numLayers,
      allContainer.getMaxCorner().getY() + spaceInterval*numLayers,
      allContainer.getMaxCorner().getZ() + spaceInterval*numLayers ));
    
      Vec v1 = grid.getMinCorner();
      Vec v2 = grid.getMaxCorner();
      Vec vspan = v2 - v1;

      boost::mpi::request *reqs = new boost::mpi::request [mpiSize - 1];
      std::vector<Particle*> tmpParticleVec;
      for (int iRank = mpiSize - 1; iRank >= 0; --iRank) {
  tmpParticleVec.clear(); // do not release memory!
  int ndim = 3;
  int coords[3];
  MPI_Cart_coords(cartComm, iRank, ndim, coords);
  Rectangle container(v1.getX(),
          v1.getY(),
          v1.getZ(),
          v2.getX(),
          v2.getY(),
          v2.getZ() );
  findParticleInRectangle(container, allParticleVec, tmpParticleVec);
  if (iRank != 0)
    reqs[iRank - 1] = boostWorld.isend(iRank, mpiTag, tmpParticleVec); // non-blocking send
    // before send, the SPHParticle.demParticle == NULL, since NULL is assigned when SPHParticle is created
  if (iRank == 0) {
    particleVec.resize(tmpParticleVec.size());
    for (int i = 0; i < particleVec.size(); ++i){
      // default synthesized copy constructor
      particleVec[i] = new Particle(*tmpParticleVec[i]);  // at this point, particleVec and allParticleVec are pointing to the same SPHGhoastParticle
      particleVec[i]->SPHGhostParticleVec.clear();  // at this point, particleVec is pointing to nothing
      for(std::vector<sph::SPHParticle*>::iterator st=tmpParticleVec[i]->SPHGhostParticleVec.begin(); st!=tmpParticleVec[i]->SPHGhostParticleVec.end(); st++){
      sph::SPHParticle* tmp_sph = new sph::SPHParticle(**st);  // create a new SPHGhost particle, which is the same as the one in allParticleVec
      particleVec[i]->SPHGhostParticleVec.push_back(tmp_sph);  // now particleVec points to the new SPHGhostParticle
      }
    }
  } // now particleVec do not share memeory with allParticleVec
      }
      boost::mpi::wait_all(reqs, reqs + mpiSize - 1); // for non-blocking send
      delete [] reqs;

    } else { // other processes except 0
      boostWorld.recv(0, mpiTag, particleVec);
    }

    // this is important to pointer the demParticle in SPHParticle to the correct dem particle
    for(std::vector<Particle*>::iterator it=particleVec.begin(); it!=particleVec.end(); it++){
  (*it)->setDemParticleInSPHParticle();  // set SPHParticle.demParticle to (*it)
    }

    // content of allParticleVec may need to be printed, so do not clear it. 
    //if (mpiRank == 0) releaseGatheredParticle();


    ///////////////////////////////////////////////////////////////////////////////
    // partition SPH particles (free and boundary) and send to each process
    if (mpiRank == 0) { // process 0
      int numLayers = util::getParam<>("numLayers"];
      setGrid(Rectangle(allContainer.getMinCorner().getX() - spaceInterval*numLayers,
      allContainer.getMinCorner().getY() - spaceInterval*numLayers,
      allContainer.getMinCorner().getZ() - spaceInterval*numLayers,
      allContainer.getMaxCorner().getX() + spaceInterval*numLayers,
      allContainer.getMaxCorner().getY() + spaceInterval*numLayers,
      allContainer.getMaxCorner().getZ() + spaceInterval*numLayers ));
    
      Vec v1 = grid.getMinCorner();
      Vec v2 = grid.getMaxCorner();
      Vec vspan = v2 - v1;

      boost::mpi::request *reqs = new boost::mpi::request [mpiSize - 1];
      std::vector<sph::SPHParticle*> tmpSPHParticleVec;
      for (int iRank = mpiSize - 1; iRank >= 0; --iRank) {
  tmpSPHParticleVec.clear(); // do not release memory!
  int ndim = 3;
  int coords[3];
  MPI_Cart_coords(cartComm, iRank, ndim, coords);
  Rectangle container(v1.getX() + vspan.getX() / mpiProcX * coords[0],
          v1.getY() + vspan.getY() / mpiProcY * coords[1],
          v1.getZ() + vspan.getZ() / mpiProcZ * coords[2],
          v1.getX() + vspan.getX() / mpiProcX * (coords[0] + 1),
          v1.getY() + vspan.getY() / mpiProcY * (coords[1] + 1),
          v1.getZ() + vspan.getZ() / mpiProcZ * (coords[2] + 1));
  findSPHParticleInRectangle(container, allSPHParticleVec, tmpSPHParticleVec);
  if (iRank != 0)
    reqs[iRank - 1] = boostWorld.isend(iRank, mpiTag, tmpSPHParticleVec); // non-blocking send
  if (iRank == 0) {
    SPHParticleVec.resize(tmpSPHParticleVec.size());
    for (int i = 0; i < SPHParticleVec.size(); ++i)
      SPHParticleVec[i] = new sph::SPHParticle(*tmpSPHParticleVec[i]); // default synthesized copy constructor
  } // now particleVec do not share memeory with allParticleVec
      }
      boost::mpi::wait_all(reqs, reqs + mpiSize - 1); // for non-blocking send
      delete [] reqs;

    } else { // other processes except 0
      boostWorld.recv(0, mpiTag, SPHParticleVec);
    }


    // broadcast necessary info
    broadcast(boostWorld, gradation, 0);
    broadcast(boostWorld, boundaryVec, 0);
    broadcast(boostWorld, allContainer, 0);
    broadcast(boostWorld, grid, 0);
  } // scatterDEMSPHParticleCopyDEM

  void SmoothParticleHydrodynamics::commuParticle()   // the communication of sph ghost particles are implemented here, July 27, 2015
  {
    // determine container of each process
    Vec v1 = grid.getMinCorner();
    Vec v2 = grid.getMaxCorner();
    Vec vspan = v2 - v1;
    container = Rectangle(v1.getX() + vspan.getX() / mpiProcX * mpiCoords[0],
        v1.getY() + vspan.getY() / mpiProcY * mpiCoords[1],
        v1.getZ() + vspan.getZ() / mpiProcZ * mpiCoords[2],
        v1.getX() + vspan.getX() / mpiProcX * (mpiCoords[0] + 1),
        v1.getY() + vspan.getY() / mpiProcY * (mpiCoords[1] + 1),
        v1.getZ() + vspan.getZ() / mpiProcZ * (mpiCoords[2] + 1));

    // if found, communicate with neighboring blocks
    std::vector<Particle*> particleX1, particleX2;
    std::vector<Particle*> particleY1, particleY2;
    std::vector<Particle*> particleZ1, particleZ2;
    std::vector<Particle*> particleX1Y1, particleX1Y2, particleX1Z1, particleX1Z2; 
    std::vector<Particle*> particleX2Y1, particleX2Y2, particleX2Z1, particleX2Z2; 
    std::vector<Particle*> particleY1Z1, particleY1Z2, particleY2Z1, particleY2Z2; 
    std::vector<Particle*> particleX1Y1Z1, particleX1Y1Z2, particleX1Y2Z1, particleX1Y2Z2; 
    std::vector<Particle*> particleX2Y1Z1, particleX2Y1Z2, particleX2Y2Z1, particleX2Y2Z2; 
    boost::mpi::request reqX1[2], reqX2[2];
    boost::mpi::request reqY1[2], reqY2[2];
    boost::mpi::request reqZ1[2], reqZ2[2];
    boost::mpi::request reqX1Y1[2], reqX1Y2[2], reqX1Z1[2], reqX1Z2[2];
    boost::mpi::request reqX2Y1[2], reqX2Y2[2], reqX2Z1[2], reqX2Z2[2];
    boost::mpi::request reqY1Z1[2], reqY1Z2[2], reqY2Z1[2], reqY2Z2[2];
    boost::mpi::request reqX1Y1Z1[2], reqX1Y1Z2[2], reqX1Y2Z1[2], reqX1Y2Z2[2];
    boost::mpi::request reqX2Y1Z1[2], reqX2Y1Z2[2], reqX2Y2Z1[2], reqX2Y2Z2[2];
    v1 = container.getMinCorner(); // redefine v1, v2 in terms of process
    v2 = container.getMaxCorner();   
    //debugInf << "rank=" << mpiRank << ' ' << v1.getX() << ' ' << v1.getY() << ' ' << v1.getZ() << ' '  << v2.getX() << ' ' << v2.getY() << ' ' << v2.getZ() << std::endl;
    REAL cellSize = gradation.getPtclMaxRadius() * 2;
    // 6 surfaces
    if (rankX1 >= 0) { // surface x1
      Rectangle containerX1(v1.getX(), v1.getY(), v1.getZ(), 
          v1.getX() + cellSize, v2.getY(), v2.getZ());
      findParticleInRectangle(containerX1, particleVec, particleX1);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX1.begin(); it!=particleX1.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();  // now the SPHGhostParticles in part of particleVec are pointing to NULL, needed to be back after send
      reqX1[0] = boostWorld.isend(rankX1, mpiTag,  particleX1);
      reqX1[1] = boostWorld.irecv(rankX1, mpiTag, rParticleX1);
    }
    if (rankX2 >= 0) { // surface x2
      Rectangle containerX2(v2.getX() - cellSize, v1.getY(), v1.getZ(),
          v2.getX(), v2.getY(), v2.getZ());
      findParticleInRectangle(containerX2, particleVec, particleX2);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX2.begin(); it!=particleX2.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX2[0] = boostWorld.isend(rankX2, mpiTag,  particleX2);
      reqX2[1] = boostWorld.irecv(rankX2, mpiTag, rParticleX2);
    }
    if (rankY1 >= 0) {  // surface y1
      Rectangle containerY1(v1.getX(), v1.getY(), v1.getZ(), 
          v2.getX(), v1.getY() + cellSize, v2.getZ());
      findParticleInRectangle(containerY1, particleVec, particleY1);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleY1.begin(); it!=particleY1.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqY1[0] = boostWorld.isend(rankY1, mpiTag,  particleY1);
      reqY1[1] = boostWorld.irecv(rankY1, mpiTag, rParticleY1);
    }
    if (rankY2 >= 0) {  // surface y2
      Rectangle containerY2(v1.getX(), v2.getY() - cellSize, v1.getZ(),
          v2.getX(), v2.getY(), v2.getZ());
      findParticleInRectangle(containerY2, particleVec, particleY2);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleY2.begin(); it!=particleY2.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqY2[0] = boostWorld.isend(rankY2, mpiTag,  particleY2);
      reqY2[1] = boostWorld.irecv(rankY2, mpiTag, rParticleY2);
    }
    if (rankZ1 >= 0) {  // surface z1
      Rectangle containerZ1(v1.getX(), v1.getY(), v1.getZ(),
          v2.getX(), v2.getY(), v1.getZ() + cellSize);
      findParticleInRectangle(containerZ1, particleVec, particleZ1);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleZ1.begin(); it!=particleZ1.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqZ1[0] = boostWorld.isend(rankZ1, mpiTag,  particleZ1);
      reqZ1[1] = boostWorld.irecv(rankZ1, mpiTag, rParticleZ1);
    }
    if (rankZ2 >= 0) {  // surface z2
      Rectangle containerZ2(v1.getX(), v1.getY(), v2.getZ() - cellSize,
          v2.getX(), v2.getY(), v2.getZ());
      findParticleInRectangle(containerZ2, particleVec, particleZ2);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleZ2.begin(); it!=particleZ2.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqZ2[0] = boostWorld.isend(rankZ2, mpiTag,  particleZ2);
      reqZ2[1] = boostWorld.irecv(rankZ2, mpiTag, rParticleZ2);
    }
    // 12 edges
    if (rankX1Y1 >= 0) { // edge x1y1
      Rectangle containerX1Y1(v1.getX(), v1.getY(), v1.getZ(),
            v1.getX() + cellSize, v1.getY() + cellSize, v2.getZ());
      findParticleInRectangle(containerX1Y1, particleVec, particleX1Y1);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX1Y1.begin(); it!=particleX1Y1.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX1Y1[0] = boostWorld.isend(rankX1Y1, mpiTag,  particleX1Y1);
      reqX1Y1[1] = boostWorld.irecv(rankX1Y1, mpiTag, rParticleX1Y1);
    }
    if (rankX1Y2 >= 0) { // edge x1y2
      Rectangle containerX1Y2(v1.getX(), v2.getY() - cellSize, v1.getZ(),
            v1.getX() + cellSize, v2.getY(), v2.getZ());
      findParticleInRectangle(containerX1Y2, particleVec, particleX1Y2);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX1Y2.begin(); it!=particleX1Y2.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX1Y2[0] = boostWorld.isend(rankX1Y2, mpiTag,  particleX1Y2);
      reqX1Y2[1] = boostWorld.irecv(rankX1Y2, mpiTag, rParticleX1Y2);
    }
    if (rankX1Z1 >= 0) { // edge x1z1
      Rectangle containerX1Z1(v1.getX(), v1.getY(), v1.getZ(),
            v1.getX() + cellSize, v2.getY(), v1.getZ() + cellSize);
      findParticleInRectangle(containerX1Z1, particleVec, particleX1Z1);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX1Z1.begin(); it!=particleX1Z1.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX1Z1[0] = boostWorld.isend(rankX1Z1, mpiTag,  particleX1Z1);
      reqX1Z1[1] = boostWorld.irecv(rankX1Z1, mpiTag, rParticleX1Z1);
    }
    if (rankX1Z2 >= 0) { // edge x1z2
      Rectangle containerX1Z2(v1.getX(), v1.getY(), v2.getZ() - cellSize,
            v1.getX() + cellSize, v2.getY(), v2.getZ());
      findParticleInRectangle(containerX1Z2, particleVec, particleX1Z2);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX1Z2.begin(); it!=particleX1Z2.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX1Z2[0] = boostWorld.isend(rankX1Z2, mpiTag,  particleX1Z2);
      reqX1Z2[1] = boostWorld.irecv(rankX1Z2, mpiTag, rParticleX1Z2);
    }
    if (rankX2Y1 >= 0) { // edge x2y1
      Rectangle containerX2Y1(v2.getX() - cellSize, v1.getY(), v1.getZ(),
            v2.getX(), v1.getY() + cellSize, v2.getZ());
      findParticleInRectangle(containerX2Y1, particleVec, particleX2Y1);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX2Y1.begin(); it!=particleX2Y1.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX2Y1[0] = boostWorld.isend(rankX2Y1, mpiTag,  particleX2Y1);
      reqX2Y1[1] = boostWorld.irecv(rankX2Y1, mpiTag, rParticleX2Y1);
    }
    if (rankX2Y2 >= 0) { // edge x2y2
      Rectangle containerX2Y2(v2.getX() - cellSize, v2.getY() - cellSize, v1.getZ(),
            v2.getX(), v2.getY(), v2.getZ());
      findParticleInRectangle(containerX2Y2, particleVec, particleX2Y2);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX2Y2.begin(); it!=particleX2Y2.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX2Y2[0] = boostWorld.isend(rankX2Y2, mpiTag,  particleX2Y2);
      reqX2Y2[1] = boostWorld.irecv(rankX2Y2, mpiTag, rParticleX2Y2);
    }
    if (rankX2Z1 >= 0) { // edge x2z1
      Rectangle containerX2Z1(v2.getX() - cellSize, v1.getY(), v1.getZ(),
            v2.getX(), v2.getY(), v1.getZ() + cellSize);
      findParticleInRectangle(containerX2Z1, particleVec, particleX2Z1);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX2Z1.begin(); it!=particleX2Z1.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX2Z1[0] = boostWorld.isend(rankX2Z1, mpiTag,  particleX2Z1);
      reqX2Z1[1] = boostWorld.irecv(rankX2Z1, mpiTag, rParticleX2Z1);
    }
    if (rankX2Z2 >= 0) { // edge x2z2
      Rectangle containerX2Z2(v2.getX() - cellSize, v1.getY(), v2.getZ() - cellSize,
            v2.getX(), v2.getY(), v2.getZ());
      findParticleInRectangle(containerX2Z2, particleVec, particleX2Z2);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX2Z2.begin(); it!=particleX2Z2.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX2Z2[0] = boostWorld.isend(rankX2Z2, mpiTag,  particleX2Z2);
      reqX2Z2[1] = boostWorld.irecv(rankX2Z2, mpiTag, rParticleX2Z2);
    }
    if (rankY1Z1 >= 0) { // edge y1z1
      Rectangle containerY1Z1(v1.getX(), v1.getY(), v1.getZ(),
            v2.getX(), v1.getY() + cellSize, v1.getZ() + cellSize);
      findParticleInRectangle(containerY1Z1, particleVec, particleY1Z1);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleY1Z1.begin(); it!=particleY1Z1.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqY1Z1[0] = boostWorld.isend(rankY1Z1, mpiTag,  particleY1Z1);
      reqY1Z1[1] = boostWorld.irecv(rankY1Z1, mpiTag, rParticleY1Z1);
    }
    if (rankY1Z2 >= 0) { // edge y1z2
      Rectangle containerY1Z2(v1.getX(), v1.getY(), v2.getZ() - cellSize,
            v2.getX(), v1.getY() + cellSize, v2.getZ());
      findParticleInRectangle(containerY1Z2, particleVec, particleY1Z2);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleY1Z2.begin(); it!=particleY1Z2.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqY1Z2[0] = boostWorld.isend(rankY1Z2, mpiTag,  particleY1Z2);
      reqY1Z2[1] = boostWorld.irecv(rankY1Z2, mpiTag, rParticleY1Z2);
    }
    if (rankY2Z1 >= 0) { // edge y2z1
      Rectangle containerY2Z1(v1.getX(), v2.getY() - cellSize, v1.getZ(),
            v2.getX(), v2.getY(), v1.getZ() + cellSize);
      findParticleInRectangle(containerY2Z1, particleVec, particleY2Z1);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleY2Z1.begin(); it!=particleY2Z1.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqY2Z1[0] = boostWorld.isend(rankY2Z1, mpiTag,  particleY2Z1);
      reqY2Z1[1] = boostWorld.irecv(rankY2Z1, mpiTag, rParticleY2Z1);
    }
    if (rankY2Z2 >= 0) { // edge y2z2
      Rectangle containerY2Z2(v1.getX(), v2.getY() - cellSize, v2.getZ() - cellSize,
            v2.getX(), v2.getY(), v2.getZ());
      findParticleInRectangle(containerY2Z2, particleVec, particleY2Z2);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleY2Z2.begin(); it!=particleY2Z2.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqY2Z2[0] = boostWorld.isend(rankY2Z2, mpiTag,  particleY2Z2);
      reqY2Z2[1] = boostWorld.irecv(rankY2Z2, mpiTag, rParticleY2Z2);
    }
    // 8 vertices
    if (rankX1Y1Z1 >= 0) { // edge x1y1z1
      Rectangle containerX1Y1Z1(v1.getX(), v1.getY(), v1.getZ(),
        v1.getX() + cellSize, v1.getY() + cellSize, v1.getZ() + cellSize);
      findParticleInRectangle(containerX1Y1Z1, particleVec, particleX1Y1Z1);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX1Y1Z1.begin(); it!=particleX1Y1Z1.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX1Y1Z1[0] = boostWorld.isend(rankX1Y1Z1, mpiTag,  particleX1Y1Z1);
      reqX1Y1Z1[1] = boostWorld.irecv(rankX1Y1Z1, mpiTag, rParticleX1Y1Z1);
    }
    if (rankX1Y1Z2 >= 0) { // edge x1y1z2
      Rectangle containerX1Y1Z2(v1.getX(), v1.getY(), v2.getZ() - cellSize,
        v1.getX() + cellSize, v1.getY() + cellSize, v2.getZ());
      findParticleInRectangle(containerX1Y1Z2, particleVec, particleX1Y1Z2);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX1Y1Z2.begin(); it!=particleX1Y1Z2.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX1Y1Z2[0] = boostWorld.isend(rankX1Y1Z2, mpiTag,  particleX1Y1Z2);
      reqX1Y1Z2[1] = boostWorld.irecv(rankX1Y1Z2, mpiTag, rParticleX1Y1Z2);
    }
    if (rankX1Y2Z1 >= 0) { // edge x1y2z1
      Rectangle containerX1Y2Z1(v1.getX(), v2.getY() - cellSize, v1.getZ(),
        v1.getX() + cellSize, v2.getY(), v1.getZ() + cellSize);
      findParticleInRectangle(containerX1Y2Z1, particleVec, particleX1Y2Z1);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX1Y2Z1.begin(); it!=particleX1Y2Z1.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX1Y2Z1[0] = boostWorld.isend(rankX1Y2Z1, mpiTag,  particleX1Y2Z1);
      reqX1Y2Z1[1] = boostWorld.irecv(rankX1Y2Z1, mpiTag, rParticleX1Y2Z1);
    }
    if (rankX1Y2Z2 >= 0) { // edge x1y2z2
      Rectangle containerX1Y2Z2(v1.getX(), v2.getY() - cellSize, v2.getZ() - cellSize,
        v1.getX() + cellSize, v2.getY() + cellSize, v2.getZ());
      findParticleInRectangle(containerX1Y2Z2, particleVec, particleX1Y2Z2);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX1Y2Z2.begin(); it!=particleX1Y2Z2.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX1Y2Z2[0] = boostWorld.isend(rankX1Y2Z2, mpiTag,  particleX1Y2Z2);
      reqX1Y2Z2[1] = boostWorld.irecv(rankX1Y2Z2, mpiTag, rParticleX1Y2Z2);
    }
    if (rankX2Y1Z1 >= 0) { // edge x2y1z1
      Rectangle containerX2Y1Z1(v2.getX() - cellSize, v1.getY(), v1.getZ(),
        v2.getX(), v1.getY() + cellSize, v1.getZ() + cellSize);
      findParticleInRectangle(containerX2Y1Z1, particleVec, particleX2Y1Z1);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX2Y1Z1.begin(); it!=particleX2Y1Z1.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX2Y1Z1[0] = boostWorld.isend(rankX2Y1Z1, mpiTag,  particleX2Y1Z1);
      reqX2Y1Z1[1] = boostWorld.irecv(rankX2Y1Z1, mpiTag, rParticleX2Y1Z1);
    }
    if (rankX2Y1Z2 >= 0) { // edge x2y1z2
      Rectangle containerX2Y1Z2(v2.getX() - cellSize, v1.getY(), v2.getZ() - cellSize,
        v2.getX(), v1.getY() + cellSize, v2.getZ());
      findParticleInRectangle(containerX2Y1Z2, particleVec, particleX2Y1Z2);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX2Y1Z2.begin(); it!=particleX2Y1Z2.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX2Y1Z2[0] = boostWorld.isend(rankX2Y1Z2, mpiTag,  particleX2Y1Z2);
      reqX2Y1Z2[1] = boostWorld.irecv(rankX2Y1Z2, mpiTag, rParticleX2Y1Z2);
    }
    if (rankX2Y2Z1 >= 0) { // edge x2y2z1
      Rectangle containerX2Y2Z1(v2.getX() - cellSize, v2.getY() - cellSize, v1.getZ(),
        v2.getX(), v2.getY(), v1.getZ() + cellSize);
      findParticleInRectangle(containerX2Y2Z1, particleVec, particleX2Y2Z1);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX2Y2Z1.begin(); it!=particleX2Y2Z1.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX2Y2Z1[0] = boostWorld.isend(rankX2Y2Z1, mpiTag,  particleX2Y2Z1);
      reqX2Y2Z1[1] = boostWorld.irecv(rankX2Y2Z1, mpiTag, rParticleX2Y2Z1);
    }
    if (rankX2Y2Z2 >= 0) { // edge x2y2z2
      Rectangle containerX2Y2Z2(v2.getX() - cellSize, v2.getY() - cellSize, v2.getZ() - cellSize,
        v2.getX(), v2.getY(), v2.getZ());
      findParticleInRectangle(containerX2Y2Z2, particleVec, particleX2Y2Z2);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX2Y2Z2.begin(); it!=particleX2Y2Z2.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX2Y2Z2[0] = boostWorld.isend(rankX2Y2Z2, mpiTag,  particleX2Y2Z2);
      reqX2Y2Z2[1] = boostWorld.irecv(rankX2Y2Z2, mpiTag, rParticleX2Y2Z2);
    }

    // 6 surfaces
    if (rankX1 >= 0) boost::mpi::wait_all(reqX1, reqX1 + 2);
    if (rankX2 >= 0) boost::mpi::wait_all(reqX2, reqX2 + 2);
    if (rankY1 >= 0) boost::mpi::wait_all(reqY1, reqY1 + 2);
    if (rankY2 >= 0) boost::mpi::wait_all(reqY2, reqY2 + 2);
    if (rankZ1 >= 0) boost::mpi::wait_all(reqZ1, reqZ1 + 2);
    if (rankZ2 >= 0) boost::mpi::wait_all(reqZ2, reqZ2 + 2);
    // 12 edges
    if (rankX1Y1 >= 0) boost::mpi::wait_all(reqX1Y1, reqX1Y1 + 2);
    if (rankX1Y2 >= 0) boost::mpi::wait_all(reqX1Y2, reqX1Y2 + 2);  
    if (rankX1Z1 >= 0) boost::mpi::wait_all(reqX1Z1, reqX1Z1 + 2);
    if (rankX1Z2 >= 0) boost::mpi::wait_all(reqX1Z2, reqX1Z2 + 2);
    if (rankX2Y1 >= 0) boost::mpi::wait_all(reqX2Y1, reqX2Y1 + 2);
    if (rankX2Y2 >= 0) boost::mpi::wait_all(reqX2Y2, reqX2Y2 + 2);  
    if (rankX2Z1 >= 0) boost::mpi::wait_all(reqX2Z1, reqX2Z1 + 2);
    if (rankX2Z2 >= 0) boost::mpi::wait_all(reqX2Z2, reqX2Z2 + 2); 
    if (rankY1Z1 >= 0) boost::mpi::wait_all(reqY1Z1, reqY1Z1 + 2);
    if (rankY1Z2 >= 0) boost::mpi::wait_all(reqY1Z2, reqY1Z2 + 2);
    if (rankY2Z1 >= 0) boost::mpi::wait_all(reqY2Z1, reqY2Z1 + 2);
    if (rankY2Z2 >= 0) boost::mpi::wait_all(reqY2Z2, reqY2Z2 + 2); 
    // 8 vertices
    if (rankX1Y1Z1 >= 0) boost::mpi::wait_all(reqX1Y1Z1, reqX1Y1Z1 + 2);
    if (rankX1Y1Z2 >= 0) boost::mpi::wait_all(reqX1Y1Z2, reqX1Y1Z2 + 2);
    if (rankX1Y2Z1 >= 0) boost::mpi::wait_all(reqX1Y2Z1, reqX1Y2Z1 + 2);
    if (rankX1Y2Z2 >= 0) boost::mpi::wait_all(reqX1Y2Z2, reqX1Y2Z2 + 2);
    if (rankX2Y1Z1 >= 0) boost::mpi::wait_all(reqX2Y1Z1, reqX2Y1Z1 + 2);
    if (rankX2Y1Z2 >= 0) boost::mpi::wait_all(reqX2Y1Z2, reqX2Y1Z2 + 2);
    if (rankX2Y2Z1 >= 0) boost::mpi::wait_all(reqX2Y2Z1, reqX2Y2Z1 + 2);
    if (rankX2Y2Z2 >= 0) boost::mpi::wait_all(reqX2Y2Z2, reqX2Y2Z2 + 2);  

    // merge: particles inside container (at front) + particles from neighoring blocks (at end)
    recvParticleVec.clear();
    // 6 surfaces
    if (rankX1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1.begin(), rParticleX1.end());
    if (rankX2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2.begin(), rParticleX2.end());
    if (rankY1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY1.begin(), rParticleY1.end());
    if (rankY2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY2.begin(), rParticleY2.end());
    if (rankZ1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleZ1.begin(), rParticleZ1.end());
    if (rankZ2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleZ2.begin(), rParticleZ2.end());
    // 12 edges
    if (rankX1Y1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1.begin(), rParticleX1Y1.end());
    if (rankX1Y2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2.begin(), rParticleX1Y2.end());
    if (rankX1Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Z1.begin(), rParticleX1Z1.end());
    if (rankX1Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Z2.begin(), rParticleX1Z2.end());
    if (rankX2Y1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1.begin(), rParticleX2Y1.end());
    if (rankX2Y2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2.begin(), rParticleX2Y2.end());
    if (rankX2Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Z1.begin(), rParticleX2Z1.end());
    if (rankX2Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Z2.begin(), rParticleX2Z2.end());
    if (rankY1Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY1Z1.begin(), rParticleY1Z1.end());
    if (rankY1Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY1Z2.begin(), rParticleY1Z2.end());
    if (rankY2Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY2Z1.begin(), rParticleY2Z1.end());
    if (rankY2Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY2Z2.begin(), rParticleY2Z2.end());
    // 8 vertices
    if (rankX1Y1Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1Z1.begin(), rParticleX1Y1Z1.end());
    if (rankX1Y1Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1Z2.begin(), rParticleX1Y1Z2.end());
    if (rankX1Y2Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2Z1.begin(), rParticleX1Y2Z1.end());
    if (rankX1Y2Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2Z2.begin(), rParticleX1Y2Z2.end());
    if (rankX2Y1Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1Z1.begin(), rParticleX2Y1Z1.end());
    if (rankX2Y1Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1Z2.begin(), rParticleX2Y1Z2.end());
    if (rankX2Y2Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2Z1.begin(), rParticleX2Y2Z1.end());
    if (rankX2Y2Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2Z2.begin(), rParticleX2Y2Z2.end());

    // after receive, set SPHParticle.demParticle
    for(std::vector<Particle*>::iterator it=recvParticleVec.begin(); it!=recvParticleVec.end(); it++){
  (*it)->setDemParticleInSPHParticle();
    }

    mergeParticleVec.clear();
    mergeParticleVec = particleVec; // duplicate pointers, pointing to the same memory
    mergeParticleVec.insert(mergeParticleVec.end(), recvParticleVec.begin(), recvParticleVec.end());

    // at this point, the SPHGhostParticles in part of particleVec are pointing to NULL, needed to assign back
    // this must be done after the receive is complete

    // 6 surfaces
    if (rankX1 >= 0)
      for(std::vector<Particle*>::iterator it=particleX1.begin(); it!=particleX1.end(); it++)
  (*it)->setDemParticleInSPHParticle();
    if (rankX2 >= 0) 
      for(std::vector<Particle*>::iterator it=particleX2.begin(); it!=particleX2.end(); it++)
  (*it)->setDemParticleInSPHParticle();
    if (rankY1 >= 0) 
      for(std::vector<Particle*>::iterator it=particleY1.begin(); it!=particleY1.end(); it++)
  (*it)->setDemParticleInSPHParticle();
    if (rankY2 >= 0) 
      for(std::vector<Particle*>::iterator it=particleY2.begin(); it!=particleY2.end(); it++)
  (*it)->setDemParticleInSPHParticle();
    if (rankZ1 >= 0) 
      for(std::vector<Particle*>::iterator it=particleZ1.begin(); it!=particleZ1.end(); it++)
  (*it)->setDemParticleInSPHParticle();
    if (rankZ2 >= 0) 
      for(std::vector<Particle*>::iterator it=particleZ2.begin(); it!=particleZ2.end(); it++)
  (*it)->setDemParticleInSPHParticle();
    // 12 edges
    if (rankX1Y1 >= 0) 
      for(std::vector<Particle*>::iterator it=particleX1Y1.begin(); it!=particleX1Y1.end(); it++)
  (*it)->setDemParticleInSPHParticle();
    if (rankX1Y2 >= 0) 
      for(std::vector<Particle*>::iterator it=particleX1Y2.begin(); it!=particleX1Y2.end(); it++)
  (*it)->setDemParticleInSPHParticle();
    if (rankX1Z1 >= 0) 
      for(std::vector<Particle*>::iterator it=particleX1Z1.begin(); it!=particleX1Z1.end(); it++)
  (*it)->setDemParticleInSPHParticle();
    if (rankX1Z2 >= 0) 
      for(std::vector<Particle*>::iterator it=particleX1Z2.begin(); it!=particleX1Z2.end(); it++)
  (*it)->setDemParticleInSPHParticle();
    if (rankX2Y1 >= 0) 
      for(std::vector<Particle*>::iterator it=particleX2Y1.begin(); it!=particleX2Y1.end(); it++)
  (*it)->setDemParticleInSPHParticle();
    if (rankX2Y2 >= 0) 
      for(std::vector<Particle*>::iterator it=particleX2Y2.begin(); it!=particleX2Y2.end(); it++)
  (*it)->setDemParticleInSPHParticle();
    if (rankX2Z1 >= 0) 
      for(std::vector<Particle*>::iterator it=particleX2Z1.begin(); it!=particleX2Z1.end(); it++)
  (*it)->setDemParticleInSPHParticle();
    if (rankX2Z2 >= 0) 
      for(std::vector<Particle*>::iterator it=particleX2Z2.begin(); it!=particleX2Z2.end(); it++)
  (*it)->setDemParticleInSPHParticle();
    if (rankY1Z1 >= 0)
      for(std::vector<Particle*>::iterator it=particleY1Z1.begin(); it!=particleY1Z1.end(); it++)
  (*it)->setDemParticleInSPHParticle();
    if (rankY1Z2 >= 0) 
      for(std::vector<Particle*>::iterator it=particleY1Z2.begin(); it!=particleY1Z2.end(); it++)
  (*it)->setDemParticleInSPHParticle();
    if (rankY2Z1 >= 0) 
      for(std::vector<Particle*>::iterator it=particleY2Z1.begin(); it!=particleY2Z1.end(); it++)
  (*it)->setDemParticleInSPHParticle();
    if (rankY2Z2 >= 0) 
      for(std::vector<Particle*>::iterator it=particleY2Z2.begin(); it!=particleY2Z2.end(); it++)
  (*it)->setDemParticleInSPHParticle();
    // 8 vertices
    if (rankX1Y1Z1 >= 0) 
      for(std::vector<Particle*>::iterator it=particleX1Y1Z1.begin(); it!=particleX1Y1Z1.end(); it++)
  (*it)->setDemParticleInSPHParticle();
    if (rankX1Y1Z2 >= 0) 
      for(std::vector<Particle*>::iterator it=particleX1Y1Z2.begin(); it!=particleX1Y1Z2.end(); it++)
  (*it)->setDemParticleInSPHParticle();
    if (rankX1Y2Z1 >= 0) 
      for(std::vector<Particle*>::iterator it=particleX1Y2Z1.begin(); it!=particleX1Y2Z1.end(); it++)
  (*it)->setDemParticleInSPHParticle();
    if (rankX1Y2Z2 >= 0) 
      for(std::vector<Particle*>::iterator it=particleX1Y2Z2.begin(); it!=particleX1Y2Z2.end(); it++)
  (*it)->setDemParticleInSPHParticle();
    if (rankX2Y1Z1 >= 0) 
      for(std::vector<Particle*>::iterator it=particleX2Y1Z1.begin(); it!=particleX2Y1Z1.end(); it++)
  (*it)->setDemParticleInSPHParticle();
    if (rankX2Y1Z2 >= 0) 
      for(std::vector<Particle*>::iterator it=particleX2Y1Z2.begin(); it!=particleX2Y1Z2.end(); it++)
  (*it)->setDemParticleInSPHParticle();
    if (rankX2Y2Z1 >= 0) 
      for(std::vector<Particle*>::iterator it=particleX2Y2Z1.begin(); it!=particleX2Y2Z1.end(); it++)
  (*it)->setDemParticleInSPHParticle();
    if (rankX2Y2Z2 >= 0) 
      for(std::vector<Particle*>::iterator it=particleX2Y2Z2.begin(); it!=particleX2Y2Z2.end(); it++)
  (*it)->setDemParticleInSPHParticle();
    

    /*
      std::vector<Particle*> testParticleVec;
      testParticleVec.insert(testParticleVec.end(), rParticleX1.begin(), rParticleX1.end());
      testParticleVec.insert(testParticleVec.end(), rParticleX2.begin(), rParticleX2.end());
      testParticleVec.insert(testParticleVec.end(), rParticleY1.begin(), rParticleY1.end());
      testParticleVec.insert(testParticleVec.end(), rParticleY2.begin(), rParticleY2.end());
      testParticleVec.insert(testParticleVec.end(), rParticleZ1.begin(), rParticleZ1.end());
      testParticleVec.insert(testParticleVec.end(), rParticleZ2.begin(), rParticleZ2.end());
      debugInf << "iter=" << std::setw(4) << iteration << " rank=" << std::setw(4) << mpiRank 
      << " ptclNum=" << std::setw(4) << particleVec.size() 
      << " surface="
      << std::setw(4) << particleX1.size()  << std::setw(4) << particleX2.size()
      << std::setw(4) << particleY1.size()  << std::setw(4) << particleY2.size()
      << std::setw(4) << particleZ1.size()  << std::setw(4) << particleZ2.size()  
      << " recv="
      << std::setw(4) << rParticleX1.size() << std::setw(4) << rParticleX2.size()
      << std::setw(4) << rParticleY1.size() << std::setw(4) << rParticleY2.size()
      << std::setw(4) << rParticleZ1.size() << std::setw(4) << rParticleZ2.size() 
      << " rNum="    
      << std::setw(4) << recvParticleVec.size() << ": ";   

      for (std::vector<Particle*>::const_iterator it = testParticleVec.begin(); it != testParticleVec.end();++it)
      debugInf << (*it)->getId() << ' ';
      debugInf << std::endl;
      testParticleVec.clear();
    */
  }

  // the communication of ghost sph particles are not implemented, the communication of ghost sph should be
  // corresponding to their dem particles. July 16, 2015
  void SmoothParticleHydrodynamics::commuSPHParticle() 
  {
    // determine container of each process
    Vec v1 = grid.getMinCorner();
    Vec v2 = grid.getMaxCorner();
    Vec vspan = v2 - v1;
    container = Rectangle(v1.getX() + vspan.getX() / mpiProcX * mpiCoords[0],
        v1.getY() + vspan.getY() / mpiProcY * mpiCoords[1],
        v1.getZ() + vspan.getZ() / mpiProcZ * mpiCoords[2],
        v1.getX() + vspan.getX() / mpiProcX * (mpiCoords[0] + 1),
        v1.getY() + vspan.getY() / mpiProcY * (mpiCoords[1] + 1),
        v1.getZ() + vspan.getZ() / mpiProcZ * (mpiCoords[2] + 1));

    // if found, communicate with neighboring blocks
    std::vector<sph::SPHParticle*> sphParticleX1, sphParticleX2;
    std::vector<sph::SPHParticle*> sphParticleY1, sphParticleY2;
    std::vector<sph::SPHParticle*> sphParticleZ1, sphParticleZ2;
    std::vector<sph::SPHParticle*> sphParticleX1Y1, sphParticleX1Y2, sphParticleX1Z1, sphParticleX1Z2; 
    std::vector<sph::SPHParticle*> sphParticleX2Y1, sphParticleX2Y2, sphParticleX2Z1, sphParticleX2Z2; 
    std::vector<sph::SPHParticle*> sphParticleY1Z1, sphParticleY1Z2, sphParticleY2Z1, sphParticleY2Z2; 
    std::vector<sph::SPHParticle*> sphParticleX1Y1Z1, sphParticleX1Y1Z2, sphParticleX1Y2Z1, sphParticleX1Y2Z2; 
    std::vector<sph::SPHParticle*> sphParticleX2Y1Z1, sphParticleX2Y1Z2, sphParticleX2Y2Z1, sphParticleX2Y2Z2; 
    boost::mpi::request reqX1[2], reqX2[2];
    boost::mpi::request reqY1[2], reqY2[2];
    boost::mpi::request reqZ1[2], reqZ2[2];
    boost::mpi::request reqX1Y1[2], reqX1Y2[2], reqX1Z1[2], reqX1Z2[2];
    boost::mpi::request reqX2Y1[2], reqX2Y2[2], reqX2Z1[2], reqX2Z2[2];
    boost::mpi::request reqY1Z1[2], reqY1Z2[2], reqY2Z1[2], reqY2Z2[2];
    boost::mpi::request reqX1Y1Z1[2], reqX1Y1Z2[2], reqX1Y2Z1[2], reqX1Y2Z2[2];
    boost::mpi::request reqX2Y1Z1[2], reqX2Y1Z2[2], reqX2Y2Z1[2], reqX2Y2Z2[2];
    v1 = container.getMinCorner(); // redefine v1, v2 in terms of process
    v2 = container.getMaxCorner();   
    //debugInf << "rank=" << mpiRank << ' ' << v1.getX() << ' ' << v1.getY() << ' ' << v1.getZ() << ' '  << v2.getX() << ' ' << v2.getY() << ' ' << v2.getZ() << std::endl;
    REAL cellSize = sphCellSize;
    // 6 surfaces
    if (rankX1 >= 0) { // surface x1
      Rectangle containerX1(v1.getX(), v1.getY(), v1.getZ(), 
          v1.getX() + cellSize, v2.getY(), v2.getZ());
      findSPHParticleInRectangle(containerX1, SPHParticleVec, sphParticleX1);
      reqX1[0] = boostWorld.isend(rankX1, mpiTag,  sphParticleX1);
      reqX1[1] = boostWorld.irecv(rankX1, mpiTag, rsphParticleX1);
    }
    if (rankX2 >= 0) { // surface x2
      Rectangle containerX2(v2.getX() - cellSize, v1.getY(), v1.getZ(),
          v2.getX(), v2.getY(), v2.getZ());
      findSPHParticleInRectangle(containerX2, SPHParticleVec, sphParticleX2);
      reqX2[0] = boostWorld.isend(rankX2, mpiTag,  sphParticleX2);
      reqX2[1] = boostWorld.irecv(rankX2, mpiTag, rsphParticleX2);
    }
    if (rankY1 >= 0) {  // surface y1
      Rectangle containerY1(v1.getX(), v1.getY(), v1.getZ(), 
          v2.getX(), v1.getY() + cellSize, v2.getZ());
      findSPHParticleInRectangle(containerY1, SPHParticleVec, sphParticleY1);
      reqY1[0] = boostWorld.isend(rankY1, mpiTag,  sphParticleY1);
      reqY1[1] = boostWorld.irecv(rankY1, mpiTag, rsphParticleY1);
    }
    if (rankY2 >= 0) {  // surface y2
      Rectangle containerY2(v1.getX(), v2.getY() - cellSize, v1.getZ(),
          v2.getX(), v2.getY(), v2.getZ());
      findSPHParticleInRectangle(containerY2, SPHParticleVec, sphParticleY2);
      reqY2[0] = boostWorld.isend(rankY2, mpiTag,  sphParticleY2);
      reqY2[1] = boostWorld.irecv(rankY2, mpiTag, rsphParticleY2);
    }
    if (rankZ1 >= 0) {  // surface z1
      Rectangle containerZ1(v1.getX(), v1.getY(), v1.getZ(),
          v2.getX(), v2.getY(), v1.getZ() + cellSize);
      findSPHParticleInRectangle(containerZ1, SPHParticleVec, sphParticleZ1);
      reqZ1[0] = boostWorld.isend(rankZ1, mpiTag,  sphParticleZ1);
      reqZ1[1] = boostWorld.irecv(rankZ1, mpiTag, rsphParticleZ1);
    }
    if (rankZ2 >= 0) {  // surface z2
      Rectangle containerZ2(v1.getX(), v1.getY(), v2.getZ() - cellSize,
          v2.getX(), v2.getY(), v2.getZ());
      findSPHParticleInRectangle(containerZ2, SPHParticleVec, sphParticleZ2);
      reqZ2[0] = boostWorld.isend(rankZ2, mpiTag,  sphParticleZ2);
      reqZ2[1] = boostWorld.irecv(rankZ2, mpiTag, rsphParticleZ2);
    }
    // 12 edges
    if (rankX1Y1 >= 0) { // edge x1y1
      Rectangle containerX1Y1(v1.getX(), v1.getY(), v1.getZ(),
            v1.getX() + cellSize, v1.getY() + cellSize, v2.getZ());
      findSPHParticleInRectangle(containerX1Y1, SPHParticleVec, sphParticleX1Y1);
      reqX1Y1[0] = boostWorld.isend(rankX1Y1, mpiTag,  sphParticleX1Y1);
      reqX1Y1[1] = boostWorld.irecv(rankX1Y1, mpiTag, rsphParticleX1Y1);
    }
    if (rankX1Y2 >= 0) { // edge x1y2
      Rectangle containerX1Y2(v1.getX(), v2.getY() - cellSize, v1.getZ(),
            v1.getX() + cellSize, v2.getY(), v2.getZ());
      findSPHParticleInRectangle(containerX1Y2, SPHParticleVec, sphParticleX1Y2);
      reqX1Y2[0] = boostWorld.isend(rankX1Y2, mpiTag,  sphParticleX1Y2);
      reqX1Y2[1] = boostWorld.irecv(rankX1Y2, mpiTag, rsphParticleX1Y2);
    }
    if (rankX1Z1 >= 0) { // edge x1z1
      Rectangle containerX1Z1(v1.getX(), v1.getY(), v1.getZ(),
            v1.getX() + cellSize, v2.getY(), v1.getZ() + cellSize);
      findSPHParticleInRectangle(containerX1Z1, SPHParticleVec, sphParticleX1Z1);
      reqX1Z1[0] = boostWorld.isend(rankX1Z1, mpiTag,  sphParticleX1Z1);
      reqX1Z1[1] = boostWorld.irecv(rankX1Z1, mpiTag, rsphParticleX1Z1);
    }
    if (rankX1Z2 >= 0) { // edge x1z2
      Rectangle containerX1Z2(v1.getX(), v1.getY(), v2.getZ() - cellSize,
            v1.getX() + cellSize, v2.getY(), v2.getZ());
      findSPHParticleInRectangle(containerX1Z2, SPHParticleVec, sphParticleX1Z2);
      reqX1Z2[0] = boostWorld.isend(rankX1Z2, mpiTag,  sphParticleX1Z2);
      reqX1Z2[1] = boostWorld.irecv(rankX1Z2, mpiTag, rsphParticleX1Z2);
    }
    if (rankX2Y1 >= 0) { // edge x2y1
      Rectangle containerX2Y1(v2.getX() - cellSize, v1.getY(), v1.getZ(),
            v2.getX(), v1.getY() + cellSize, v2.getZ());
      findSPHParticleInRectangle(containerX2Y1, SPHParticleVec, sphParticleX2Y1);
      reqX2Y1[0] = boostWorld.isend(rankX2Y1, mpiTag,  sphParticleX2Y1);
      reqX2Y1[1] = boostWorld.irecv(rankX2Y1, mpiTag, rsphParticleX2Y1);
    }
    if (rankX2Y2 >= 0) { // edge x2y2
      Rectangle containerX2Y2(v2.getX() - cellSize, v2.getY() - cellSize, v1.getZ(),
            v2.getX(), v2.getY(), v2.getZ());
      findSPHParticleInRectangle(containerX2Y2, SPHParticleVec, sphParticleX2Y2);
      reqX2Y2[0] = boostWorld.isend(rankX2Y2, mpiTag,  sphParticleX2Y2);
      reqX2Y2[1] = boostWorld.irecv(rankX2Y2, mpiTag, rsphParticleX2Y2);
    }
    if (rankX2Z1 >= 0) { // edge x2z1
      Rectangle containerX2Z1(v2.getX() - cellSize, v1.getY(), v1.getZ(),
            v2.getX(), v2.getY(), v1.getZ() + cellSize);
      findSPHParticleInRectangle(containerX2Z1, SPHParticleVec, sphParticleX2Z1);
      reqX2Z1[0] = boostWorld.isend(rankX2Z1, mpiTag,  sphParticleX2Z1);
      reqX2Z1[1] = boostWorld.irecv(rankX2Z1, mpiTag, rsphParticleX2Z1);
    }
    if (rankX2Z2 >= 0) { // edge x2z2
      Rectangle containerX2Z2(v2.getX() - cellSize, v1.getY(), v2.getZ() - cellSize,
            v2.getX(), v2.getY(), v2.getZ());
      findSPHParticleInRectangle(containerX2Z2, SPHParticleVec, sphParticleX2Z2);
      reqX2Z2[0] = boostWorld.isend(rankX2Z2, mpiTag,  sphParticleX2Z2);
      reqX2Z2[1] = boostWorld.irecv(rankX2Z2, mpiTag, rsphParticleX2Z2);
    }
    if (rankY1Z1 >= 0) { // edge y1z1
      Rectangle containerY1Z1(v1.getX(), v1.getY(), v1.getZ(),
            v2.getX(), v1.getY() + cellSize, v1.getZ() + cellSize);
      findSPHParticleInRectangle(containerY1Z1, SPHParticleVec, sphParticleY1Z1);
      reqY1Z1[0] = boostWorld.isend(rankY1Z1, mpiTag,  sphParticleY1Z1);
      reqY1Z1[1] = boostWorld.irecv(rankY1Z1, mpiTag, rsphParticleY1Z1);
    }
    if (rankY1Z2 >= 0) { // edge y1z2
      Rectangle containerY1Z2(v1.getX(), v1.getY(), v2.getZ() - cellSize,
            v2.getX(), v1.getY() + cellSize, v2.getZ());
      findSPHParticleInRectangle(containerY1Z2, SPHParticleVec, sphParticleY1Z2);
      reqY1Z2[0] = boostWorld.isend(rankY1Z2, mpiTag,  sphParticleY1Z2);
      reqY1Z2[1] = boostWorld.irecv(rankY1Z2, mpiTag, rsphParticleY1Z2);
    }
    if (rankY2Z1 >= 0) { // edge y2z1
      Rectangle containerY2Z1(v1.getX(), v2.getY() - cellSize, v1.getZ(),
            v2.getX(), v2.getY(), v1.getZ() + cellSize);
      findSPHParticleInRectangle(containerY2Z1, SPHParticleVec, sphParticleY2Z1);
      reqY2Z1[0] = boostWorld.isend(rankY2Z1, mpiTag,  sphParticleY2Z1);
      reqY2Z1[1] = boostWorld.irecv(rankY2Z1, mpiTag, rsphParticleY2Z1);
    }
    if (rankY2Z2 >= 0) { // edge y2z2
      Rectangle containerY2Z2(v1.getX(), v2.getY() - cellSize, v2.getZ() - cellSize,
            v2.getX(), v2.getY(), v2.getZ());
      findSPHParticleInRectangle(containerY2Z2, SPHParticleVec, sphParticleY2Z2);
      reqY2Z2[0] = boostWorld.isend(rankY2Z2, mpiTag,  sphParticleY2Z2);
      reqY2Z2[1] = boostWorld.irecv(rankY2Z2, mpiTag, rsphParticleY2Z2);
    }
    // 8 vertices
    if (rankX1Y1Z1 >= 0) { // edge x1y1z1
      Rectangle containerX1Y1Z1(v1.getX(), v1.getY(), v1.getZ(),
        v1.getX() + cellSize, v1.getY() + cellSize, v1.getZ() + cellSize);
      findSPHParticleInRectangle(containerX1Y1Z1, SPHParticleVec, sphParticleX1Y1Z1);
      reqX1Y1Z1[0] = boostWorld.isend(rankX1Y1Z1, mpiTag,  sphParticleX1Y1Z1);
      reqX1Y1Z1[1] = boostWorld.irecv(rankX1Y1Z1, mpiTag, rsphParticleX1Y1Z1);
    }
    if (rankX1Y1Z2 >= 0) { // edge x1y1z2
      Rectangle containerX1Y1Z2(v1.getX(), v1.getY(), v2.getZ() - cellSize,
        v1.getX() + cellSize, v1.getY() + cellSize, v2.getZ());
      findSPHParticleInRectangle(containerX1Y1Z2, SPHParticleVec, sphParticleX1Y1Z2);
      reqX1Y1Z2[0] = boostWorld.isend(rankX1Y1Z2, mpiTag,  sphParticleX1Y1Z2);
      reqX1Y1Z2[1] = boostWorld.irecv(rankX1Y1Z2, mpiTag, rsphParticleX1Y1Z2);
    }
    if (rankX1Y2Z1 >= 0) { // edge x1y2z1
      Rectangle containerX1Y2Z1(v1.getX(), v2.getY() - cellSize, v1.getZ(),
        v1.getX() + cellSize, v2.getY(), v1.getZ() + cellSize);
      findSPHParticleInRectangle(containerX1Y2Z1, SPHParticleVec, sphParticleX1Y2Z1);
      reqX1Y2Z1[0] = boostWorld.isend(rankX1Y2Z1, mpiTag,  sphParticleX1Y2Z1);
      reqX1Y2Z1[1] = boostWorld.irecv(rankX1Y2Z1, mpiTag, rsphParticleX1Y2Z1);
    }
    if (rankX1Y2Z2 >= 0) { // edge x1y2z2
      Rectangle containerX1Y2Z2(v1.getX(), v2.getY() - cellSize, v2.getZ() - cellSize,
        v1.getX() + cellSize, v2.getY() + cellSize, v2.getZ());
      findSPHParticleInRectangle(containerX1Y2Z2, SPHParticleVec, sphParticleX1Y2Z2);
      reqX1Y2Z2[0] = boostWorld.isend(rankX1Y2Z2, mpiTag,  sphParticleX1Y2Z2);
      reqX1Y2Z2[1] = boostWorld.irecv(rankX1Y2Z2, mpiTag, rsphParticleX1Y2Z2);
    }
    if (rankX2Y1Z1 >= 0) { // edge x2y1z1
      Rectangle containerX2Y1Z1(v2.getX() - cellSize, v1.getY(), v1.getZ(),
        v2.getX(), v1.getY() + cellSize, v1.getZ() + cellSize);
      findSPHParticleInRectangle(containerX2Y1Z1, SPHParticleVec, sphParticleX2Y1Z1);
      reqX2Y1Z1[0] = boostWorld.isend(rankX2Y1Z1, mpiTag,  sphParticleX2Y1Z1);
      reqX2Y1Z1[1] = boostWorld.irecv(rankX2Y1Z1, mpiTag, rsphParticleX2Y1Z1);
    }
    if (rankX2Y1Z2 >= 0) { // edge x2y1z2
      Rectangle containerX2Y1Z2(v2.getX() - cellSize, v1.getY(), v2.getZ() - cellSize,
        v2.getX(), v1.getY() + cellSize, v2.getZ());
      findSPHParticleInRectangle(containerX2Y1Z2, SPHParticleVec, sphParticleX2Y1Z2);
      reqX2Y1Z2[0] = boostWorld.isend(rankX2Y1Z2, mpiTag,  sphParticleX2Y1Z2);
      reqX2Y1Z2[1] = boostWorld.irecv(rankX2Y1Z2, mpiTag, rsphParticleX2Y1Z2);
    }
    if (rankX2Y2Z1 >= 0) { // edge x2y2z1
      Rectangle containerX2Y2Z1(v2.getX() - cellSize, v2.getY() - cellSize, v1.getZ(),
        v2.getX(), v2.getY(), v1.getZ() + cellSize);
      findSPHParticleInRectangle(containerX2Y2Z1, SPHParticleVec, sphParticleX2Y2Z1);
      reqX2Y2Z1[0] = boostWorld.isend(rankX2Y2Z1, mpiTag,  sphParticleX2Y2Z1);
      reqX2Y2Z1[1] = boostWorld.irecv(rankX2Y2Z1, mpiTag, rsphParticleX2Y2Z1);
    }
    if (rankX2Y2Z2 >= 0) { // edge x2y2z2
      Rectangle containerX2Y2Z2(v2.getX() - cellSize, v2.getY() - cellSize, v2.getZ() - cellSize,
        v2.getX(), v2.getY(), v2.getZ());
      findSPHParticleInRectangle(containerX2Y2Z2, SPHParticleVec, sphParticleX2Y2Z2);
      reqX2Y2Z2[0] = boostWorld.isend(rankX2Y2Z2, mpiTag,  sphParticleX2Y2Z2);
      reqX2Y2Z2[1] = boostWorld.irecv(rankX2Y2Z2, mpiTag, rsphParticleX2Y2Z2);
    }

    // 6 surfaces
    if (rankX1 >= 0) boost::mpi::wait_all(reqX1, reqX1 + 2);
    if (rankX2 >= 0) boost::mpi::wait_all(reqX2, reqX2 + 2);
    if (rankY1 >= 0) boost::mpi::wait_all(reqY1, reqY1 + 2);
    if (rankY2 >= 0) boost::mpi::wait_all(reqY2, reqY2 + 2);
    if (rankZ1 >= 0) boost::mpi::wait_all(reqZ1, reqZ1 + 2);
    if (rankZ2 >= 0) boost::mpi::wait_all(reqZ2, reqZ2 + 2);
    // 12 edges
    if (rankX1Y1 >= 0) boost::mpi::wait_all(reqX1Y1, reqX1Y1 + 2);
    if (rankX1Y2 >= 0) boost::mpi::wait_all(reqX1Y2, reqX1Y2 + 2);  
    if (rankX1Z1 >= 0) boost::mpi::wait_all(reqX1Z1, reqX1Z1 + 2);
    if (rankX1Z2 >= 0) boost::mpi::wait_all(reqX1Z2, reqX1Z2 + 2);
    if (rankX2Y1 >= 0) boost::mpi::wait_all(reqX2Y1, reqX2Y1 + 2);
    if (rankX2Y2 >= 0) boost::mpi::wait_all(reqX2Y2, reqX2Y2 + 2);  
    if (rankX2Z1 >= 0) boost::mpi::wait_all(reqX2Z1, reqX2Z1 + 2);
    if (rankX2Z2 >= 0) boost::mpi::wait_all(reqX2Z2, reqX2Z2 + 2); 
    if (rankY1Z1 >= 0) boost::mpi::wait_all(reqY1Z1, reqY1Z1 + 2);
    if (rankY1Z2 >= 0) boost::mpi::wait_all(reqY1Z2, reqY1Z2 + 2);
    if (rankY2Z1 >= 0) boost::mpi::wait_all(reqY2Z1, reqY2Z1 + 2);
    if (rankY2Z2 >= 0) boost::mpi::wait_all(reqY2Z2, reqY2Z2 + 2); 
    // 8 vertices
    if (rankX1Y1Z1 >= 0) boost::mpi::wait_all(reqX1Y1Z1, reqX1Y1Z1 + 2);
    if (rankX1Y1Z2 >= 0) boost::mpi::wait_all(reqX1Y1Z2, reqX1Y1Z2 + 2);
    if (rankX1Y2Z1 >= 0) boost::mpi::wait_all(reqX1Y2Z1, reqX1Y2Z1 + 2);
    if (rankX1Y2Z2 >= 0) boost::mpi::wait_all(reqX1Y2Z2, reqX1Y2Z2 + 2);
    if (rankX2Y1Z1 >= 0) boost::mpi::wait_all(reqX2Y1Z1, reqX2Y1Z1 + 2);
    if (rankX2Y1Z2 >= 0) boost::mpi::wait_all(reqX2Y1Z2, reqX2Y1Z2 + 2);
    if (rankX2Y2Z1 >= 0) boost::mpi::wait_all(reqX2Y2Z1, reqX2Y2Z1 + 2);
    if (rankX2Y2Z2 >= 0) boost::mpi::wait_all(reqX2Y2Z2, reqX2Y2Z2 + 2);  

    // merge: sphParticles inside container (at front) + sphParticles from neighoring blocks (at end)
    recvSPHParticleVec.clear();
    // 6 surfaces
    if (rankX1 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX1.begin(), rsphParticleX1.end());
    if (rankX2 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX2.begin(), rsphParticleX2.end());
    if (rankY1 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleY1.begin(), rsphParticleY1.end());
    if (rankY2 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleY2.begin(), rsphParticleY2.end());
    if (rankZ1 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleZ1.begin(), rsphParticleZ1.end());
    if (rankZ2 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleZ2.begin(), rsphParticleZ2.end());
    // 12 edges
    if (rankX1Y1 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX1Y1.begin(), rsphParticleX1Y1.end());
    if (rankX1Y2 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX1Y2.begin(), rsphParticleX1Y2.end());
    if (rankX1Z1 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX1Z1.begin(), rsphParticleX1Z1.end());
    if (rankX1Z2 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX1Z2.begin(), rsphParticleX1Z2.end());
    if (rankX2Y1 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX2Y1.begin(), rsphParticleX2Y1.end());
    if (rankX2Y2 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX2Y2.begin(), rsphParticleX2Y2.end());
    if (rankX2Z1 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX2Z1.begin(), rsphParticleX2Z1.end());
    if (rankX2Z2 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX2Z2.begin(), rsphParticleX2Z2.end());
    if (rankY1Z1 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleY1Z1.begin(), rsphParticleY1Z1.end());
    if (rankY1Z2 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleY1Z2.begin(), rsphParticleY1Z2.end());
    if (rankY2Z1 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleY2Z1.begin(), rsphParticleY2Z1.end());
    if (rankY2Z2 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleY2Z2.begin(), rsphParticleY2Z2.end());
    // 8 vertices
    if (rankX1Y1Z1 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX1Y1Z1.begin(), rsphParticleX1Y1Z1.end());
    if (rankX1Y1Z2 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX1Y1Z2.begin(), rsphParticleX1Y1Z2.end());
    if (rankX1Y2Z1 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX1Y2Z1.begin(), rsphParticleX1Y2Z1.end());
    if (rankX1Y2Z2 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX1Y2Z2.begin(), rsphParticleX1Y2Z2.end());
    if (rankX2Y1Z1 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX2Y1Z1.begin(), rsphParticleX2Y1Z1.end());
    if (rankX2Y1Z2 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX2Y1Z2.begin(), rsphParticleX2Y1Z2.end());
    if (rankX2Y2Z1 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX2Y2Z1.begin(), rsphParticleX2Y2Z1.end());
    if (rankX2Y2Z2 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX2Y2Z2.begin(), rsphParticleX2Y2Z2.end());

    mergeSPHParticleVec.clear();
    mergeSPHParticleVec = SPHParticleVec; // duplicate pointers, pointing to the same memory
    mergeSPHParticleVec.insert(mergeSPHParticleVec.end(), recvSPHParticleVec.begin(), recvSPHParticleVec.end());

    /*
      std::vector<Particle*> testParticleVec;
      testParticleVec.insert(testParticleVec.end(), rParticleX1.begin(), rParticleX1.end());
      testParticleVec.insert(testParticleVec.end(), rParticleX2.begin(), rParticleX2.end());
      testParticleVec.insert(testParticleVec.end(), rParticleY1.begin(), rParticleY1.end());
      testParticleVec.insert(testParticleVec.end(), rParticleY2.begin(), rParticleY2.end());
      testParticleVec.insert(testParticleVec.end(), rParticleZ1.begin(), rParticleZ1.end());
      testParticleVec.insert(testParticleVec.end(), rParticleZ2.begin(), rParticleZ2.end());
      debugInf << "iter=" << std::setw(4) << iteration << " rank=" << std::setw(4) << mpiRank 
      << " ptclNum=" << std::setw(4) << particleVec.size() 
      << " surface="
      << std::setw(4) << particleX1.size()  << std::setw(4) << particleX2.size()
      << std::setw(4) << particleY1.size()  << std::setw(4) << particleY2.size()
      << std::setw(4) << particleZ1.size()  << std::setw(4) << particleZ2.size()  
      << " recv="
      << std::setw(4) << rParticleX1.size() << std::setw(4) << rParticleX2.size()
      << std::setw(4) << rParticleY1.size() << std::setw(4) << rParticleY2.size()
      << std::setw(4) << rParticleZ1.size() << std::setw(4) << rParticleZ2.size() 
      << " rNum="    
      << std::setw(4) << recvParticleVec.size() << ": ";   

      for (std::vector<Particle*>::const_iterator it = testParticleVec.begin(); it != testParticleVec.end();++it)
      debugInf << (*it)->getId() << ' ';
      debugInf << std::endl;
      testParticleVec.clear();
    */
  }


  void SmoothParticleHydrodynamics::releaseRecvParticle() {
    // release memory of received particles
    for (std::vector<Particle*>::iterator it = recvParticleVec.begin(); it != recvParticleVec.end(); ++it){
      for(std::vector<sph::SPHParticle*>::iterator st=(*it)->SPHGhostParticleVec.begin(); st!=(*it)->SPHGhostParticleVec.end(); st++){
  delete (*st);  // release memory of these sph ghost particles
      }
      (*it)->SPHGhostParticleVec.clear();
      delete (*it);
    }
    recvParticleVec.clear();
    // 6 surfaces
    rParticleX1.clear();
    rParticleX2.clear();
    rParticleY1.clear();
    rParticleY2.clear();
    rParticleZ1.clear();
    rParticleZ2.clear();
    // 12 edges
    rParticleX1Y1.clear();
    rParticleX1Y2.clear();
    rParticleX1Z1.clear();
    rParticleX1Z2.clear();
    rParticleX2Y1.clear();
    rParticleX2Y2.clear();
    rParticleX2Z1.clear();
    rParticleX2Z2.clear();
    rParticleY1Z1.clear();
    rParticleY1Z2.clear();
    rParticleY2Z1.clear();
    rParticleY2Z2.clear();
    // 8 vertices
    rParticleX1Y1Z1.clear();
    rParticleX1Y1Z2.clear();
    rParticleX1Y2Z1.clear();
    rParticleX1Y2Z2.clear();
    rParticleX2Y1Z1.clear();
    rParticleX2Y1Z2.clear();
    rParticleX2Y2Z1.clear();
    rParticleX2Y2Z2.clear();
  }


  void SmoothParticleHydrodynamics::releaseRecvSPHParticle() {
    // release memory of received particles
    for (std::vector<sph::SPHParticle*>::iterator it = recvSPHParticleVec.begin(); it != recvSPHParticleVec.end(); ++it)
      delete (*it);
    recvSPHParticleVec.clear();
    // 6 surfaces
    rsphParticleX1.clear();
    rsphParticleX2.clear();
    rsphParticleY1.clear();
    rsphParticleY2.clear();
    rsphParticleZ1.clear();
    rsphParticleZ2.clear();
    // 12 edges
    rsphParticleX1Y1.clear();
    rsphParticleX1Y2.clear();
    rsphParticleX1Z1.clear();
    rsphParticleX1Z2.clear();
    rsphParticleX2Y1.clear();
    rsphParticleX2Y2.clear();
    rsphParticleX2Z1.clear();
    rsphParticleX2Z2.clear();
    rsphParticleY1Z1.clear();
    rsphParticleY1Z2.clear();
    rsphParticleY2Z1.clear();
    rsphParticleY2Z2.clear();
    // 8 vertices
    rsphParticleX1Y1Z1.clear();
    rsphParticleX1Y1Z2.clear();
    rsphParticleX1Y2Z1.clear();
    rsphParticleX1Y2Z2.clear();
    rsphParticleX2Y1Z1.clear();
    rsphParticleX2Y1Z2.clear();
    rsphParticleX2Y2Z1.clear();
    rsphParticleX2Y2Z2.clear();
  }

  void SmoothParticleHydrodynamics::migrateParticle() 
  {
    Vec vspan = grid.getMaxCorner() - grid.getMinCorner();
    REAL segX = vspan.getX() / mpiProcX;
    REAL segY = vspan.getY() / mpiProcY;
    REAL segZ = vspan.getZ() / mpiProcZ;
    Vec v1 = container.getMinCorner(); // v1, v2 in terms of process
    Vec v2 = container.getMaxCorner();  

    // if a neighbor exists, transfer particles crossing the boundary in between.
    std::vector<Particle*> particleX1, particleX2;
    std::vector<Particle*> particleY1, particleY2;
    std::vector<Particle*> particleZ1, particleZ2;
    std::vector<Particle*> particleX1Y1, particleX1Y2, particleX1Z1, particleX1Z2; 
    std::vector<Particle*> particleX2Y1, particleX2Y2, particleX2Z1, particleX2Z2; 
    std::vector<Particle*> particleY1Z1, particleY1Z2, particleY2Z1, particleY2Z2; 
    std::vector<Particle*> particleX1Y1Z1, particleX1Y1Z2, particleX1Y2Z1, particleX1Y2Z2; 
    std::vector<Particle*> particleX2Y1Z1, particleX2Y1Z2, particleX2Y2Z1, particleX2Y2Z2; 
    boost::mpi::request reqX1[2], reqX2[2];
    boost::mpi::request reqY1[2], reqY2[2];
    boost::mpi::request reqZ1[2], reqZ2[2];
    boost::mpi::request reqX1Y1[2], reqX1Y2[2], reqX1Z1[2], reqX1Z2[2];
    boost::mpi::request reqX2Y1[2], reqX2Y2[2], reqX2Z1[2], reqX2Z2[2];
    boost::mpi::request reqY1Z1[2], reqY1Z2[2], reqY2Z1[2], reqY2Z2[2];
    boost::mpi::request reqX1Y1Z1[2], reqX1Y1Z2[2], reqX1Y2Z1[2], reqX1Y2Z2[2];
    boost::mpi::request reqX2Y1Z1[2], reqX2Y1Z2[2], reqX2Y2Z1[2], reqX2Y2Z2[2];

    // 6 surfaces
    if (rankX1 >= 0) { // surface x1
      Rectangle containerX1(v1.getX() - segX, v1.getY(), v1.getZ(), 
          v1.getX(), v2.getY(), v2.getZ());
      findParticleInRectangle(containerX1, particleVec, particleX1);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX1.begin(); it!=particleX1.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();  // now the SPHGhostParticles in part of particleVec are pointing to NULL, don't needed to be back in migrate
      reqX1[0] = boostWorld.isend(rankX1, mpiTag,  particleX1);
      reqX1[1] = boostWorld.irecv(rankX1, mpiTag, rParticleX1);
    }
    if (rankX2 >= 0) { // surface x2
      Rectangle containerX2(v2.getX(), v1.getY(), v1.getZ(),
          v2.getX() + segX, v2.getY(), v2.getZ());
      findParticleInRectangle(containerX2, particleVec, particleX2);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX2.begin(); it!=particleX2.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX2[0] = boostWorld.isend(rankX2, mpiTag,  particleX2);
      reqX2[1] = boostWorld.irecv(rankX2, mpiTag, rParticleX2);
    }
    if (rankY1 >= 0) {  // surface y1
      Rectangle containerY1(v1.getX(), v1.getY() - segY, v1.getZ(), 
          v2.getX(), v1.getY(), v2.getZ());
      findParticleInRectangle(containerY1, particleVec, particleY1);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleY1.begin(); it!=particleY1.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqY1[0] = boostWorld.isend(rankY1, mpiTag,  particleY1);
      reqY1[1] = boostWorld.irecv(rankY1, mpiTag, rParticleY1);
    }
    if (rankY2 >= 0) {  // surface y2
      Rectangle containerY2(v1.getX(), v2.getY(), v1.getZ(),
          v2.getX(), v2.getY() + segY, v2.getZ());
      findParticleInRectangle(containerY2, particleVec, particleY2);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleY2.begin(); it!=particleY2.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqY2[0] = boostWorld.isend(rankY2, mpiTag,  particleY2);
      reqY2[1] = boostWorld.irecv(rankY2, mpiTag, rParticleY2);
    }
    if (rankZ1 >= 0) {  // surface z1
      Rectangle containerZ1(v1.getX(), v1.getY(), v1.getZ() - segZ,
          v2.getX(), v2.getY(), v1.getZ());
      findParticleInRectangle(containerZ1, particleVec, particleZ1);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleZ1.begin(); it!=particleZ1.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqZ1[0] = boostWorld.isend(rankZ1, mpiTag,  particleZ1);
      reqZ1[1] = boostWorld.irecv(rankZ1, mpiTag, rParticleZ1);
    }
    if (rankZ2 >= 0) {  // surface z2
      Rectangle containerZ2(v1.getX(), v1.getY(), v2.getZ(),
          v2.getX(), v2.getY(), v2.getZ() + segZ);
      findParticleInRectangle(containerZ2, particleVec, particleZ2);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleZ2.begin(); it!=particleZ2.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqZ2[0] = boostWorld.isend(rankZ2, mpiTag,  particleZ2);
      reqZ2[1] = boostWorld.irecv(rankZ2, mpiTag, rParticleZ2);
    }
    // 12 edges
    if (rankX1Y1 >= 0) { // edge x1y1
      Rectangle containerX1Y1(v1.getX() - segX, v1.getY() - segY, v1.getZ(),
            v1.getX(), v1.getY(), v2.getZ());
      findParticleInRectangle(containerX1Y1, particleVec, particleX1Y1);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX1Y1.begin(); it!=particleX1Y1.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX1Y1[0] = boostWorld.isend(rankX1Y1, mpiTag,  particleX1Y1);
      reqX1Y1[1] = boostWorld.irecv(rankX1Y1, mpiTag, rParticleX1Y1);
    }
    if (rankX1Y2 >= 0) { // edge x1y2
      Rectangle containerX1Y2(v1.getX() - segX, v2.getY(), v1.getZ(),
            v1.getX(), v2.getY() + segY, v2.getZ());
      findParticleInRectangle(containerX1Y2, particleVec, particleX1Y2);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX1Y2.begin(); it!=particleX1Y2.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX1Y2[0] = boostWorld.isend(rankX1Y2, mpiTag,  particleX1Y2);
      reqX1Y2[1] = boostWorld.irecv(rankX1Y2, mpiTag, rParticleX1Y2);
    }
    if (rankX1Z1 >= 0) { // edge x1z1
      Rectangle containerX1Z1(v1.getX() - segX, v1.getY(), v1.getZ() -segZ,
            v1.getX(), v2.getY(), v1.getZ());
      findParticleInRectangle(containerX1Z1, particleVec, particleX1Z1);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX1Z1.begin(); it!=particleX1Z1.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX1Z1[0] = boostWorld.isend(rankX1Z1, mpiTag,  particleX1Z1);
      reqX1Z1[1] = boostWorld.irecv(rankX1Z1, mpiTag, rParticleX1Z1);
    }
    if (rankX1Z2 >= 0) { // edge x1z2
      Rectangle containerX1Z2(v1.getX() - segX, v1.getY(), v2.getZ(),
            v1.getX(), v2.getY(), v2.getZ() + segZ);
      findParticleInRectangle(containerX1Z2, particleVec, particleX1Z2);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX1Z2.begin(); it!=particleX1Z2.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX1Z2[0] = boostWorld.isend(rankX1Z2, mpiTag,  particleX1Z2);
      reqX1Z2[1] = boostWorld.irecv(rankX1Z2, mpiTag, rParticleX1Z2);
    }
    if (rankX2Y1 >= 0) { // edge x2y1
      Rectangle containerX2Y1(v2.getX(), v1.getY() - segY, v1.getZ(),
            v2.getX() + segX, v1.getY(), v2.getZ());
      findParticleInRectangle(containerX2Y1, particleVec, particleX2Y1);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX2Y1.begin(); it!=particleX2Y1.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX2Y1[0] = boostWorld.isend(rankX2Y1, mpiTag,  particleX2Y1);
      reqX2Y1[1] = boostWorld.irecv(rankX2Y1, mpiTag, rParticleX2Y1);
    }
    if (rankX2Y2 >= 0) { // edge x2y2
      Rectangle containerX2Y2(v2.getX(), v2.getY(), v1.getZ(),
            v2.getX() + segX, v2.getY() + segY, v2.getZ());
      findParticleInRectangle(containerX2Y2, particleVec, particleX2Y2);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX2Y2.begin(); it!=particleX2Y2.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX2Y2[0] = boostWorld.isend(rankX2Y2, mpiTag,  particleX2Y2);
      reqX2Y2[1] = boostWorld.irecv(rankX2Y2, mpiTag, rParticleX2Y2);
    }
    if (rankX2Z1 >= 0) { // edge x2z1
      Rectangle containerX2Z1(v2.getX(), v1.getY(), v1.getZ() - segZ,
            v2.getX() + segX, v2.getY(), v1.getZ());
      findParticleInRectangle(containerX2Z1, particleVec, particleX2Z1);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX2Z1.begin(); it!=particleX2Z1.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX2Z1[0] = boostWorld.isend(rankX2Z1, mpiTag,  particleX2Z1);
      reqX2Z1[1] = boostWorld.irecv(rankX2Z1, mpiTag, rParticleX2Z1);
    }
    if (rankX2Z2 >= 0) { // edge x2z2
      Rectangle containerX2Z2(v2.getX(), v1.getY(), v2.getZ(),
            v2.getX() + segX, v2.getY(), v2.getZ() + segZ);
      findParticleInRectangle(containerX2Z2, particleVec, particleX2Z2);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX2Z2.begin(); it!=particleX2Z2.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX2Z2[0] = boostWorld.isend(rankX2Z2, mpiTag,  particleX2Z2);
      reqX2Z2[1] = boostWorld.irecv(rankX2Z2, mpiTag, rParticleX2Z2);
    }
    if (rankY1Z1 >= 0) { // edge y1z1
      Rectangle containerY1Z1(v1.getX(), v1.getY() - segY, v1.getZ() - segZ,
            v2.getX(), v1.getY(), v1.getZ());
      findParticleInRectangle(containerY1Z1, particleVec, particleY1Z1);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleY1Z1.begin(); it!=particleY1Z1.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqY1Z1[0] = boostWorld.isend(rankY1Z1, mpiTag,  particleY1Z1);
      reqY1Z1[1] = boostWorld.irecv(rankY1Z1, mpiTag, rParticleY1Z1);
    }
    if (rankY1Z2 >= 0) { // edge y1z2
      Rectangle containerY1Z2(v1.getX(), v1.getY() - segY, v2.getZ(),
            v2.getX(), v1.getY(), v2.getZ() + segZ);
      findParticleInRectangle(containerY1Z2, particleVec, particleY1Z2);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleY1Z2.begin(); it!=particleY1Z2.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqY1Z2[0] = boostWorld.isend(rankY1Z2, mpiTag,  particleY1Z2);
      reqY1Z2[1] = boostWorld.irecv(rankY1Z2, mpiTag, rParticleY1Z2);
    }
    if (rankY2Z1 >= 0) { // edge y2z1
      Rectangle containerY2Z1(v1.getX(), v2.getY(), v1.getZ() - segZ,
            v2.getX(), v2.getY() + segY, v1.getZ());
      findParticleInRectangle(containerY2Z1, particleVec, particleY2Z1);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleY2Z1.begin(); it!=particleY2Z1.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqY2Z1[0] = boostWorld.isend(rankY2Z1, mpiTag,  particleY2Z1);
      reqY2Z1[1] = boostWorld.irecv(rankY2Z1, mpiTag, rParticleY2Z1);
    }
    if (rankY2Z2 >= 0) { // edge y2z2
      Rectangle containerY2Z2(v1.getX(), v2.getY(), v2.getZ(),
            v2.getX(), v2.getY() + segY, v2.getZ() + segZ);
      findParticleInRectangle(containerY2Z2, particleVec, particleY2Z2);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleY2Z2.begin(); it!=particleY2Z2.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqY2Z2[0] = boostWorld.isend(rankY2Z2, mpiTag,  particleY2Z2);
      reqY2Z2[1] = boostWorld.irecv(rankY2Z2, mpiTag, rParticleY2Z2);
    }
    // 8 vertices
    if (rankX1Y1Z1 >= 0) { // edge x1y1z1
      Rectangle containerX1Y1Z1(v1.getX() - segX, v1.getY() - segY, v1.getZ() - segZ,
        v1.getX(), v1.getY(), v1.getZ());
      findParticleInRectangle(containerX1Y1Z1, particleVec, particleX1Y1Z1);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX1Y1Z1.begin(); it!=particleX1Y1Z1.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX1Y1Z1[0] = boostWorld.isend(rankX1Y1Z1, mpiTag,  particleX1Y1Z1);
      reqX1Y1Z1[1] = boostWorld.irecv(rankX1Y1Z1, mpiTag, rParticleX1Y1Z1);
    }
    if (rankX1Y1Z2 >= 0) { // edge x1y1z2
      Rectangle containerX1Y1Z2(v1.getX() - segX, v1.getY() - segY, v2.getZ(),
        v1.getX(), v1.getY(), v2.getZ() + segZ);
      findParticleInRectangle(containerX1Y1Z2, particleVec, particleX1Y1Z2);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX1Y1Z2.begin(); it!=particleX1Y1Z2.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX1Y1Z2[0] = boostWorld.isend(rankX1Y1Z2, mpiTag,  particleX1Y1Z2);
      reqX1Y1Z2[1] = boostWorld.irecv(rankX1Y1Z2, mpiTag, rParticleX1Y1Z2);
    }
    if (rankX1Y2Z1 >= 0) { // edge x1y2z1
      Rectangle containerX1Y2Z1(v1.getX() - segX, v2.getY(), v1.getZ() - segZ,
        v1.getX(), v2.getY() + segY, v1.getZ());
      findParticleInRectangle(containerX1Y2Z1, particleVec, particleX1Y2Z1);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX1Y2Z1.begin(); it!=particleX1Y2Z1.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX1Y2Z1[0] = boostWorld.isend(rankX1Y2Z1, mpiTag,  particleX1Y2Z1);
      reqX1Y2Z1[1] = boostWorld.irecv(rankX1Y2Z1, mpiTag, rParticleX1Y2Z1);
    }
    if (rankX1Y2Z2 >= 0) { // edge x1y2z2
      Rectangle containerX1Y2Z2(v1.getX() - segX, v2.getY(), v2.getZ(),
        v1.getX(), v2.getY() + segY, v2.getZ() + segZ);
      findParticleInRectangle(containerX1Y2Z2, particleVec, particleX1Y2Z2);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX1Y2Z2.begin(); it!=particleX1Y2Z2.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX1Y2Z2[0] = boostWorld.isend(rankX1Y2Z2, mpiTag,  particleX1Y2Z2);
      reqX1Y2Z2[1] = boostWorld.irecv(rankX1Y2Z2, mpiTag, rParticleX1Y2Z2);
    }
    if (rankX2Y1Z1 >= 0) { // edge x2y1z1
      Rectangle containerX2Y1Z1(v2.getX(), v1.getY() - segY, v1.getZ() - segZ,
        v2.getX() + segX, v1.getY(), v1.getZ());
      findParticleInRectangle(containerX2Y1Z1, particleVec, particleX2Y1Z1);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX2Y1Z1.begin(); it!=particleX2Y1Z1.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX2Y1Z1[0] = boostWorld.isend(rankX2Y1Z1, mpiTag,  particleX2Y1Z1);
      reqX2Y1Z1[1] = boostWorld.irecv(rankX2Y1Z1, mpiTag, rParticleX2Y1Z1);
    }
    if (rankX2Y1Z2 >= 0) { // edge x2y1z2
      Rectangle containerX2Y1Z2(v2.getX(), v1.getY() - segY, v2.getZ(),
        v2.getX() + segX, v1.getY(), v2.getZ() + segZ);
      findParticleInRectangle(containerX2Y1Z2, particleVec, particleX2Y1Z2);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX2Y1Z2.begin(); it!=particleX2Y1Z2.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX2Y1Z2[0] = boostWorld.isend(rankX2Y1Z2, mpiTag,  particleX2Y1Z2);
      reqX2Y1Z2[1] = boostWorld.irecv(rankX2Y1Z2, mpiTag, rParticleX2Y1Z2);
    }
    if (rankX2Y2Z1 >= 0) { // edge x2y2z1
      Rectangle containerX2Y2Z1(v2.getX(), v2.getY(), v1.getZ() - segZ,
        v2.getX() + segX, v2.getY() + segY, v1.getZ());
      findParticleInRectangle(containerX2Y2Z1, particleVec, particleX2Y2Z1);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX2Y2Z1.begin(); it!=particleX2Y2Z1.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX2Y2Z1[0] = boostWorld.isend(rankX2Y2Z1, mpiTag,  particleX2Y2Z1);
      reqX2Y2Z1[1] = boostWorld.irecv(rankX2Y2Z1, mpiTag, rParticleX2Y2Z1);
    }
    if (rankX2Y2Z2 >= 0) { // edge x2y2z2
      Rectangle containerX2Y2Z2(v2.getX(), v2.getY(), v2.getZ(),
        v2.getX() + segX, v2.getY() + segY, v2.getZ() + segZ);
      findParticleInRectangle(containerX2Y2Z2, particleVec, particleX2Y2Z2);
      // before send, SPHParticle.demParticle should be NULL
      for(std::vector<Particle*>::iterator it=particleX2Y2Z2.begin(); it!=particleX2Y2Z2.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();
      reqX2Y2Z2[0] = boostWorld.isend(rankX2Y2Z2, mpiTag,  particleX2Y2Z2);
      reqX2Y2Z2[1] = boostWorld.irecv(rankX2Y2Z2, mpiTag, rParticleX2Y2Z2);
    }
    // 6 surfaces
    if (rankX1 >= 0) boost::mpi::wait_all(reqX1, reqX1 + 2);
    if (rankX2 >= 0) boost::mpi::wait_all(reqX2, reqX2 + 2);
    if (rankY1 >= 0) boost::mpi::wait_all(reqY1, reqY1 + 2);
    if (rankY2 >= 0) boost::mpi::wait_all(reqY2, reqY2 + 2);
    if (rankZ1 >= 0) boost::mpi::wait_all(reqZ1, reqZ1 + 2);
    if (rankZ2 >= 0) boost::mpi::wait_all(reqZ2, reqZ2 + 2);
    // 12 edges
    if (rankX1Y1 >= 0) boost::mpi::wait_all(reqX1Y1, reqX1Y1 + 2);
    if (rankX1Y2 >= 0) boost::mpi::wait_all(reqX1Y2, reqX1Y2 + 2);  
    if (rankX1Z1 >= 0) boost::mpi::wait_all(reqX1Z1, reqX1Z1 + 2);
    if (rankX1Z2 >= 0) boost::mpi::wait_all(reqX1Z2, reqX1Z2 + 2);
    if (rankX2Y1 >= 0) boost::mpi::wait_all(reqX2Y1, reqX2Y1 + 2);
    if (rankX2Y2 >= 0) boost::mpi::wait_all(reqX2Y2, reqX2Y2 + 2);  
    if (rankX2Z1 >= 0) boost::mpi::wait_all(reqX2Z1, reqX2Z1 + 2);
    if (rankX2Z2 >= 0) boost::mpi::wait_all(reqX2Z2, reqX2Z2 + 2); 
    if (rankY1Z1 >= 0) boost::mpi::wait_all(reqY1Z1, reqY1Z1 + 2);
    if (rankY1Z2 >= 0) boost::mpi::wait_all(reqY1Z2, reqY1Z2 + 2);
    if (rankY2Z1 >= 0) boost::mpi::wait_all(reqY2Z1, reqY2Z1 + 2);
    if (rankY2Z2 >= 0) boost::mpi::wait_all(reqY2Z2, reqY2Z2 + 2); 
    // 8 vertices
    if (rankX1Y1Z1 >= 0) boost::mpi::wait_all(reqX1Y1Z1, reqX1Y1Z1 + 2);
    if (rankX1Y1Z2 >= 0) boost::mpi::wait_all(reqX1Y1Z2, reqX1Y1Z2 + 2);
    if (rankX1Y2Z1 >= 0) boost::mpi::wait_all(reqX1Y2Z1, reqX1Y2Z1 + 2);
    if (rankX1Y2Z2 >= 0) boost::mpi::wait_all(reqX1Y2Z2, reqX1Y2Z2 + 2);
    if (rankX2Y1Z1 >= 0) boost::mpi::wait_all(reqX2Y1Z1, reqX2Y1Z1 + 2);
    if (rankX2Y1Z2 >= 0) boost::mpi::wait_all(reqX2Y1Z2, reqX2Y1Z2 + 2);
    if (rankX2Y2Z1 >= 0) boost::mpi::wait_all(reqX2Y2Z1, reqX2Y2Z1 + 2);
    if (rankX2Y2Z2 >= 0) boost::mpi::wait_all(reqX2Y2Z2, reqX2Y2Z2 + 2);  

    // do not need to assign SPHGhostParticle.demParticle back in particleX1..., since these particles will be removed
    // delete outgoing particles
    removeParticleOutRectangle();

    // add incoming particles
    recvParticleVec.clear(); // new use of recvParticleVec
    // 6 surfaces
    if (rankX1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1.begin(), rParticleX1.end());
    if (rankX2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2.begin(), rParticleX2.end());
    if (rankY1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY1.begin(), rParticleY1.end());
    if (rankY2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY2.begin(), rParticleY2.end());
    if (rankZ1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleZ1.begin(), rParticleZ1.end());
    if (rankZ2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleZ2.begin(), rParticleZ2.end());
    // 12 edges
    if (rankX1Y1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1.begin(), rParticleX1Y1.end());
    if (rankX1Y2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2.begin(), rParticleX1Y2.end());
    if (rankX1Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Z1.begin(), rParticleX1Z1.end());
    if (rankX1Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Z2.begin(), rParticleX1Z2.end());
    if (rankX2Y1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1.begin(), rParticleX2Y1.end());
    if (rankX2Y2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2.begin(), rParticleX2Y2.end());
    if (rankX2Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Z1.begin(), rParticleX2Z1.end());
    if (rankX2Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Z2.begin(), rParticleX2Z2.end());
    if (rankY1Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY1Z1.begin(), rParticleY1Z1.end());
    if (rankY1Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY1Z2.begin(), rParticleY1Z2.end());
    if (rankY2Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY2Z1.begin(), rParticleY2Z1.end());
    if (rankY2Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY2Z2.begin(), rParticleY2Z2.end());
    // 8 vertices
    if (rankX1Y1Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1Z1.begin(), rParticleX1Y1Z1.end());
    if (rankX1Y1Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1Z2.begin(), rParticleX1Y1Z2.end());
    if (rankX1Y2Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2Z1.begin(), rParticleX1Y2Z1.end());
    if (rankX1Y2Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2Z2.begin(), rParticleX1Y2Z2.end());
    if (rankX2Y1Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1Z1.begin(), rParticleX2Y1Z1.end());
    if (rankX2Y1Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1Z2.begin(), rParticleX2Y1Z2.end());
    if (rankX2Y2Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2Z1.begin(), rParticleX2Y2Z1.end());
    if (rankX2Y2Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2Z2.begin(), rParticleX2Y2Z2.end());

    // after receive, set SPHParticle.demParticle
    for(std::vector<Particle*>::iterator it=recvParticleVec.begin(); it!=recvParticleVec.end(); it++){
  (*it)->setDemParticleInSPHParticle();
    }

    particleVec.insert(particleVec.end(), recvParticleVec.begin(), recvParticleVec.end());

    /*
      if (recvParticleVec.size() > 0) {    
      debugInf << "iter=" << std::setw(8) << iteration << " rank=" << std::setw(2) << mpiRank 
      << "   added=";
      for (std::vector<Particle*>::const_iterator it = recvParticleVec.begin(); it != recvParticleVec.end(); ++it)
      debugInf << std::setw(3) << (*it)->getId();
      debugInf << " now " << particleVec.size() << ": ";
      for (std::vector<Particle*>::const_iterator it = particleVec.begin(); it != particleVec.end(); ++it)
      debugInf << std::setw(3) << (*it)->getId();
      debugInf << std::endl;
      }
    */

    // do not release memory of received particles because they are part of and managed by particleVec
    // 6 surfaces
    rParticleX1.clear();
    rParticleX2.clear();
    rParticleY1.clear();
    rParticleY2.clear();
    rParticleZ1.clear();
    rParticleZ2.clear();
    // 12 edges
    rParticleX1Y1.clear();
    rParticleX1Y2.clear();
    rParticleX1Z1.clear();
    rParticleX1Z2.clear();
    rParticleX2Y1.clear();
    rParticleX2Y2.clear();
    rParticleX2Z1.clear();
    rParticleX2Z2.clear();
    rParticleY1Z1.clear();
    rParticleY1Z2.clear();
    rParticleY2Z1.clear();
    rParticleY2Z2.clear();
    // 8 vertices
    rParticleX1Y1Z1.clear();
    rParticleX1Y1Z2.clear();
    rParticleX1Y2Z1.clear();
    rParticleX1Y2Z2.clear();
    rParticleX2Y1Z1.clear();
    rParticleX2Y1Z2.clear();
    rParticleX2Y2Z1.clear();
    rParticleX2Y2Z2.clear();

    recvParticleVec.clear();
  }


  void SmoothParticleHydrodynamics::migrateSPHParticle() 
  {   
    Vec vspan = grid.getMaxCorner() - grid.getMinCorner();
    REAL segX = vspan.getX() / mpiProcX;
    REAL segY = vspan.getY() / mpiProcY;
    REAL segZ = vspan.getZ() / mpiProcZ;
    Vec v1 = container.getMinCorner(); // v1, v2 in terms of process
    Vec v2 = container.getMaxCorner();  

    // if a neighbor exists, transfer particles crossing the boundary in between.
    std::vector<sph::SPHParticle*> sphParticleX1, sphParticleX2;
    std::vector<sph::SPHParticle*> sphParticleY1, sphParticleY2;
    std::vector<sph::SPHParticle*> sphParticleZ1, sphParticleZ2;
    std::vector<sph::SPHParticle*> sphParticleX1Y1, sphParticleX1Y2, sphParticleX1Z1, sphParticleX1Z2; 
    std::vector<sph::SPHParticle*> sphParticleX2Y1, sphParticleX2Y2, sphParticleX2Z1, sphParticleX2Z2; 
    std::vector<sph::SPHParticle*> sphParticleY1Z1, sphParticleY1Z2, sphParticleY2Z1, sphParticleY2Z2; 
    std::vector<sph::SPHParticle*> sphParticleX1Y1Z1, sphParticleX1Y1Z2, sphParticleX1Y2Z1, sphParticleX1Y2Z2; 
    std::vector<sph::SPHParticle*> sphParticleX2Y1Z1, sphParticleX2Y1Z2, sphParticleX2Y2Z1, sphParticleX2Y2Z2; 
    boost::mpi::request reqX1[2], reqX2[2];
    boost::mpi::request reqY1[2], reqY2[2];
    boost::mpi::request reqZ1[2], reqZ2[2];
    boost::mpi::request reqX1Y1[2], reqX1Y2[2], reqX1Z1[2], reqX1Z2[2];
    boost::mpi::request reqX2Y1[2], reqX2Y2[2], reqX2Z1[2], reqX2Z2[2];
    boost::mpi::request reqY1Z1[2], reqY1Z2[2], reqY2Z1[2], reqY2Z2[2];
    boost::mpi::request reqX1Y1Z1[2], reqX1Y1Z2[2], reqX1Y2Z1[2], reqX1Y2Z2[2];
    boost::mpi::request reqX2Y1Z1[2], reqX2Y1Z2[2], reqX2Y2Z1[2], reqX2Y2Z2[2];

    // 6 surfaces
    if (rankX1 >= 0) { // surface x1
      Rectangle containerX1(v1.getX() - segX, v1.getY(), v1.getZ(), 
          v1.getX(), v2.getY(), v2.getZ());
      findSPHParticleInRectangle(containerX1, SPHParticleVec, sphParticleX1);
      reqX1[0] = boostWorld.isend(rankX1, mpiTag,  sphParticleX1);
      reqX1[1] = boostWorld.irecv(rankX1, mpiTag, rsphParticleX1);
    }
    if (rankX2 >= 0) { // surface x2
      Rectangle containerX2(v2.getX(), v1.getY(), v1.getZ(),
          v2.getX() + segX, v2.getY(), v2.getZ());
      findSPHParticleInRectangle(containerX2, SPHParticleVec, sphParticleX2);
      reqX2[0] = boostWorld.isend(rankX2, mpiTag,  sphParticleX2);
      reqX2[1] = boostWorld.irecv(rankX2, mpiTag, rsphParticleX2);
    }
    if (rankY1 >= 0) {  // surface y1
      Rectangle containerY1(v1.getX(), v1.getY() - segY, v1.getZ(), 
          v2.getX(), v1.getY(), v2.getZ());
      findSPHParticleInRectangle(containerY1, SPHParticleVec, sphParticleY1);
      reqY1[0] = boostWorld.isend(rankY1, mpiTag,  sphParticleY1);
      reqY1[1] = boostWorld.irecv(rankY1, mpiTag, rsphParticleY1);
    }
    if (rankY2 >= 0) {  // surface y2
      Rectangle containerY2(v1.getX(), v2.getY(), v1.getZ(),
          v2.getX(), v2.getY() + segY, v2.getZ());
      findSPHParticleInRectangle(containerY2, SPHParticleVec, sphParticleY2);
      reqY2[0] = boostWorld.isend(rankY2, mpiTag,  sphParticleY2);
      reqY2[1] = boostWorld.irecv(rankY2, mpiTag, rsphParticleY2);
    }
    if (rankZ1 >= 0) {  // surface z1
      Rectangle containerZ1(v1.getX(), v1.getY(), v1.getZ() - segZ,
          v2.getX(), v2.getY(), v1.getZ());
      findSPHParticleInRectangle(containerZ1, SPHParticleVec, sphParticleZ1);
      reqZ1[0] = boostWorld.isend(rankZ1, mpiTag,  sphParticleZ1);
      reqZ1[1] = boostWorld.irecv(rankZ1, mpiTag, rsphParticleZ1);
    }
    if (rankZ2 >= 0) {  // surface z2
      Rectangle containerZ2(v1.getX(), v1.getY(), v2.getZ(),
          v2.getX(), v2.getY(), v2.getZ() + segZ);
      findSPHParticleInRectangle(containerZ2, SPHParticleVec, sphParticleZ2);
      reqZ2[0] = boostWorld.isend(rankZ2, mpiTag,  sphParticleZ2);
      reqZ2[1] = boostWorld.irecv(rankZ2, mpiTag, rsphParticleZ2);
    }
    // 12 edges
    if (rankX1Y1 >= 0) { // edge x1y1
      Rectangle containerX1Y1(v1.getX() - segX, v1.getY() - segY, v1.getZ(),
            v1.getX(), v1.getY(), v2.getZ());
      findSPHParticleInRectangle(containerX1Y1, SPHParticleVec, sphParticleX1Y1);
      reqX1Y1[0] = boostWorld.isend(rankX1Y1, mpiTag,  sphParticleX1Y1);
      reqX1Y1[1] = boostWorld.irecv(rankX1Y1, mpiTag, rsphParticleX1Y1);
    }
    if (rankX1Y2 >= 0) { // edge x1y2
      Rectangle containerX1Y2(v1.getX() - segX, v2.getY(), v1.getZ(),
            v1.getX(), v2.getY() + segY, v2.getZ());
      findSPHParticleInRectangle(containerX1Y2, SPHParticleVec, sphParticleX1Y2);
      reqX1Y2[0] = boostWorld.isend(rankX1Y2, mpiTag,  sphParticleX1Y2);
      reqX1Y2[1] = boostWorld.irecv(rankX1Y2, mpiTag, rsphParticleX1Y2);
    }
    if (rankX1Z1 >= 0) { // edge x1z1
      Rectangle containerX1Z1(v1.getX() - segX, v1.getY(), v1.getZ() -segZ,
            v1.getX(), v2.getY(), v1.getZ());
      findSPHParticleInRectangle(containerX1Z1, SPHParticleVec, sphParticleX1Z1);
      reqX1Z1[0] = boostWorld.isend(rankX1Z1, mpiTag,  sphParticleX1Z1);
      reqX1Z1[1] = boostWorld.irecv(rankX1Z1, mpiTag, rsphParticleX1Z1);
    }
    if (rankX1Z2 >= 0) { // edge x1z2
      Rectangle containerX1Z2(v1.getX() - segX, v1.getY(), v2.getZ(),
            v1.getX(), v2.getY(), v2.getZ() + segZ);
      findSPHParticleInRectangle(containerX1Z2, SPHParticleVec, sphParticleX1Z2);
      reqX1Z2[0] = boostWorld.isend(rankX1Z2, mpiTag,  sphParticleX1Z2);
      reqX1Z2[1] = boostWorld.irecv(rankX1Z2, mpiTag, rsphParticleX1Z2);
    }
    if (rankX2Y1 >= 0) { // edge x2y1
      Rectangle containerX2Y1(v2.getX(), v1.getY() - segY, v1.getZ(),
            v2.getX() + segX, v1.getY(), v2.getZ());
      findSPHParticleInRectangle(containerX2Y1, SPHParticleVec, sphParticleX2Y1);
      reqX2Y1[0] = boostWorld.isend(rankX2Y1, mpiTag,  sphParticleX2Y1);
      reqX2Y1[1] = boostWorld.irecv(rankX2Y1, mpiTag, rsphParticleX2Y1);
    }
    if (rankX2Y2 >= 0) { // edge x2y2
      Rectangle containerX2Y2(v2.getX(), v2.getY(), v1.getZ(),
            v2.getX() + segX, v2.getY() + segY, v2.getZ());
      findSPHParticleInRectangle(containerX2Y2, SPHParticleVec, sphParticleX2Y2);
      reqX2Y2[0] = boostWorld.isend(rankX2Y2, mpiTag,  sphParticleX2Y2);
      reqX2Y2[1] = boostWorld.irecv(rankX2Y2, mpiTag, rsphParticleX2Y2);
    }
    if (rankX2Z1 >= 0) { // edge x2z1
      Rectangle containerX2Z1(v2.getX(), v1.getY(), v1.getZ() - segZ,
            v2.getX() + segX, v2.getY(), v1.getZ());
      findSPHParticleInRectangle(containerX2Z1, SPHParticleVec, sphParticleX2Z1);
      reqX2Z1[0] = boostWorld.isend(rankX2Z1, mpiTag,  sphParticleX2Z1);
      reqX2Z1[1] = boostWorld.irecv(rankX2Z1, mpiTag, rsphParticleX2Z1);
    }
    if (rankX2Z2 >= 0) { // edge x2z2
      Rectangle containerX2Z2(v2.getX(), v1.getY(), v2.getZ(),
            v2.getX() + segX, v2.getY(), v2.getZ() + segZ);
      findSPHParticleInRectangle(containerX2Z2, SPHParticleVec, sphParticleX2Z2);
      reqX2Z2[0] = boostWorld.isend(rankX2Z2, mpiTag,  sphParticleX2Z2);
      reqX2Z2[1] = boostWorld.irecv(rankX2Z2, mpiTag, rsphParticleX2Z2);
    }
    if (rankY1Z1 >= 0) { // edge y1z1
      Rectangle containerY1Z1(v1.getX(), v1.getY() - segY, v1.getZ() - segZ,
            v2.getX(), v1.getY(), v1.getZ());
      findSPHParticleInRectangle(containerY1Z1, SPHParticleVec, sphParticleY1Z1);
      reqY1Z1[0] = boostWorld.isend(rankY1Z1, mpiTag,  sphParticleY1Z1);
      reqY1Z1[1] = boostWorld.irecv(rankY1Z1, mpiTag, rsphParticleY1Z1);
    }
    if (rankY1Z2 >= 0) { // edge y1z2
      Rectangle containerY1Z2(v1.getX(), v1.getY() - segY, v2.getZ(),
            v2.getX(), v1.getY(), v2.getZ() + segZ);
      findSPHParticleInRectangle(containerY1Z2, SPHParticleVec, sphParticleY1Z2);
      reqY1Z2[0] = boostWorld.isend(rankY1Z2, mpiTag,  sphParticleY1Z2);
      reqY1Z2[1] = boostWorld.irecv(rankY1Z2, mpiTag, rsphParticleY1Z2);
    }
    if (rankY2Z1 >= 0) { // edge y2z1
      Rectangle containerY2Z1(v1.getX(), v2.getY(), v1.getZ() - segZ,
            v2.getX(), v2.getY() + segY, v1.getZ());
      findSPHParticleInRectangle(containerY2Z1, SPHParticleVec, sphParticleY2Z1);
      reqY2Z1[0] = boostWorld.isend(rankY2Z1, mpiTag,  sphParticleY2Z1);
      reqY2Z1[1] = boostWorld.irecv(rankY2Z1, mpiTag, rsphParticleY2Z1);
    }
    if (rankY2Z2 >= 0) { // edge y2z2
      Rectangle containerY2Z2(v1.getX(), v2.getY(), v2.getZ(),
            v2.getX(), v2.getY() + segY, v2.getZ() + segZ);
      findSPHParticleInRectangle(containerY2Z2, SPHParticleVec, sphParticleY2Z2);
      reqY2Z2[0] = boostWorld.isend(rankY2Z2, mpiTag,  sphParticleY2Z2);
      reqY2Z2[1] = boostWorld.irecv(rankY2Z2, mpiTag, rsphParticleY2Z2);
    }
    // 8 vertices
    if (rankX1Y1Z1 >= 0) { // edge x1y1z1
      Rectangle containerX1Y1Z1(v1.getX() - segX, v1.getY() - segY, v1.getZ() - segZ,
        v1.getX(), v1.getY(), v1.getZ());
      findSPHParticleInRectangle(containerX1Y1Z1, SPHParticleVec, sphParticleX1Y1Z1);
      reqX1Y1Z1[0] = boostWorld.isend(rankX1Y1Z1, mpiTag,  sphParticleX1Y1Z1);
      reqX1Y1Z1[1] = boostWorld.irecv(rankX1Y1Z1, mpiTag, rsphParticleX1Y1Z1);
    }
    if (rankX1Y1Z2 >= 0) { // edge x1y1z2
      Rectangle containerX1Y1Z2(v1.getX() - segX, v1.getY() - segY, v2.getZ(),
        v1.getX(), v1.getY(), v2.getZ() + segZ);
      findSPHParticleInRectangle(containerX1Y1Z2, SPHParticleVec, sphParticleX1Y1Z2);
      reqX1Y1Z2[0] = boostWorld.isend(rankX1Y1Z2, mpiTag,  sphParticleX1Y1Z2);
      reqX1Y1Z2[1] = boostWorld.irecv(rankX1Y1Z2, mpiTag, rsphParticleX1Y1Z2);
    }
    if (rankX1Y2Z1 >= 0) { // edge x1y2z1
      Rectangle containerX1Y2Z1(v1.getX() - segX, v2.getY(), v1.getZ() - segZ,
        v1.getX(), v2.getY() + segY, v1.getZ());
      findSPHParticleInRectangle(containerX1Y2Z1, SPHParticleVec, sphParticleX1Y2Z1);
      reqX1Y2Z1[0] = boostWorld.isend(rankX1Y2Z1, mpiTag,  sphParticleX1Y2Z1);
      reqX1Y2Z1[1] = boostWorld.irecv(rankX1Y2Z1, mpiTag, rsphParticleX1Y2Z1);
    }
    if (rankX1Y2Z2 >= 0) { // edge x1y2z2
      Rectangle containerX1Y2Z2(v1.getX() - segX, v2.getY(), v2.getZ(),
        v1.getX(), v2.getY() + segY, v2.getZ() + segZ);
      findSPHParticleInRectangle(containerX1Y2Z2, SPHParticleVec, sphParticleX1Y2Z2);
      reqX1Y2Z2[0] = boostWorld.isend(rankX1Y2Z2, mpiTag,  sphParticleX1Y2Z2);
      reqX1Y2Z2[1] = boostWorld.irecv(rankX1Y2Z2, mpiTag, rsphParticleX1Y2Z2);
    }
    if (rankX2Y1Z1 >= 0) { // edge x2y1z1
      Rectangle containerX2Y1Z1(v2.getX(), v1.getY() - segY, v1.getZ() - segZ,
        v2.getX() + segX, v1.getY(), v1.getZ());
      findSPHParticleInRectangle(containerX2Y1Z1, SPHParticleVec, sphParticleX2Y1Z1);
      reqX2Y1Z1[0] = boostWorld.isend(rankX2Y1Z1, mpiTag,  sphParticleX2Y1Z1);
      reqX2Y1Z1[1] = boostWorld.irecv(rankX2Y1Z1, mpiTag, rsphParticleX2Y1Z1);
    }
    if (rankX2Y1Z2 >= 0) { // edge x2y1z2
      Rectangle containerX2Y1Z2(v2.getX(), v1.getY() - segY, v2.getZ(),
        v2.getX() + segX, v1.getY(), v2.getZ() + segZ);
      findSPHParticleInRectangle(containerX2Y1Z2, SPHParticleVec, sphParticleX2Y1Z2);
      reqX2Y1Z2[0] = boostWorld.isend(rankX2Y1Z2, mpiTag,  sphParticleX2Y1Z2);
      reqX2Y1Z2[1] = boostWorld.irecv(rankX2Y1Z2, mpiTag, rsphParticleX2Y1Z2);
    }
    if (rankX2Y2Z1 >= 0) { // edge x2y2z1
      Rectangle containerX2Y2Z1(v2.getX(), v2.getY(), v1.getZ() - segZ,
        v2.getX() + segX, v2.getY() + segY, v1.getZ());
      findSPHParticleInRectangle(containerX2Y2Z1, SPHParticleVec, sphParticleX2Y2Z1);
      reqX2Y2Z1[0] = boostWorld.isend(rankX2Y2Z1, mpiTag,  sphParticleX2Y2Z1);
      reqX2Y2Z1[1] = boostWorld.irecv(rankX2Y2Z1, mpiTag, rsphParticleX2Y2Z1);
    }
    if (rankX2Y2Z2 >= 0) { // edge x2y2z2
      Rectangle containerX2Y2Z2(v2.getX(), v2.getY(), v2.getZ(),
        v2.getX() + segX, v2.getY() + segY, v2.getZ() + segZ);
      findSPHParticleInRectangle(containerX2Y2Z2, SPHParticleVec, sphParticleX2Y2Z2);
      reqX2Y2Z2[0] = boostWorld.isend(rankX2Y2Z2, mpiTag,  sphParticleX2Y2Z2);
      reqX2Y2Z2[1] = boostWorld.irecv(rankX2Y2Z2, mpiTag, rsphParticleX2Y2Z2);
    }
    // 6 surfaces
    if (rankX1 >= 0) boost::mpi::wait_all(reqX1, reqX1 + 2);
    if (rankX2 >= 0) boost::mpi::wait_all(reqX2, reqX2 + 2);
    if (rankY1 >= 0) boost::mpi::wait_all(reqY1, reqY1 + 2);
    if (rankY2 >= 0) boost::mpi::wait_all(reqY2, reqY2 + 2);
    if (rankZ1 >= 0) boost::mpi::wait_all(reqZ1, reqZ1 + 2);
    if (rankZ2 >= 0) boost::mpi::wait_all(reqZ2, reqZ2 + 2);
    // 12 edges
    if (rankX1Y1 >= 0) boost::mpi::wait_all(reqX1Y1, reqX1Y1 + 2);
    if (rankX1Y2 >= 0) boost::mpi::wait_all(reqX1Y2, reqX1Y2 + 2);  
    if (rankX1Z1 >= 0) boost::mpi::wait_all(reqX1Z1, reqX1Z1 + 2);
    if (rankX1Z2 >= 0) boost::mpi::wait_all(reqX1Z2, reqX1Z2 + 2);
    if (rankX2Y1 >= 0) boost::mpi::wait_all(reqX2Y1, reqX2Y1 + 2);
    if (rankX2Y2 >= 0) boost::mpi::wait_all(reqX2Y2, reqX2Y2 + 2);  
    if (rankX2Z1 >= 0) boost::mpi::wait_all(reqX2Z1, reqX2Z1 + 2);
    if (rankX2Z2 >= 0) boost::mpi::wait_all(reqX2Z2, reqX2Z2 + 2); 
    if (rankY1Z1 >= 0) boost::mpi::wait_all(reqY1Z1, reqY1Z1 + 2);
    if (rankY1Z2 >= 0) boost::mpi::wait_all(reqY1Z2, reqY1Z2 + 2);
    if (rankY2Z1 >= 0) boost::mpi::wait_all(reqY2Z1, reqY2Z1 + 2);
    if (rankY2Z2 >= 0) boost::mpi::wait_all(reqY2Z2, reqY2Z2 + 2); 
    // 8 vertices
    if (rankX1Y1Z1 >= 0) boost::mpi::wait_all(reqX1Y1Z1, reqX1Y1Z1 + 2);
    if (rankX1Y1Z2 >= 0) boost::mpi::wait_all(reqX1Y1Z2, reqX1Y1Z2 + 2);
    if (rankX1Y2Z1 >= 0) boost::mpi::wait_all(reqX1Y2Z1, reqX1Y2Z1 + 2);
    if (rankX1Y2Z2 >= 0) boost::mpi::wait_all(reqX1Y2Z2, reqX1Y2Z2 + 2);
    if (rankX2Y1Z1 >= 0) boost::mpi::wait_all(reqX2Y1Z1, reqX2Y1Z1 + 2);
    if (rankX2Y1Z2 >= 0) boost::mpi::wait_all(reqX2Y1Z2, reqX2Y1Z2 + 2);
    if (rankX2Y2Z1 >= 0) boost::mpi::wait_all(reqX2Y2Z1, reqX2Y2Z1 + 2);
    if (rankX2Y2Z2 >= 0) boost::mpi::wait_all(reqX2Y2Z2, reqX2Y2Z2 + 2);  

    // delete outgoing particles
    removeSPHParticleOutRectangle();

    // add incoming particles
    recvSPHParticleVec.clear(); // new use of recvParticleVec
    // 6 surfaces
    if (rankX1 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX1.begin(), rsphParticleX1.end());
    if (rankX2 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX2.begin(), rsphParticleX2.end());
    if (rankY1 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleY1.begin(), rsphParticleY1.end());
    if (rankY2 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleY2.begin(), rsphParticleY2.end());
    if (rankZ1 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleZ1.begin(), rsphParticleZ1.end());
    if (rankZ2 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleZ2.begin(), rsphParticleZ2.end());
    // 12 edges
    if (rankX1Y1 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX1Y1.begin(), rsphParticleX1Y1.end());
    if (rankX1Y2 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX1Y2.begin(), rsphParticleX1Y2.end());
    if (rankX1Z1 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX1Z1.begin(), rsphParticleX1Z1.end());
    if (rankX1Z2 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX1Z2.begin(), rsphParticleX1Z2.end());
    if (rankX2Y1 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX2Y1.begin(), rsphParticleX2Y1.end());
    if (rankX2Y2 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX2Y2.begin(), rsphParticleX2Y2.end());
    if (rankX2Z1 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX2Z1.begin(), rsphParticleX2Z1.end());
    if (rankX2Z2 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX2Z2.begin(), rsphParticleX2Z2.end());
    if (rankY1Z1 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleY1Z1.begin(), rsphParticleY1Z1.end());
    if (rankY1Z2 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleY1Z2.begin(), rsphParticleY1Z2.end());
    if (rankY2Z1 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleY2Z1.begin(), rsphParticleY2Z1.end());
    if (rankY2Z2 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleY2Z2.begin(), rsphParticleY2Z2.end());
    // 8 vertices
    if (rankX1Y1Z1 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX1Y1Z1.begin(), rsphParticleX1Y1Z1.end());
    if (rankX1Y1Z2 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX1Y1Z2.begin(), rsphParticleX1Y1Z2.end());
    if (rankX1Y2Z1 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX1Y2Z1.begin(), rsphParticleX1Y2Z1.end());
    if (rankX1Y2Z2 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX1Y2Z2.begin(), rsphParticleX1Y2Z2.end());
    if (rankX2Y1Z1 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX2Y1Z1.begin(), rsphParticleX2Y1Z1.end());
    if (rankX2Y1Z2 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX2Y1Z2.begin(), rsphParticleX2Y1Z2.end());
    if (rankX2Y2Z1 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX2Y2Z1.begin(), rsphParticleX2Y2Z1.end());
    if (rankX2Y2Z2 >= 0) recvSPHParticleVec.insert(recvSPHParticleVec.end(), rsphParticleX2Y2Z2.begin(), rsphParticleX2Y2Z2.end());

    SPHParticleVec.insert(SPHParticleVec.end(), recvSPHParticleVec.begin(), recvSPHParticleVec.end());

    /*
      if (recvParticleVec.size() > 0) {    
      debugInf << "iter=" << std::setw(8) << iteration << " rank=" << std::setw(2) << mpiRank 
      << "   added=";
      for (std::vector<Particle*>::const_iterator it = recvParticleVec.begin(); it != recvParticleVec.end(); ++it)
      debugInf << std::setw(3) << (*it)->getId();
      debugInf << " now " << particleVec.size() << ": ";
      for (std::vector<Particle*>::const_iterator it = particleVec.begin(); it != particleVec.end(); ++it)
      debugInf << std::setw(3) << (*it)->getId();
      debugInf << std::endl;
      }
    */

    // do not release memory of received particles because they are part of and managed by particleVec
    // 6 surfaces
    rsphParticleX1.clear();
    rsphParticleX2.clear();
    rsphParticleY1.clear();
    rsphParticleY2.clear();
    rsphParticleZ1.clear();
    rsphParticleZ2.clear();
    // 12 edges
    rsphParticleX1Y1.clear();
    rsphParticleX1Y2.clear();
    rsphParticleX1Z1.clear();
    rsphParticleX1Z2.clear();
    rsphParticleX2Y1.clear();
    rsphParticleX2Y2.clear();
    rsphParticleX2Z1.clear();
    rsphParticleX2Z2.clear();
    rsphParticleY1Z1.clear();
    rsphParticleY1Z2.clear();
    rsphParticleY2Z1.clear();
    rsphParticleY2Z2.clear();
    // 8 vertices
    rsphParticleX1Y1Z1.clear();
    rsphParticleX1Y1Z2.clear();
    rsphParticleX1Y2Z1.clear();
    rsphParticleX1Y2Z2.clear();
    rsphParticleX2Y1Z1.clear();
    rsphParticleX2Y1Z2.clear();
    rsphParticleX2Y2Z1.clear();
    rsphParticleX2Y2Z2.clear();

    recvSPHParticleVec.clear();
  }

  void SmoothParticleHydrodynamics::gatherParticle() {
    // before send, SPHParticle.demParticle should be NULL
    for(std::vector<Particle*>::iterator it=particleVec.begin(); it!=particleVec.end(); it++)
  (*it)->setNULLDemParticleInSPHParticle();  // at this point, SPHGhostParticle.demParticle is not pointing to particleVec

    // update allParticleVec: process 0 collects all updated particles from each other process  
    if (mpiRank != 0) {// each process except 0
      boostWorld.send(0, mpiTag, particleVec);
    }
    else { // process 0
      // allParticleVec is cleared before filling with new data
      releaseGatheredParticle();
  
      // duplicate particleVec so that it is not destroyed by allParticleVec in next iteration,
      // otherwise it causes memory error.
      std::vector<Particle*> dupParticleVec(particleVec.size());
      for (std::size_t i = 0; i < dupParticleVec.size(); ++i){
  dupParticleVec[i] = new Particle(*particleVec[i]);  // at this point, dupParticleVec and particleVec are pointint to the same SPHGhoastParticle
  dupParticleVec[i]->SPHGhostParticleVec.clear();    // at this point, dupParticleVec is pointing to nothing
  for(std::vector<sph::SPHParticle*>::iterator st=particleVec[i]->SPHGhostParticleVec.begin(); st!=particleVec[i]->SPHGhostParticleVec.end(); st++){
    sph::SPHParticle* tmp_sph = new sph::SPHParticle(**st);  // create a new SPHGhost particle, which is the same as the one in particleVec
    dupParticleVec[i]->SPHGhostParticleVec.push_back(tmp_sph);  // now dupParticleVec points to the new SPHGhostParticle
  }
      }

      // fill allParticleVec with dupParticleVec and received particles
      allParticleVec.insert(allParticleVec.end(), dupParticleVec.begin(), dupParticleVec.end());

      std::vector<Particle*> tmpParticleVec;
      long gatherRam = 0;
      for (int iRank = 1; iRank < mpiSize; ++iRank) {

  tmpParticleVec.clear();// do not destroy particles!
  boostWorld.recv(iRank, mpiTag, tmpParticleVec);
  allParticleVec.insert(allParticleVec.end(), tmpParticleVec.begin(), tmpParticleVec.end());
  gatherRam += tmpParticleVec.size();

      }
      //debugInf << "gather: particleNum = " << gatherRam <<  " particleRam = " << gatherRam * sizeof(Particle) << std::endl;
    }
    // after receive, set SPHParticle.demParticle
    for(std::vector<Particle*>::iterator it=particleVec.begin(); it!=particleVec.end(); it++){
  (*it)->setDemParticleInSPHParticle();
    }
  }


  void SmoothParticleHydrodynamics::releaseGatheredParticle() {
    // clear allParticleVec, avoid long time memory footprint.
    for (std::vector<Particle*>::iterator it = allParticleVec.begin(); it != allParticleVec.end(); ++it){
      for(std::vector<sph::SPHParticle*>::iterator st=(*it)->SPHGhostParticleVec.begin(); st!=(*it)->SPHGhostParticleVec.end(); ++st){
  delete (*st);  // this is important to free the memories of sph ghost particles
      }
      (*it)->SPHGhostParticleVec.clear();
      std::vector<sph::SPHParticle*>().swap((*it)->SPHGhostParticleVec); // actual memory release
      delete (*it);
    }
    allParticleVec.clear();
    std::vector<Particle*>().swap(allParticleVec); // actual memory release
  }


  void SmoothParticleHydrodynamics::gatherSPHParticle() {
    // update allSPHParticleVec: process 0 collects all updated particles from each other process  
    if (mpiRank != 0) {// each process except 0
      boostWorld.send(0, mpiTag, SPHParticleVec);
    }
    else { // process 0
      // allSPHParticleVec is cleared before filling with new data
      releaseGatheredSPHParticle();

      // duplicate SPHParticleVec so that it is not destroyed by allSPHParticleVec in next iteration,
      // otherwise it causes memory error.
      std::vector<sph::SPHParticle*> dupSPHParticleVec(SPHParticleVec.size());
      for (std::size_t i = 0; i < dupSPHParticleVec.size(); ++i)
  dupSPHParticleVec[i] = new sph::SPHParticle(*SPHParticleVec[i]);

      // fill allParticleVec with dupParticleVec and received particles
      allSPHParticleVec.insert(allSPHParticleVec.end(), dupSPHParticleVec.begin(), dupSPHParticleVec.end());

      std::vector<sph::SPHParticle*> tmpSPHParticleVec;
      long gatherRam = 0;
      for (int iRank = 1; iRank < mpiSize; ++iRank) {

  tmpSPHParticleVec.clear();// do not destroy particles!
  boostWorld.recv(iRank, mpiTag, tmpSPHParticleVec);
  allSPHParticleVec.insert(allSPHParticleVec.end(), tmpSPHParticleVec.begin(), tmpSPHParticleVec.end());
  gatherRam += tmpSPHParticleVec.size();

      }
      //debugInf << "gather: particleNum = " << gatherRam <<  " particleRam = " << gatherRam * sizeof(Particle) << std::endl;
    }
  }

  void SmoothParticleHydrodynamics::releaseGatheredSPHParticle() {
    // clear allParticleVec, avoid long time memory footprint.
    for (std::vector<sph::SPHParticle*>::iterator it = allSPHParticleVec.begin(); it != allSPHParticleVec.end(); ++it)
      delete (*it);
    allSPHParticleVec.clear();
    std::vector<sph::SPHParticle*>().swap(allSPHParticleVec); // actual memory release
  }

  void SmoothParticleHydrodynamics::gatherBdryContact() {
    if (isBdryProcess()) {
      if (mpiRank != 0)
  boostWorld.send(0, mpiTag, boundaryVec);
    }

    if (mpiRank == 0) {
      mergeBoundaryVec.clear();
      std::vector<Boundary*>().swap(mergeBoundaryVec); // actual memory release
      mergeBoundaryVec = boundaryVec; 

      std::vector<Boundary*> tmpBoundaryVec;   
      for (std::size_t it = 0; it < bdryProcess.size(); ++it) {
  if (bdryProcess[it] != 0) {// not root process
    tmpBoundaryVec.clear();  // do not destroy particles!
    boostWorld.recv(bdryProcess[it], mpiTag, tmpBoundaryVec);
    // merge tmpBoundaryVec into mergeBoundaryVec
    assert(tmpBoundaryVec.size() == mergeBoundaryVec.size());
    for (std::size_t jt = 0; jt < tmpBoundaryVec.size(); ++jt)
      mergeBoundaryVec[jt]->getContactInfo().insert(   \
                mergeBoundaryVec[jt]->getContactInfo().end(), \
                tmpBoundaryVec[jt]->getContactInfo().begin(), \
                tmpBoundaryVec[jt]->getContactInfo().end() );
  }    
      }

      // must update after collecting all boundary contact info
      for(std::vector<Boundary*>::iterator it = mergeBoundaryVec.begin(); it != mergeBoundaryVec.end(); ++it)
  (*it)->updateStatForce();
    }
  }
  

  void  SmoothParticleHydrodynamics::openSPHTecplot(std::ofstream &ofs, const char *str) {
    ofs.open(str);
    if(!ofs) { debugInf << "stream error: openSPHTecplot" << std::endl; exit(-1); }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);

    ofs << "Title = \"SPH Particle Information\"" << std::endl;
    ofs << "VARIABLES = \"x\", \"y\",\"z\" \"Ux\" \"Uy\" \"Uz\" \"Vx\" \"Vy\" \"Vz\" \"Pressure\" \"a_x\" \"a_y\" \"a_z\" \"density_dot\" \"density\" "
  << std::endl;
  }
  
  void SmoothParticleHydrodynamics::printSPHTecplot(std::ofstream &ofs, int iframe) {
  ofs << "ZONE T =\" " << iframe << "-th Load Step\" "<< std::endl;
  // Output the coordinates and the array information
  for(std::vector<sph::SPHParticle*>::iterator pt = allSPHParticleVec.begin(); pt!= allSPHParticleVec.end(); pt++) {

//// not print the most right layer of SPH free particles, 2015/05/19
//if((*pt)->getInitPosition().getx()==25){
//continue;
//}
//      if((*pt)->getType()==3) continue;  // not print out boundary particles

      ofs << std::setw(20) << (*pt)->getCurrPosition().getX()
    << std::setw(20) << (*pt)->getCurrPosition().getY()
    << std::setw(20) << (*pt)->getCurrPosition().getZ()
    << std::setw(20) << (*pt)->getDisplacement().getX()
    << std::setw(20) << (*pt)->getDisplacement().getY() 
    << std::setw(20) << (*pt)->getDisplacement().getZ() 
    << std::setw(20) << (*pt)->getVelocity().getX()
    << std::setw(20) << (*pt)->getVelocity().getY() 
    << std::setw(20) << (*pt)->getVelocity().getZ();

      (*pt)->calculateParticlePressure();
       ofs << std::setw(20) <<(*pt)->getParticlePressure();

      ofs << std::setw(20) << (*pt)->getVelocityDot().getX()
    << std::setw(20) << (*pt)->getVelocityDot().getY() 
    << std::setw(20) << (*pt)->getVelocityDot().getZ()
    << std::setw(20) << (*pt)->getDensityDot() 
    << std::setw(20) << (*pt)->getParticleDensity() << std::endl;

  }
  }
  void SmoothParticleHydrodynamics::printSPHParticle(const char *str) const{
      std::ofstream ofs(str);
      if(!ofs) { debugInf << "stream error: printSPHParticle" << std::endl; exit(-1); }
      ofs.setf(std::ios::scientific, std::ios::floatfield);
      ofs.precision(OPREC);

        ofs << "Title = \"SPH Particle Information\"" << std::endl;
        ofs << "VARIABLES = \"x\", \"y\",\"z\" \"Ux\" \"Uy\" \"Uz\" \"Vx\" \"Vy\" \"Vz\" \"Pressure\" \"a_x\" \"a_y\" \"a_z\" \"density_dot\" \"density\" "
      << std::endl;
 
  // Output the coordinates and the array information
  for(std::vector<sph::SPHParticle*>::const_iterator pt = allSPHParticleVec.begin(); pt!= allSPHParticleVec.end(); pt++) {

//// not print the most right layer of SPH free particles, 2015/05/19
//if((*pt)->getInitPosition().getx()==25){
//continue;
//}
      if((*pt)->getType()==3) continue;  // not print out boundary particles

      ofs << std::setw(20) << (*pt)->getCurrPosition().getX()
    << std::setw(20) << (*pt)->getCurrPosition().getY()
    << std::setw(20) << (*pt)->getCurrPosition().getZ()
    << std::setw(20) << (*pt)->getDisplacement().getX()
    << std::setw(20) << (*pt)->getDisplacement().getY() 
    << std::setw(20) << (*pt)->getDisplacement().getZ() 
    << std::setw(20) << (*pt)->getVelocity().getX()
    << std::setw(20) << (*pt)->getVelocity().getY() 
    << std::setw(20) << (*pt)->getVelocity().getZ();

      (*pt)->calculateParticlePressure();
       ofs << std::setw(20) <<(*pt)->getParticlePressure();

      ofs << std::setw(20) << (*pt)->getVelocityDot().getX()
    << std::setw(20) << (*pt)->getVelocityDot().getY() 
    << std::setw(20) << (*pt)->getVelocityDot().getZ()
    << std::setw(20) << (*pt)->getDensityDot() 
    << std::setw(20) << (*pt)->getParticleDensity() << std::endl;

  }
  
    ofs.close();
  }

    void SmoothParticleHydrodynamics::initialSPHVelocity2D(){
  commuParticle();
  commuSPHParticle();  // this will update container and mergeSPHParticleVec, both are needed for divideSPHDomain
  calculateSPHDensityDotVelocityDotLinkedList2D();  // calculate velocityDot and densityDot
  // fix the y for all SPH points
  for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
      switch((*pt)->getType()){
        case 1: // free sph particle
    (*pt)->fixY();
    break;
        case 2: // ghost sph particle
    break;
        case 3: // boundary sph particle
    (*pt)->fixXYZ();
    break;
        default:
    std::cout << "SPH particle type of pta should be 1, 2 or 3!" << std::endl;
    exit(-1);
      } // switch
  }
  initialSPHLeapFrogVelocity();  // initial velocity only for free SPH particles based on equation (4.3)

        releaseRecvParticle(); 
        releaseRecvSPHParticle(); 
    } // initialSPHVelocity()

    void SmoothParticleHydrodynamics::initialSPHVelocity3D(){
  commuParticle();
  commuSPHParticle();  // this will update container and mergeSPHParticleVec, both are needed for divideSPHDomain
  calculateSPHDensityDotVelocityDotLinkedList3D();  // calculate velocityDot and densityDot
  initialSPHLeapFrogVelocity();  // initial velocity only for free SPH particles based on equation (4.3)

        releaseRecvParticle(); 
        releaseRecvSPHParticle();
    } // initialSPHVelocity()

    void SmoothParticleHydrodynamics::initialSPHVelocityCopyDEM3D(){
  commuSPHParticle();  // this will update container and mergeSPHParticleVec, both are needed for divideSPHDomain
  calculateSPHDensityDotVelocityDotLinkedList3D();  // calculate velocityDot and densityDot
  initialSPHLeapFrogVelocity();  // initial velocity only for free SPH particles based on equation (4.3)

        releaseRecvSPHParticle();
    } // initialSPHVelocity()


    void SmoothParticleHydrodynamics::initialSPHLeapFrogVelocity(){  // initial particle velocity only for free SPH particles based on equation (4.3)
  for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
      if((*pt)->getType() == 1) (*pt)->initialParticleVelocityLeapFrog();
  }
    } // end initialSPHLeapFrogVelocity

    void SmoothParticleHydrodynamics::updateSPHLeapFrogPositionDensity(){  // update particle position and density based on equation (4.1)
  for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
      switch((*pt)->getType()){
        case 1: // free sph particle
    (*pt)->updateParticlePositionDensityLeapFrog();
    break;
        case 2: // ghost sph particle
    break;
        case 3: // boundary sph particle
    (*pt)->updateParticleDensity();
    break;
        default:
    std::cout << "SPH particle type of pta should be 1, 2 or 3!" << std::endl;
    exit(-1);
      } // switch
  }

  std::vector<sph::SPHParticle*>::iterator gt;
  for(std::vector<Particle*>::iterator demt=particleVec.begin(); demt!=particleVec.end(); demt++){
      (*demt)->updateSPHGhostParticle();  // update position and density by dem particle and sph particle density, and velocity
  }
    } // end updateSPHLeapFrogPositionDensity


    void SmoothParticleHydrodynamics::updateSPHLeapFrogVelocity(){  // update particle velocity only for free SPH particles based on equation (4.2)
  for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
      if((*pt)->getType()==1) (*pt)->updateParticleVelocityLeapFrog();
  }
    } // end updateSPHLeapFrogVelocity


//  REAL factor = 1.0/(120.0*dem::PI*h*h*h);  // 3D quintic kernel factor
//  REAL factor = 7.0/(478.0*dem::PI*h*h);    // 2D quintic kernel factor
    // kernel function, vec is the position of a b, h is smoothing length  
    inline REAL SmoothParticleHydrodynamics::kernelFunction(const dem::Vec& a, const dem::Vec& b){
  REAL rab = dem::vfabs(a-b);
  REAL s = rab*one_devide_h;
  REAL item1_5 = (3-s)*(3-s)*(3-s)*(3-s)*(3-s);  // (3-s)^5
  REAL item2_5 = (2-s)*(2-s)*(2-s)*(2-s)*(2-s);   // (2-s)^5
  REAL item3_5 = (1-s)*(1-s)*(1-s)*(1-s)*(1-s);   // (1-s)^5
  if(s<1){
      return factor_kernel*(item1_5-6*item2_5+15*item3_5);
  }
  else if(s<2){
      return factor_kernel*(item1_5-6*item2_5);
  }
  else if(s<3){
      return factor_kernel*item1_5;
  }
  else{
      return 0.0;
  }

    } // end kernelFunction

    // kernel function, s is rab/h, h is smoothing length  
    inline REAL SmoothParticleHydrodynamics::kernelFunction(REAL s){
  REAL item1_5 = (3-s)*(3-s)*(3-s)*(3-s)*(3-s);  // (3-s)^5
  REAL item2_5 = (2-s)*(2-s)*(2-s)*(2-s)*(2-s);   // (2-s)^5
  REAL item3_5 = (1-s)*(1-s)*(1-s)*(1-s)*(1-s);   // (1-s)^5
  if(s<1){
      return factor_kernel*(item1_5-6*item2_5+15*item3_5);
  }
  else if(s<2){
      return factor_kernel*(item1_5-6*item2_5);
  }
  else if(s<3){
      return factor_kernel*item1_5;
  }
  else{
      return 0.0;
  }

    } // end kernelFunction

//  REAL factor = 1.0/(120.0*dem::PI*h*h*h*h);  // 3D quintic kernel factor
//  REAL factor = 1.0/(120.0*h*h);  // 1D quintic kernel factor
//  REAL factor = 7.0/(478.0*dem::PI*h*h*h);  // 2D quintic kernel factor
    // to calculate delta_aWab, where a is the position of the first particle
    inline dem::Vec SmoothParticleHydrodynamics::gradientKernelFunction(const dem::Vec& a, const dem::Vec& b){
  REAL rab = dem::vfabs(a-b);
  REAL s = rab*one_devide_h;
  REAL item1_4 = (3-s)*(3-s)*(3-s)*(3-s);  // (3-s)^4
  REAL item2_4 = (2-s)*(2-s)*(2-s)*(2-s); // (2-s)^4
  REAL item3_4 = (1-s)*(1-s)*(1-s)*(1-s); // (1-s)^4
  if(s<1){
      return factor_kernel_gradient *(-5*item1_4+30*item2_4-75*item3_4)*(a-b)/rab;  // it is strange that s is devided here. compare to Liu's SPH code (National University of Singapore) 
  }
  else if(s<2){
      return factor_kernel_gradient *(-5*item1_4+30*item2_4)*(a-b)/rab;
  }
  else if(s<3){
      return factor_kernel_gradient *(-5*item1_4)*(a-b)/rab;
  }
  else{
      return dem::Vec(0.0);
  }

    } // end gradientKernelFunction

    // to calculate partial differential dWab_dra
    inline REAL SmoothParticleHydrodynamics::partialKernelFunction(const dem::Vec& a, const dem::Vec& b){
  REAL rab = dem::vfabs(a-b);
  REAL s = rab*one_devide_h;
  REAL item1_4 = (3-s)*(3-s)*(3-s)*(3-s);  // (3-s)^4
  REAL item2_4 = (2-s)*(2-s)*(2-s)*(2-s); // (2-s)^4
  REAL item3_4 = (1-s)*(1-s)*(1-s)*(1-s); // (1-s)^4
  if(s<1){
      return factor_kernel_gradient *(-5*item1_4+30*item2_4-75*item3_4);
  }
  else if(s<2){
      return factor_kernel_gradient *(-5*item1_4+30*item2_4);
  }
  else if(s<3){
      return factor_kernel_gradient *(-5*item1_4);
  }
  else{
      return 0.0;
  }

    } // end partialKernelFunction


    // divide SPH domain in each cpu into different cells in 2D in xz plane.
    // see the notes 5/20/2015 and 5/21/2015 or Simpson's paper "Numerical techniques for three-dimensional Smoothed Particle Hydrodynamics"
    void SmoothParticleHydrodynamics::divideSPHDomain2D(){

  // clear the std::vector< std::vector<sph::SPHParticle*> > SPHParticleCellVec
  for(int pvec=0; pvec<SPHParticleCellVec.size(); pvec++){
      SPHParticleCellVec[pvec].clear();      
  }
  SPHParticleCellVec.clear();

      Vec  vmin = container.getMinCorner();
      Vec  vmax = container.getMaxCorner();
  REAL small_value = 0.01*spaceInterval;
      REAL xmin = vmin.getX()-sphCellSize-small_value; REAL zmin = vmin.getZ()-sphCellSize-small_value;  // expand the container by sphCellSize for cells domain is necessary
      REAL xmax = vmax.getX()+sphCellSize+small_value; REAL zmax = vmax.getZ()+sphCellSize+small_value; // since the sph domain that we divide is mergeSPHParticleVec
  Nx = (xmax-xmin)/kernelSize+1;  // (xmax-xmin)/(3h)+1
  Nz = (zmax-zmin)/kernelSize+1;  // (zmax-zmin)/(3h)+1
  numCell = Nx*Nz;
  SPHParticleCellVec.resize(numCell);  // at this point, the cellVec contains numCell vectors,
            // each vector is empty but ready to store the pointer of SPH particles (this process will be in the calculation of SPH forces)

  dem::Vec tmp_xyz;
  int num;
  for(std::vector<sph::SPHParticle*>::iterator pt=mergeSPHParticleVec.begin(); pt!=mergeSPHParticleVec.end(); pt++){ // mergeSPHParticle contains also the sph particles from neighboring cpus
      tmp_xyz = (*pt)->getCurrPosition();
      num = int( (tmp_xyz.getZ()-zmin)/kernelSize )*Nx + int( (tmp_xyz.getX()-xmin)/kernelSize );
      SPHParticleCellVec[num].push_back(*pt);
  }

  REAL maxRadius = gradation.getPtclMaxRadius();
  std::vector<sph::SPHParticle*>::iterator gt;
  for(std::vector<Particle*>::iterator pdt=particleVec.begin(); pdt!=particleVec.end(); ++pdt){  // the sph ghost particles in the current cpu's dem particles
                          // will be definitely inside the [xmin, ymin, zmin, xmax ...]
      tmp_xyz = (*pdt)->getCurrPos();
      if(tmp_xyz.getX()>=xmax+maxRadius || tmp_xyz.getZ()>=zmax+maxRadius
      || tmp_xyz.getX()<=xmin-maxRadius || tmp_xyz.getZ()<=zmin-maxRadius )  // dem particle is outside of the 
    continue;  
      for(gt=(*pdt)->SPHGhostParticleVec.begin(); gt!=(*pdt)->SPHGhostParticleVec.end(); gt++){  // set all SPH ghost particles into their cells
        tmp_xyz = (*gt)->getCurrPosition();
    if(tmp_xyz.getX()>=xmax || tmp_xyz.getZ()>=zmax
    || tmp_xyz.getX()<=xmin || tmp_xyz.getZ()<=zmin )
       continue;  // this sph ghost particle is outside of the container, go to next sph ghost particle
        num = int( (tmp_xyz.getZ()-zmin)/kernelSize )*Nx + int( (tmp_xyz.getX()-xmin)/kernelSize );
        SPHParticleCellVec[num].push_back(*gt);
      }
  }

  for(std::vector<Particle*>::iterator pdt=recvParticleVec.begin(); pdt!=recvParticleVec.end(); pdt++){  // the this time, recvParticleVec are the dem particles that are communicated by the commuParticle
      tmp_xyz = (*pdt)->getCurrPos();
      if(tmp_xyz.getX()>=xmax+maxRadius || tmp_xyz.getZ()>=zmax+maxRadius
      || tmp_xyz.getX()<=xmin-maxRadius || tmp_xyz.getZ()<=zmin-maxRadius )  // dem particle is outside of the 
    continue;                        // expanded domain 

      for(gt=(*pdt)->SPHGhostParticleVec.begin(); gt!=(*pdt)->SPHGhostParticleVec.end(); gt++){  // set part SPH ghost particles into their cells
        tmp_xyz = (*gt)->getCurrPosition();
    // not all sph ghost particles in these received particles are inside the container
    if(tmp_xyz.getX()>=xmax || tmp_xyz.getZ()>=zmax
    || tmp_xyz.getX()<=xmin || tmp_xyz.getZ()<=zmin )
       continue;  // this sph ghost particle is outside of the container, go to next sph ghost particle
        num = int( (tmp_xyz.getZ()-zmin)/kernelSize )*Nx + int( (tmp_xyz.getX()-xmin)/kernelSize );
        SPHParticleCellVec[num].push_back(*gt);
      }
  }

    } // divideSPHDomain2D

    //// the momentum equilibrium equation and state equation are implemented as in Monaghan's paper (1994), simulate free surface flow using sph
    //// here the neighboring list of SPH particles is searched by the cells, vector< vector<sph::SPHParticle*> > SPHParticleCellVec;
    void SmoothParticleHydrodynamics::calculateSPHDensityDotVelocityDotLinkedList2D(){

  // divide the SPH domain into different cells, each cell will contain SPH particles within it
  divideSPHDomain2D();

//  checkDivision();  // pass, May 22, 2015
//  checkNeighborCells2D();  // pass, May 22, 2015

  // initialize the densityDot and velocityDot of all the SPH particles 
  dem::Vec tmp_vec = dem::Vec(0,0, -(util::getParam<>("gravAccel"])*(util::getParam<>("gravScale"]) );
  dem::Vec zero_vec = dem::Vec(0,0,0);
  // in the calculation of forces between sph particles, we need to use the mergeSPHParticleVec, 
  // since mergeVec contains the sph particles from neighboring cpus. While we can only update and migrate and communicate SPHParticleVec
  for(std::vector<sph::SPHParticle*>::iterator pt=mergeSPHParticleVec.begin(); pt!=mergeSPHParticleVec.end(); pt++){
      (*pt)->setDensityDotVelocityDotZero();
      (*pt)->calculateParticleViscosity();  // everytime when use getParticleViscolisity(), make sure that pressure has been calculated!!!!
      (*pt)->calculateParticlePressure();    // everytime when use getParticlePressure(), make sure that pressure has been calculated!!!!
      switch((*pt)->getType()){
        case 1: // free sph particle
    (*pt)->addVelocityDot(tmp_vec);
    break;
        case 2: // ghost sph particle
    break;
        case 3: // boundary sph particle
    (*pt)->addVelocityDot(zero_vec);
    break;
        default:
    std::cout << "SPH particle type of pta should be 1, 2 or 3!" << std::endl;
    exit(-1);
      } // switch
  }

  std::vector<sph::SPHParticle*>::iterator gt;
  for(std::vector<Particle*>::iterator pdt=mergeParticleVec.begin(); pdt!=mergeParticleVec.end(); pdt++){  // all sph ghost particles
      for(gt=(*pdt)->SPHGhostParticleVec.begin(); gt!=(*pdt)->SPHGhostParticleVec.end(); gt++){
        (*gt)->setDensityDotVelocityDotZero();
        (*gt)->calculateParticleViscosity();  // everytime when use getParticleViscolisity(), make sure that pressure has been calculated!!!!
        (*gt)->calculateParticlePressure();  // everytime when use getParticlePressure(), make sure that pressure has been calculated!!!!
      }
  }

  // temporary variables used in the loop
  dem::Vec pta_position;
  dem::Vec ptb_position;
  dem::Vec delta_aWab;
  dem::Vec delta_bWba;
  dem::Vec vab;
  dem::Vec vba;
  dem::Vec vdem;
  dem::Vec dva_dt;
  dem::Vec dvb_dt;
  dem::Vec delta_a;
  dem::Vec delta_b;
  REAL pa, pb, rhoa, rhob, mua, mub;
  REAL rab;
        REAL dWab_dra;
  REAL dWba_drb;
  REAL Wab, Wba;
  REAL da, dB;
  REAL beta;
  REAL xa, ya, xB, yB, k, sqrt_xaya;  // variables for Morris' method to calculate da/dB
  std::vector<sph::SPHParticle*>::iterator ptb;
  REAL ra, rb, rc;  // the geometry of the dem particle
  dem::Vec pta_local, ptb_local;  // local position of pta and ptb in the dem particle
  dem::Vec demt_curr;
  REAL Gamma_ab, mu_ab, vr_dot;
  REAL alpha = util::getParam<>("alpha"];
    REAL alpha_zero = 0;  // the viscous between free and ghost/boundary
  REAL epsilon = util::getParam<>("epsilon"];  // parameter for velocity correction
  REAL Wq, Ra, Rb, phi_4, coefficient_a, coefficient_b;
  dem::Particle* demt;

        int pnum;
  std::vector<sph::SPHParticle*> tmp_particleVec;  // sph particles in neighboring cells
  for(int pvec=0; pvec<SPHParticleCellVec.size(); ++pvec){

      // store all SPH particles in pvec's neighboring cells
      tmp_particleVec.clear();  // it is the same for the particles in the same cell
      if(pvec+1<numCell){
    pnum = pvec+1;
            for(std::vector<sph::SPHParticle*>::iterator pt = SPHParticleCellVec[pnum].begin(); pt!= SPHParticleCellVec[pnum].end(); pt++) {
        tmp_particleVec.push_back(*pt);
            }
      }
      if(pvec+Nx-1<numCell){
    pnum = pvec+Nx-1;
            for(std::vector<sph::SPHParticle*>::iterator pt = SPHParticleCellVec[pnum].begin(); pt!= SPHParticleCellVec[pnum].end(); pt++) {
        tmp_particleVec.push_back(*pt);
            }
      }
      if(pvec+Nx<numCell){
    pnum = pvec+Nx;
            for(std::vector<sph::SPHParticle*>::iterator pt = SPHParticleCellVec[pnum].begin(); pt!= SPHParticleCellVec[pnum].end(); pt++) {
        tmp_particleVec.push_back(*pt);
            }
      }
      if(pvec+Nx+1<numCell){
    pnum = pvec+Nx+1;
            for(std::vector<sph::SPHParticle*>::iterator pt = SPHParticleCellVec[pnum].begin(); pt!= SPHParticleCellVec[pnum].end(); pt++) {
        tmp_particleVec.push_back(*pt);
            }
      }

      for(std::vector<sph::SPHParticle*>::iterator pta=SPHParticleCellVec[pvec].begin(); pta!=SPHParticleCellVec[pvec].end(); ++pta){  // SPH particles in cell pvec
    pa = (*pta)->getParticlePressure();
        if(pa>=0){
        Ra = 0.006;
        }
        else{
        Ra = 0.6;
        }
//        mua = (*pta)->getParticleViscosity();
        rhoa = (*pta)->getParticleDensity();
        pta_position = (*pta)->getCurrPosition();

//    if((*pta)->getType()!=1){  // pta is not free sph particles, i.e. pta is ghost or boundary particles
//        continue;  // do not consider pta, as we treat before
//    }    

    for(ptb=pta+1; ptb!=SPHParticleCellVec[pvec].end(); ++ptb){  // sum over the SPH particles in the same cell pvec
        ptb_position = (*ptb)->getCurrPosition();
        rab = dem::vfabs(pta_position-ptb_position);
        if(rab<=kernelSize){  // ptb is in the smooth kernel
          pb = (*ptb)->getParticlePressure();
//          mub = (*ptb)->getParticleViscosity();
          rhob = (*ptb)->getParticleDensity();

          Wq = kernelFunction(rab/smoothLength);
          phi_4 = pow(Wq/Wqmin,4);
          if(pb>=0){
          Rb = 0.006;
          }
          else{
          Rb = 0.6;
          }
          coefficient_a = 1+Ra*phi_4;
          coefficient_b = 1+Rb*phi_4;

      // we have three types of SPH particles: 1, free particle; 2, ghost particle; 3, boundary particle
      // Then we have 3x3 = 9 different types of interactions with the three types of particles
      // so we cannot judge the interaction type only by the type of particle ptb, we need to consider pta also
      // pta          ptb          need to consider or not
      // 1            1              V      free with free
      // 1            2              V      free with ghost
      // 1            3              V      free with boundary

      // 2            1              V      ghost with free
      // 2            2              X      ghost with ghost
      // 2            3              X      ghost with boundary

      // 3            1              V      boundary with free
      // 3            2              X      boundary with ghost
      // 3            3              X       boundary with boundary


      // add the density dot for pta and ptb
          delta_aWab = gradientKernelFunction(pta_position, ptb_position);  // this is to add SPH pta 
          delta_bWba = - delta_aWab;
          dWab_dra = partialKernelFunction(pta_position, ptb_position);  // this is to add SPH pta
          dWba_drb = dWab_dra;
          Wab = kernelFunction(pta_position, ptb_position);  // this is to add SPH pta
          Wba = Wab; 

      switch((*ptb)->getType()){
        case 1:  // ptb is free SPH particle
          switch((*pta)->getType()){
            case 1:  // free with free
                vab = (*pta)->getVelocity()-(*ptb)->getVelocity();
               vba = -vab;

                (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab*delta_aWab) );
                (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba*delta_bWba) );

                vr_dot = vab*(pta_position-ptb_position);
                if(vr_dot<0){
            mu_ab = smoothLength*vr_dot/(rab*rab+0.01*smoothLength*smoothLength);
                Gamma_ab = (-alpha*(util::getParam<>("soundSpeed"])*mu_ab)/(rhoa+rhob)*2;
                }
                else{
             Gamma_ab = 0;
                }


                dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
                (*pta)->addVelocityDot(dva_dt);
                dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)*coefficient_b+pa/(rhoa*rhoa)*coefficient_a+Gamma_ab)*delta_bWba;
                (*ptb)->addVelocityDot(dvb_dt);

                delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
                (*pta)->addVelocityCorrection(delta_a);
                delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
                (*ptb)->addVelocityCorrection(delta_b);
              break;
            case 2:  // ghost with free
        demt = (*pta)->getDemParticle();
        demt_curr = demt->getCurrPos();
            ptb_local = demt->globalToLocal(ptb_position-demt_curr);  // the local position of sph point ptb
            ra = demt->getA(); rc = demt->getC();
            k = 1.0/(sqrt(ptb_local.getX()*ptb_local.getX()/(ra*ra)+ptb_local.getZ()*ptb_local.getZ()/(rc*rc) ) );
            da = dem::vfabs(ptb_local-k*ptb_local);  // the distance is the same in rotated coordinates

                // calculate Vab as the method shown in Morris's paper, 1996

                // (1) here I wanna use the distance from the point a/b to the surface of the ellipsoid to simplify the problem
            pta_local = (*pta)->getLocalPosition();
            k = 1.0/(sqrt(pta_local.getX()*pta_local.getX()/(ra*ra)+pta_local.getZ()*pta_local.getZ()/(rc*rc) ) );
                dB = dem::vfabs(pta_local-k*pta_local);

//              // (2) here the Morris's method is used 
//              xa = pta_position.getx(); ya = pta_position.gety();
//              xB = ptb_position.getx(); yB = ptb_position.gety();
//              if(ya==0) {da = xa-radius; dB = radius-xB;}
//              else if(xa==0) {da = ya-radius; dB = radius-yB;}
//              else {
//          sqrt_xaya = sqrt(xa*xa+ya*ya);
//          k = radius/sqrt_xaya;
//          da = radius/k - radius;
//          dB = fabs(xa*xB+ya*yB-k*radius*radius)/sqrt_xaya;
//              }


                beta = 1+dB/da;
                if(beta>2.0 || beta<0 || isnan(beta)){ beta = 2.0; }

            vdem = (*pta)->getVelocity();
                vba = beta*((*ptb)->getVelocity()-vdem);
                vab = -vba;

                (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab*delta_aWab) );
                (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba*delta_bWba) );

                vr_dot = vab*(pta_position-ptb_position);
                if(vr_dot<0){
                mu_ab = smoothLength*vr_dot/(rab*rab+0.01*smoothLength*smoothLength);
                    Gamma_ab = (-alpha*(util::getParam<>("soundSpeed"])*mu_ab)/(rhoa+rhob)*2;
                }
                else{
                 Gamma_ab = 0;
                }

                dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
            demt->addForce((*pta)->getParticleMass()*dva_dt);
            demt->addMoment( (pta_position-demt_curr) % ((*pta)->getParticleMass()*dva_dt) );
//                (*pta)->addVelocityDot(dva_dt);

                dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)*coefficient_b+pa/(rhoa*rhoa)*coefficient_a+Gamma_ab)*delta_bWba;
                (*ptb)->addVelocityDot(dvb_dt);  // the velocities of the ghost particles will not envolve as the same way as others

//                  delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
//                  (*pta)->addVelocityCorrection(delta_a);
                  delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
                  (*ptb)->addVelocityCorrection(delta_b);

        break;

            case 3:  // boundary with free

        // calculate Vab as the method shown in Morris's paper, 1996
                // interact with boundary particles
                da = ptb_position.getZ()-allContainer.getMinCorner().getZ();  // assume with the bottom boundary
        dB = allContainer.getMinCorner().getZ()-pta_position.getX();  // assume with the bottom boundary
        if(pta_position.getX()<allContainer.getMinCorner().getX()){  // with left boundary
            da = ptb_position.getX()-allContainer.getMinCorner().getX();
            dB = allContainer.getMinCorner().getX()-pta_position.getX();
        }
        else if(pta_position.getX()>allContainer.getMaxCorner().getX()){ // with right boundary
            da = allContainer.getMaxCorner().getX()-ptb_position.getX();
            dB = pta_position.getX()-allContainer.getMaxCorner().getX();
        }
        else if(pta_position.getZ()>allContainer.getMaxCorner().getZ()){ // with top boundary
            da = allContainer.getMaxCorner().getZ()-ptb_position.getZ();
            dB = pta_position.getZ()-allContainer.getMaxCorner().getZ();
        }

                beta = 1+dB/da;
                if(beta>2 || beta<0 || isnan(beta)){ beta = 2; }

                vba = beta*(*ptb)->getVelocity();
                vab = -vba;

                (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab*delta_aWab) );
                (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba*delta_bWba) );

                vr_dot = vab*(pta_position-ptb_position);
                if(vr_dot<0){
            mu_ab = smoothLength*vr_dot/(rab*rab+0.01*smoothLength*smoothLength);
                Gamma_ab = (-alpha*(util::getParam<>("soundSpeed"])*mu_ab)/(rhoa+rhob)*2;
                }
                else{
             Gamma_ab = 0;
                }

//                dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
//                (*pta)->addVelocityDot(dva_dt);
                dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)+pa/(rhoa*rhoa)+Gamma_ab)*delta_bWba;
                (*ptb)->addVelocityDot(dvb_dt);  // the velocities of the ghost particles will not envolve as the same way as others

//            delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
//                (*pta)->addVelocityCorrection(delta_a);
                delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
                (*ptb)->addVelocityCorrection(delta_b);
      
//            // apply the boundary forces by Lennard-Jones potential as in Monaghan's paper(1994)
//                if(rab<=spaceInterval){ // ptb is in the smooth kernel
//                dvb_dt = D*(pow(spaceInterval/rab, p1)-pow(spaceInterval/rab, p2))*(ptb_position-pta_position)/(rab*rab);
//                (*ptb)->addVelocityDot(dvb_dt);
//                } // end if

        break;
            default:
            std::cout << "SPH particle type of pta should be 1, 2 or 3!" << std::endl;
            exit(-1);
            } // end switch pta

          break;
        case 2:  // ptb is ghost particle

          if((*pta)->getType()!=1){  // pta is not free sph particles, i.e. pta is ghost or boundary particles
            break;  // do not consider pta, as we treat before
          }
          demt = (*ptb)->getDemParticle();
          demt_curr = demt->getCurrPos();
          pta_local = demt->globalToLocal(pta_position-demt_curr);  // the local position of sph point pta
          ra = demt->getA(); rc = demt->getC();
          k = 1.0/(sqrt(pta_local.getX()*pta_local.getX()/(ra*ra)+pta_local.getZ()*pta_local.getZ()/(rc*rc) ) );
          da = dem::vfabs(pta_local-k*pta_local);  // the distance is the same in rotated coordinates

              // calculate Vab as the method shown in Morris's paper, 1996

              // (1) here I wanna use the distance from the point a/b to the surface of the ellipsoid to simplify the problem
          ptb_local = (*ptb)->getLocalPosition();
          k = 1.0/(sqrt(ptb_local.getX()*ptb_local.getX()/(ra*ra)+ptb_local.getZ()*ptb_local.getZ()/(rc*rc) ) );
              dB = dem::vfabs(ptb_local-k*ptb_local);

//              // (2) here the Morris's method is used 
//              xa = pta_position.getx(); ya = pta_position.gety();
//              xB = ptb_position.getx(); yB = ptb_position.gety();
//              if(ya==0) {da = xa-radius; dB = radius-xB;}
//              else if(xa==0) {da = ya-radius; dB = radius-yB;}
//              else {
//          sqrt_xaya = sqrt(xa*xa+ya*ya);
//          k = radius/sqrt_xaya;
//          da = radius/k - radius;
//          dB = fabs(xa*xB+ya*yB-k*radius*radius)/sqrt_xaya;
//              }


              beta = 1+dB/da;
              if(beta>2 || beta<0 || isnan(beta)){ beta = 2; }

          vdem = (*ptb)->getVelocity();
              vab = beta*((*pta)->getVelocity()-vdem);
              vba = -vab;

              (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab*delta_aWab) );
              (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba*delta_bWba) );

              vr_dot = vab*(pta_position-ptb_position);
              if(vr_dot<0){
            mu_ab = smoothLength*vr_dot/(rab*rab+0.01*smoothLength*smoothLength);
                Gamma_ab = (-alpha*(util::getParam<>("soundSpeed"])*mu_ab)/(rhoa+rhob)*2;
              }
              else{
             Gamma_ab = 0;
              }

              dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
              (*pta)->addVelocityDot(dva_dt);
              dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)*coefficient_b+pa/(rhoa*rhoa)*coefficient_a+Gamma_ab)*delta_bWba;
          demt->addForce((*ptb)->getParticleMass()*dvb_dt);
          demt->addMoment( (ptb_position-demt_curr) % ((*ptb)->getParticleMass()*dvb_dt) );
//              (*ptb)->addVelocityDot(dvb_dt);  // the velocities of the ghost particles will not envolve as the same way as others

                delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
                (*pta)->addVelocityCorrection(delta_a);
//                delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
//                (*ptb)->addVelocityCorrection(delta_b);

          break;
        case 3:  // ptb is boundary particle

          if((*pta)->getType()!=1){  // pta is not free sph particles, i.e. pta is ghost or boundary particles
            break;  // do not consider pta, as we treat before
          }

          // calculate Vab as the method shown in Morris's paper, 1996
              // interact with boundary particles
              da = pta_position.getZ()-allContainer.getMinCorner().getZ();  // assume with the bottom boundary
          dB = allContainer.getMinCorner().getZ()-ptb_position.getZ();  // assume with the bottom boundary
          if(ptb_position.getX() < allContainer.getMinCorner().getX()){  // with left boundary
        da = pta_position.getX() - allContainer.getMinCorner().getX();
        dB = allContainer.getMinCorner().getX() - ptb_position.getX();
          }
          else if(ptb_position.getX() > allContainer.getMaxCorner().getX()){  // with right boundary
        da = allContainer.getMaxCorner().getX() - pta_position.getX();
        dB = ptb_position.getX() - allContainer.getMaxCorner().getX();
          }
          else if(ptb_position.getZ() > allContainer.getMaxCorner().getZ()){  // with top boundary
        da = allContainer.getMaxCorner().getZ()-pta_position.getZ();
        dB = ptb_position.getZ()-allContainer.getMaxCorner().getZ();
          }

              beta = 1+dB/da;
              if(beta>2 || beta<0 || isnan(beta)){ beta = 2; }

              vab = beta*(*pta)->getVelocity();
              vba = -vab;

              (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab*delta_aWab) );
              (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba*delta_bWba) );

              vr_dot = vab*(pta_position-ptb_position);
              if(vr_dot<0){
        mu_ab = smoothLength*vr_dot/(rab*rab+0.01*smoothLength*smoothLength);
            Gamma_ab = (-alpha*(util::getParam<>("soundSpeed"])*mu_ab)/(rhoa+rhob)*2;
              }
              else{
         Gamma_ab = 0;
              }

              dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
              (*pta)->addVelocityDot(dva_dt);
//              dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)+pa/(rhoa*rhoa)+Gamma_ab)*delta_bWba;
//              (*ptb)->addVelocityDot(dvb_dt);  // the velocities of the ghost particles will not envolve as the same way as others

              delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
              (*pta)->addVelocityCorrection(delta_a);
//              delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
//              (*ptb)->addVelocityCorrection(delta_b);
      
//          // apply the boundary forces by Lennard-Jones potential as in Monaghan's paper(1994)
//              if(rab<=spaceInterval){ // ptb is in the smooth kernel
//            dva_dt = D*(pow(spaceInterval/rab, p1)-pow(spaceInterval/rab, p2))*(pta_position-ptb_position)/(rab*rab);
//            (*pta)->addVelocityDot(dva_dt);
//              } // end if
          break;
        default:
          std::cout << "SPH particle type should be 1, 2 or 3!" << std::endl;
          exit(-1);
    
      } //end swtich type
        } // end if 3h
    } // end for ptb in the same cell

    for(ptb=tmp_particleVec.begin(); ptb!=tmp_particleVec.end(); ptb++){  // all particles in pvec's neighboring cells 
        ptb_position = (*ptb)->getCurrPosition();
        rab = dem::vfabs(pta_position-ptb_position);
        if(rab<=kernelSize){  // ptb is in the smooth kernel
          pb = (*ptb)->getParticlePressure();
//          mub = (*ptb)->getParticleViscosity();
          rhob = (*ptb)->getParticleDensity();

          Wq = kernelFunction(rab/smoothLength);
          phi_4 = pow(Wq/Wqmin,4);
          if(pb>=0){
          Rb = 0.006;
          }
          else{
          Rb = 0.6;
          }
          coefficient_a = 1+Ra*phi_4;
          coefficient_b = 1+Rb*phi_4;

      // add the density dot for pta and ptb
          delta_aWab = gradientKernelFunction(pta_position, ptb_position);  // this is to add SPH pta 
          delta_bWba = - delta_aWab;
          dWab_dra = partialKernelFunction(pta_position, ptb_position);  // this is to add SPH pta
          dWba_drb = dWab_dra;
          Wab = kernelFunction(pta_position, ptb_position);  // this is to add SPH pta
          Wba = Wab; 
      switch((*ptb)->getType()){
        case 1:  // ptb is free SPH particle
          switch((*pta)->getType()){
            case 1:  // free with free
                vab = (*pta)->getVelocity()-(*ptb)->getVelocity();
               vba = -vab;

                (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab*delta_aWab) );
                (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba*delta_bWba) );

                vr_dot = vab*(pta_position-ptb_position);
                if(vr_dot<0){
            mu_ab = smoothLength*vr_dot/(rab*rab+0.01*smoothLength*smoothLength);
                Gamma_ab = (-alpha*(util::getParam<>("soundSpeed"])*mu_ab)/(rhoa+rhob)*2;
                }
                else{
             Gamma_ab = 0;
                }


                dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
                (*pta)->addVelocityDot(dva_dt);
                dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)*coefficient_b+pa/(rhoa*rhoa)*coefficient_a+Gamma_ab)*delta_bWba;
                (*ptb)->addVelocityDot(dvb_dt);

                delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
                (*pta)->addVelocityCorrection(delta_a);
                delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
                (*ptb)->addVelocityCorrection(delta_b);
              break;
            case 2:  // ghost with free
        demt = (*pta)->getDemParticle();
        demt_curr = demt->getCurrPos();
            ptb_local = demt->globalToLocal(ptb_position-demt_curr);  // the local position of sph point ptb
            ra = demt->getA(); rc = demt->getC();
            k = 1.0/(sqrt(ptb_local.getX()*ptb_local.getX()/(ra*ra)+ptb_local.getZ()*ptb_local.getZ()/(rc*rc) ) );
            da = dem::vfabs(ptb_local-k*ptb_local);  // the distance is the same in rotated coordinates

                // calculate Vab as the method shown in Morris's paper, 1996

                // (1) here I wanna use the distance from the point a/b to the surface of the ellipsoid to simplify the problem
            pta_local = (*pta)->getLocalPosition();
            k = 1.0/(sqrt(pta_local.getX()*pta_local.getX()/(ra*ra)+pta_local.getZ()*pta_local.getZ()/(rc*rc) ) );
                dB = dem::vfabs(pta_local-k*pta_local);

//              // (2) here the Morris's method is used 
//              xa = pta_position.getx(); ya = pta_position.gety();
//              xB = ptb_position.getx(); yB = ptb_position.gety();
//              if(ya==0) {da = xa-radius; dB = radius-xB;}
//              else if(xa==0) {da = ya-radius; dB = radius-yB;}
//              else {
//          sqrt_xaya = sqrt(xa*xa+ya*ya);
//          k = radius/sqrt_xaya;
//          da = radius/k - radius;
//          dB = fabs(xa*xB+ya*yB-k*radius*radius)/sqrt_xaya;
//              }


                beta = 1+dB/da;
                if(beta>2 || beta<0 || isnan(beta)){ beta = 2; }

            vdem = (*pta)->getVelocity();
                vba = beta*((*ptb)->getVelocity()-vdem);
                vab = -vba;

                (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab*delta_aWab) );
                (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba*delta_bWba) );

                vr_dot = vab*(pta_position-ptb_position);
                if(vr_dot<0){
                mu_ab = smoothLength*vr_dot/(rab*rab+0.01*smoothLength*smoothLength);
                    Gamma_ab = (-alpha*(util::getParam<>("soundSpeed"])*mu_ab)/(rhoa+rhob)*2;
                }
                else{
                 Gamma_ab = 0;
                }

                dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
            demt->addForce((*pta)->getParticleMass()*dva_dt);
            demt->addMoment( (pta_position-demt_curr) % ((*pta)->getParticleMass()*dva_dt) );
//                (*pta)->addVelocityDot(dva_dt);

                dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)*coefficient_b+pa/(rhoa*rhoa)*coefficient_a+Gamma_ab)*delta_bWba;
                (*ptb)->addVelocityDot(dvb_dt);  // the velocities of the ghost particles will not envolve as the same way as others

//                  delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
//                  (*pta)->addVelocityCorrection(delta_a);
                  delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
                  (*ptb)->addVelocityCorrection(delta_b);

        break;

            case 3:  // boundary with free

        // calculate Vab as the method shown in Morris's paper, 1996
                // interact with boundary particles
                da = ptb_position.getZ()-allContainer.getMinCorner().getZ();  // assume with the bottom boundary
        dB = allContainer.getMinCorner().getZ()-pta_position.getX();  // assume with the bottom boundary
        if(pta_position.getX()<allContainer.getMinCorner().getX()){  // with left boundary
            da = ptb_position.getX()-allContainer.getMinCorner().getX();
            dB = allContainer.getMinCorner().getX()-pta_position.getX();
        }
        else if(pta_position.getX()>allContainer.getMaxCorner().getX()){ // with right boundary
            da = allContainer.getMaxCorner().getX()-ptb_position.getX();
            dB = pta_position.getX()-allContainer.getMaxCorner().getX();
        }
        else if(pta_position.getZ()>allContainer.getMaxCorner().getZ()){ // with top boundary
            da = allContainer.getMaxCorner().getZ()-ptb_position.getZ();
            dB = pta_position.getZ()-allContainer.getMaxCorner().getZ();
        }

                beta = 1+dB/da;
                if(beta>2 || beta<0 || isnan(beta)){ beta = 2; }

                vba = beta*(*ptb)->getVelocity();
                vab = -vba;

                (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab*delta_aWab) );
                (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba*delta_bWba) );

                vr_dot = vab*(pta_position-ptb_position);
                if(vr_dot<0){
            mu_ab = smoothLength*vr_dot/(rab*rab+0.01*smoothLength*smoothLength);
                Gamma_ab = (-alpha*(util::getParam<>("soundSpeed"])*mu_ab)/(rhoa+rhob)*2;
                }
                else{
             Gamma_ab = 0;
                }

//                dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
//                (*pta)->addVelocityDot(dva_dt);
                dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)+pa/(rhoa*rhoa)+Gamma_ab)*delta_bWba;
                (*ptb)->addVelocityDot(dvb_dt);  // the velocities of the ghost particles will not envolve as the same way as others

//            delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
//                (*pta)->addVelocityCorrection(delta_a);
                delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
                (*ptb)->addVelocityCorrection(delta_b);
      
            // apply the boundary forces by Lennard-Jones potential as in Monaghan's paper(1994)
                if(rab<=spaceInterval){ // ptb is in the smooth kernel
                dvb_dt = D*(pow(spaceInterval/rab, p1)-pow(spaceInterval/rab, p2))*(ptb_position-pta_position)/(rab*rab);
                (*ptb)->addVelocityDot(dvb_dt);
                } // end if

        break;
            default:
            std::cout << "SPH particle type of pta should be 1, 2 or 3!" << std::endl;
            exit(-1);
            } // end switch pta

          break;
        case 2:  // ptb is ghost particle


          if((*pta)->getType()!=1){  // pta is not free sph particles, i.e. pta is ghost or boundary particles
            break;  // do not consider pta, as we treat before
          }
          demt = (*ptb)->getDemParticle();
          demt_curr = demt->getCurrPos();
          pta_local = demt->globalToLocal(pta_position-demt_curr);  // the local position of sph point pta
          ra = demt->getA(); rc = demt->getC();
          k = 1.0/(sqrt(pta_local.getX()*pta_local.getX()/(ra*ra)+pta_local.getZ()*pta_local.getZ()/(rc*rc) ) );
          da = dem::vfabs(pta_local-k*pta_local);  // the distance is the same in rotated coordinates

              // calculate Vab as the method shown in Morris's paper, 1996

              // (1) here I wanna use the distance from the point a/b to the surface of the ellipsoid to simplify the problem
          ptb_local = (*ptb)->getLocalPosition();
          k = 1.0/(sqrt(ptb_local.getX()*ptb_local.getX()/(ra*ra)+ptb_local.getZ()*ptb_local.getZ()/(rc*rc) ) );
              dB = dem::vfabs(ptb_local-k*ptb_local);

//              // (2) here the Morris's method is used 
//              xa = pta_position.getx(); ya = pta_position.gety();
//              xB = ptb_position.getx(); yB = ptb_position.gety();
//              if(ya==0) {da = xa-radius; dB = radius-xB;}
//              else if(xa==0) {da = ya-radius; dB = radius-yB;}
//              else {
//          sqrt_xaya = sqrt(xa*xa+ya*ya);
//          k = radius/sqrt_xaya;
//          da = radius/k - radius;
//          dB = fabs(xa*xB+ya*yB-k*radius*radius)/sqrt_xaya;
//              }


              beta = 1+dB/da;
              if(beta>2 || beta<0 || isnan(beta)){ beta = 2; }

          vdem = (*ptb)->getVelocity();
              vab = beta*((*pta)->getVelocity()-vdem);
              vba = -vab;

              (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab*delta_aWab) );
              (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba*delta_bWba) );

              vr_dot = vab*(pta_position-ptb_position);
              if(vr_dot<0){
            mu_ab = smoothLength*vr_dot/(rab*rab+0.01*smoothLength*smoothLength);
                Gamma_ab = (-alpha*(util::getParam<>("soundSpeed"])*mu_ab)/(rhoa+rhob)*2;
              }
              else{
             Gamma_ab = 0;
              }

              dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
              (*pta)->addVelocityDot(dva_dt);
              dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)*coefficient_b+pa/(rhoa*rhoa)*coefficient_a+Gamma_ab)*delta_bWba;
          demt->addForce((*ptb)->getParticleMass()*dvb_dt);
          demt->addMoment( (ptb_position-demt_curr) % ((*ptb)->getParticleMass()*dvb_dt) );
//              (*ptb)->addVelocityDot(dvb_dt);  // the velocities of the ghost particles will not envolve as the same way as others

                delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
                (*pta)->addVelocityCorrection(delta_a);
//                delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
//                (*ptb)->addVelocityCorrection(delta_b);

          break;

        case 3:  // ptb is boundary particle

          if((*pta)->getType()!=1){  // pta is not free sph particles, i.e. pta is ghost or boundary particles
            break;  // do not consider pta, as we treat before
          }

          // calculate Vab as the method shown in Morris's paper, 1996
              // interact with boundary particles
              da = pta_position.getZ()-allContainer.getMinCorner().getZ();  // assume with the bottom boundary
          dB = allContainer.getMinCorner().getZ()-ptb_position.getZ();  // assume with the bottom boundary
          if(ptb_position.getX() < allContainer.getMinCorner().getX()){  // with left boundary
        da = pta_position.getX() - allContainer.getMinCorner().getX();
        dB = allContainer.getMinCorner().getX() - ptb_position.getX();
          }
          else if(ptb_position.getX() > allContainer.getMaxCorner().getX()){  // with right boundary
        da = allContainer.getMaxCorner().getX() - pta_position.getX();
        dB = ptb_position.getX() - allContainer.getMaxCorner().getX();
          }
          else if(ptb_position.getZ() > allContainer.getMaxCorner().getZ()){  // with top boundary
        da = allContainer.getMaxCorner().getZ()-pta_position.getZ();
        dB = ptb_position.getZ()-allContainer.getMaxCorner().getZ();
          }

              beta = 1+dB/da;
              if(beta>2 || beta<0 || isnan(beta)){ beta = 2; }

              vab = beta*(*pta)->getVelocity();
              vba = -vab;

              (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab*delta_aWab) );
              (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba*delta_bWba) );

              vr_dot = vab*(pta_position-ptb_position);
              if(vr_dot<0){
        mu_ab = smoothLength*vr_dot/(rab*rab+0.01*smoothLength*smoothLength);
            Gamma_ab = (-alpha*(util::getParam<>("soundSpeed"])*mu_ab)/(rhoa+rhob)*2;
              }
              else{
         Gamma_ab = 0;
              }

              dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
              (*pta)->addVelocityDot(dva_dt);
//              dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)+pa/(rhoa*rhoa)+Gamma_ab)*delta_bWba;
//              (*ptb)->addVelocityDot(dvb_dt);  // the velocities of the ghost particles will not envolve as the same way as others

              delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
              (*pta)->addVelocityCorrection(delta_a);
//              delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
//              (*ptb)->addVelocityCorrection(delta_b);

              // apply the boundary forces by Lennard-Jones potential as in Monaghan's paper(1994)
              if(rab<=spaceInterval){ // ptb is in the smooth kernel
            dva_dt = D*(pow(spaceInterval/rab, p1)-pow(spaceInterval/rab, p2))*(pta_position-ptb_position)/(rab*rab);
            (*pta)->addVelocityDot(dva_dt);
              } // end if
      
          break;
        default:
          std::cout << "SPH particle type should be 1, 2 or 3!" << std::endl;
          exit(-1);
    
      } //end swtich type
        } // end if 3h
    } // end for ptb in neighbor cells
      } // end for pta

      tmp_particleVec.clear();  // clear elements in tmp-vector for particles neighboring cells, it is important

  } // end for pvec, different cells


//      // apply the boundary forces by Lennard-Jones potential as in Monaghan's paper(1994)
//      for(ptb=SPHBoundaryParticleVec.begin(); ptb!=SPHBoundaryParticleVec.end(); ptb++){  
//    ptb_position = (*ptb)->getCurrPosition();
//    rab = dem::vfabs(pta_position-ptb_position);
//    if(rab<=spaceInterval){ // ptb is in the smooth kernel
//        dva_dt = D*(pow(spaceInterval/rab, p1)-pow(spaceInterval/rab, p2))*(pta_position-ptb_position)/(rab*rab);
//        (*pta)->addVelocityDot(dva_dt);
//    } // end if
//      } // end ptb


    } // end calculateSPHDensityDotVelocityDotLinkedList2D()


    // divide SPH domain in each cpu into different cells in 2D in xz plane.
    // see the notes 5/20/2015 and 5/21/2015 or Simpson's paper "Numerical techniques for three-dimensional Smoothed Particle Hydrodynamics"
    void SmoothParticleHydrodynamics::divideSPHDomain3D(){

  // clear the std::vector< std::vector<sph::SPHParticle*> > SPHParticleCellVec
  for(int pvec=0; pvec<SPHParticleCellVec.size(); pvec++){
      SPHParticleCellVec[pvec].clear();      
  }
  SPHParticleCellVec.clear();

      Vec  vmin = container.getMinCorner();
      Vec  vmax = container.getMaxCorner();
  REAL small_value = 0.01*spaceInterval;
      REAL xmin = vmin.getX()-sphCellSize-small_value; REAL ymin = vmin.getY()-sphCellSize-small_value; REAL zmin = vmin.getZ()-sphCellSize-small_value; // expand the container by sphCellSize for cells domain is necessary
      REAL xmax = vmax.getX()+sphCellSize+small_value; REAL ymax = vmax.getY()+sphCellSize+small_value; REAL zmax = vmax.getZ()+sphCellSize+small_value; // since the sph domain that we divide is mergeSPHParticleVec
  Nx = (xmax-xmin)/kernelSize+1;  // (xmax-xmin)/(3h)+1
  Nz = (zmax-zmin)/kernelSize+1;  // (zmax-zmin)/(3h)+1
  Ny = (ymax-ymin)/kernelSize+1;  // (ymax-ymin)/(3h)+1
  numCell = Nx*Nz*Ny;
  SPHParticleCellVec.resize(numCell);  // at this point, the cellVec contains numCell vectors,
            // each vector is empty but ready to store the pointer of SPH particles (this process will be in the calculation of SPH forces)
        pnum_vec.clear();
  pnum_vec.push_back(1); pnum_vec.push_back(Nx-1); pnum_vec.push_back(Nx);
  pnum_vec.push_back(Nx+1); pnum_vec.push_back(Nx*Nz); pnum_vec.push_back(Nx*Nz-Nx);
  pnum_vec.push_back(Nx*Nz+Nx); pnum_vec.push_back(Nx*Nz-1); pnum_vec.push_back(Nx*Nz-1-Nx);
  pnum_vec.push_back(Nx*Nz-1+Nx); pnum_vec.push_back(Nx*Nz+1); pnum_vec.push_back(Nx*Nz+1-Nx);
  pnum_vec.push_back(Nx*Nz+1+Nx);

  dem::Vec tmp_xyz;
  int num;
  for(std::vector<sph::SPHParticle*>::iterator pt=mergeSPHParticleVec.begin(); pt!=mergeSPHParticleVec.end(); pt++){ // mergeSPHParticle contains also the sph particles from neighboring cpus
      tmp_xyz = (*pt)->getCurrPosition();
      num = int( (tmp_xyz.getY()-ymin)/kernelSize )*Nx*Nz + int( (tmp_xyz.getZ()-zmin)/kernelSize )*Nx + int( (tmp_xyz.getX()-xmin)/kernelSize );
      SPHParticleCellVec[num].push_back(*pt);
  }

  REAL maxRadius = gradation.getPtclMaxRadius();
  std::vector<sph::SPHParticle*>::iterator gt;
  for(std::vector<Particle*>::iterator pdt=particleVec.begin(); pdt!=particleVec.end(); pdt++){  // the sph ghost particles in the current cpu's dem particles
                          // will be definitely inside the [xmin, ymin, zmin, xmax ...]
      tmp_xyz = (*pdt)->getCurrPos();
      if(tmp_xyz.getX()>=xmax+maxRadius || tmp_xyz.getY()>=ymax+maxRadius || tmp_xyz.getZ()>=zmax+maxRadius
      || tmp_xyz.getX()<=xmin-maxRadius || tmp_xyz.getY()<=ymin-maxRadius || tmp_xyz.getZ()<=zmin-maxRadius )  // dem particle is outside of the 
    continue;  
      for(gt=(*pdt)->SPHGhostParticleVec.begin(); gt!=(*pdt)->SPHGhostParticleVec.end(); gt++){  // set all SPH ghost particles into their cells
        tmp_xyz = (*gt)->getCurrPosition();
    if(tmp_xyz.getX()>=xmax || tmp_xyz.getY()>=ymax || tmp_xyz.getZ()>=zmax
    || tmp_xyz.getX()<=xmin || tmp_xyz.getY()<=ymin || tmp_xyz.getZ()<=zmin )
       continue;  // this sph ghost particle is outside of the container, go to next sph ghost particle
        num = int( (tmp_xyz.getY()-ymin)/kernelSize )*Nx*Nz + int( (tmp_xyz.getZ()-zmin)/kernelSize )*Nx + int( (tmp_xyz.getX()-xmin)/kernelSize );
        SPHParticleCellVec[num].push_back(*gt);
      }
  }

  for(std::vector<Particle*>::iterator pdt=recvParticleVec.begin(); pdt!=recvParticleVec.end(); pdt++){  // the this time, recvParticleVec are the dem particles that are communicated by the commuParticle
      tmp_xyz = (*pdt)->getCurrPos();
      if(tmp_xyz.getX()>=xmax+maxRadius || tmp_xyz.getY()>=ymax+maxRadius || tmp_xyz.getZ()>=zmax+maxRadius
      || tmp_xyz.getX()<=xmin-maxRadius || tmp_xyz.getY()<=ymin-maxRadius || tmp_xyz.getZ()<=zmin-maxRadius )  // dem particle is outside of the 
    continue;                        // expanded domain 

      for(gt=(*pdt)->SPHGhostParticleVec.begin(); gt!=(*pdt)->SPHGhostParticleVec.end(); gt++){  // set part SPH ghost particles into their cells
        tmp_xyz = (*gt)->getCurrPosition();
    // not all sph ghost particles in these received particles are inside the container
    if(tmp_xyz.getX()>=xmax || tmp_xyz.getY()>=ymax || tmp_xyz.getZ()>=zmax
    || tmp_xyz.getX()<=xmin || tmp_xyz.getY()<=ymin || tmp_xyz.getZ()<=zmin )
       continue;  // this sph ghost particle is outside of the container, go to next sph ghost particle
        num = int( (tmp_xyz.getY()-ymin)/kernelSize )*Nx*Nz + int( (tmp_xyz.getZ()-zmin)/kernelSize )*Nx + int( (tmp_xyz.getX()-xmin)/kernelSize );
        SPHParticleCellVec[num].push_back(*gt);
      }
  }

    } // divideSPHDomain3D

    //// the momentum equilibrium equation and state equation are implemented as in Monaghan's paper (1994), simulate free surface flow using sph
    //// here the neighboring list of SPH particles is searched by the cells, vector< vector<sph::SPHParticle*> > SPHParticleCellVec;
    void SmoothParticleHydrodynamics::calculateSPHDensityDotVelocityDotLinkedList3D(){

  // divide the SPH domain into different cells, each cell will contain SPH particles within it
  divideSPHDomain3D();

//  checkDivision();  // pass, May 22, 2015
//  checkNeighborCells3D();  // pass, May 22, 2015

  // initialize the densityDot and velocityDot of all the SPH particles 
  dem::Vec tmp_vec = dem::Vec(0,0, -(util::getParam<>("gravAccel"])*(util::getParam<>("gravScale"]) );
  dem::Vec zero_vec = dem::Vec(0,0,0);
  // in the calculation of forces between sph particles, we need to use the mergeSPHParticleVec, 
  // since mergeVec contains the sph particles from neighboring cpus. While we can only update and migrate and communicate SPHParticleVec
  for(std::vector<sph::SPHParticle*>::iterator pt=mergeSPHParticleVec.begin(); pt!=mergeSPHParticleVec.end(); pt++){
      (*pt)->setDensityDotVelocityDotZero();
      (*pt)->calculateParticleViscosity();  // everytime when use getParticleViscolisity(), make sure that pressure has been calculated!!!!
      (*pt)->calculateParticlePressure();    // everytime when use getParticlePressure(), make sure that pressure has been calculated!!!!
      switch((*pt)->getType()){
        case 1: // free sph particle
    (*pt)->addVelocityDot(tmp_vec);
    break;
        case 2: // ghost sph particle
    break;
        case 3: // boundary sph particle
    (*pt)->addVelocityDot(zero_vec);
    break;
        default:
    std::cout << "SPH particle type of pta should be 1, 2 or 3!" << std::endl;
    exit(-1);
      } // switch
  }

  std::vector<sph::SPHParticle*>::iterator gt;
  for(std::vector<Particle*>::iterator pdt=mergeParticleVec.begin(); pdt!=mergeParticleVec.end(); pdt++){  // all sph ghost particles
      for(gt=(*pdt)->SPHGhostParticleVec.begin(); gt!=(*pdt)->SPHGhostParticleVec.end(); gt++){
        (*gt)->setDensityDotVelocityDotZero();
        (*gt)->calculateParticleViscosity();  // everytime when use getParticleViscolisity(), make sure that pressure has been calculated!!!!
        (*gt)->calculateParticlePressure();  // everytime when use getParticlePressure(), make sure that pressure has been calculated!!!!
      }
  }

  // temporary variables used in the loop
  dem::Vec pta_position;
  dem::Vec ptb_position;
  dem::Vec delta_aWab;
  dem::Vec delta_bWba;
  dem::Vec vab;
  dem::Vec vba;
  dem::Vec vdem;
  dem::Vec dva_dt;
  dem::Vec dvb_dt;
  dem::Vec delta_a;
  dem::Vec delta_b;
  REAL pa, pb, rhoa, rhob, mua, mub;
  REAL rab;
        REAL dWab_dra;
  REAL dWba_drb;
  REAL Wab, Wba;
  REAL da, dB;
  REAL beta;
  REAL xa, ya, xB, yB, k, sqrt_xaya;  // variables for Morris' method to calculate da/dB
  std::vector<sph::SPHParticle*>::iterator ptb;
  REAL ra, rb, rc;  // the geometry of the dem particle
  dem::Vec pta_local, ptb_local;  // local position of pta and ptb in the dem particle
  dem::Vec demt_curr;
  REAL Gamma_ab, mu_ab, vr_dot;
  REAL alpha = util::getParam<>("alpha"];
     REAL alpha_zero = 0;  // the viscous between free and ghost/boundary
  REAL epsilon = util::getParam<>("epsilon"];  // parameter for velocity correction
  REAL Wq, Ra, Rb, phi_4, coefficient_a, coefficient_b;
  dem::Particle* demt;
  
  int pnum;
  std::vector<sph::SPHParticle*> tmp_particleVec;  // sph particles in neighboring cells
  for(int pvec=0; pvec<SPHParticleCellVec.size(); pvec++){

      // store all SPH particles in pvec's neighboring cells
      tmp_particleVec.clear();  // it is the same for the particles in the same cell
      for(std::vector<int>::const_iterator pint=pnum_vec.begin(); pint!=pnum_vec.end(); pint++){
    pnum = pvec+(*pint);
    if(pnum<numCell){
        for(std::vector<sph::SPHParticle*>::iterator pt = SPHParticleCellVec[pnum].begin(); pt!= SPHParticleCellVec[pnum].end(); pt++) {
          tmp_particleVec.push_back(*pt);
                }
    }
      }

      for(std::vector<sph::SPHParticle*>::iterator pta=SPHParticleCellVec[pvec].begin(); pta!=SPHParticleCellVec[pvec].end(); pta++){  // SPH particles in cell pvec
    pa = (*pta)->getParticlePressure();
        if(pa>=0){
        Ra = 0.006;
        }
        else{
        Ra = 0.6;
        }
//        mua = (*pta)->getParticleViscosity();
        rhoa = (*pta)->getParticleDensity();
        pta_position = (*pta)->getCurrPosition();

//    if((*pta)->getType()!=1){  // pta is not free sph particles, i.e. pta is ghost or boundary particles
//        continue;  // do not consider pta, as we treat before
//    }    

    for(ptb=pta+1; ptb!=SPHParticleCellVec[pvec].end(); ptb++){  // sum over the SPH particles in the same cell pvec
        ptb_position = (*ptb)->getCurrPosition();
        rab = dem::vfabs(pta_position-ptb_position);
        if(rab<=kernelSize){  // ptb is in the smooth kernel
          pb = (*ptb)->getParticlePressure();
//          mub = (*ptb)->getParticleViscosity();
          rhob = (*ptb)->getParticleDensity();

          Wq = kernelFunction(rab/smoothLength);
          phi_4 = pow(Wq/Wqmin,4);
          if(pb>=0){
          Rb = 0.006;
          }
          else{
          Rb = 0.6;
          }
          coefficient_a = 1+Ra*phi_4;
          coefficient_b = 1+Rb*phi_4;

      // we have three types of SPH particles: 1, free particle; 2, ghost particle; 3, boundary particle
      // Then we have 3x3 = 9 different types of interactions with the three types of particles
      // so we cannot judge the interaction type only by the type of particle ptb, we need to consider pta also
      // pta          ptb          need to consider or not
      // 1            1              V      free with free
      // 1            2              V      free with ghost
      // 1            3              V      free with boundary

      // 2            1              V      ghost with free
      // 2            2              X      ghost with ghost
      // 2            3              X      ghost with boundary

      // 3            1              V      boundary with free
      // 3            2              X      boundary with ghost
      // 3            3              X       boundary with boundary


      // add the density dot for pta and ptb
          delta_aWab = gradientKernelFunction(pta_position, ptb_position);  // this is to add SPH pta 
          delta_bWba = - delta_aWab;
          dWab_dra = partialKernelFunction(pta_position, ptb_position);  // this is to add SPH pta
          dWba_drb = dWab_dra;
          Wab = kernelFunction(pta_position, ptb_position);  // this is to add SPH pta
          Wba = Wab; 

      switch((*ptb)->getType()){
        case 1:  // ptb is free SPH particle
          switch((*pta)->getType()){
            case 1:  // free with free
                vab = (*pta)->getVelocity()-(*ptb)->getVelocity();
               vba = -vab;

                (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab*delta_aWab) );
                (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba*delta_bWba) );

                vr_dot = vab*(pta_position-ptb_position);
                if(vr_dot<0){
            mu_ab = smoothLength*vr_dot/(rab*rab+0.01*smoothLength*smoothLength);
                Gamma_ab = (-alpha*(util::getParam<>("soundSpeed"])*mu_ab)/(rhoa+rhob)*2;
                }
                else{
             Gamma_ab = 0;
                }


                dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
                (*pta)->addVelocityDot(dva_dt);
                dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)*coefficient_b+pa/(rhoa*rhoa)*coefficient_a+Gamma_ab)*delta_bWba;
                (*ptb)->addVelocityDot(dvb_dt);

                delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
                (*pta)->addVelocityCorrection(delta_a);
                delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
                (*ptb)->addVelocityCorrection(delta_b);
              break;
            case 2:  // ghost with free
        demt = (*pta)->getDemParticle();
        demt_curr = demt->getCurrPos();
            ptb_local = demt->globalToLocal(ptb_position-demt_curr);  // the local position of sph point ptb
            ra = demt->getA(); rb = demt->getB(); rc = demt->getC();
            k = 1.0/(sqrt(ptb_local.getX()*ptb_local.getX()/(ra*ra)+ptb_local.getY()*ptb_local.getY()/(rb*rb)+ptb_local.getZ()*ptb_local.getZ()/(rc*rc) ) );
            da = dem::vfabs(ptb_local-k*ptb_local);  // the distance is the same in rotated coordinates

                // calculate Vab as the method shown in Morris's paper, 1996

                // (1) here I wanna use the distance from the point a/b to the surface of the ellipsoid to simplify the problem
            pta_local = (*pta)->getLocalPosition();
            k = 1.0/(sqrt(pta_local.getX()*pta_local.getX()/(ra*ra)+pta_local.getY()*pta_local.getY()/(rb*rb)+pta_local.getZ()*pta_local.getZ()/(rc*rc) ) );
                dB = dem::vfabs(pta_local-k*pta_local);

//              // (2) here the Morris's method is used 
//              xa = pta_position.getX(); ya = pta_position.getY();
//              xB = ptb_position.getX(); yB = ptb_position.getY();
//              if(ya==0) {da = xa-radius; dB = radius-xB;}
//              else if(xa==0) {da = ya-radius; dB = radius-yB;}
//              else {
//          sqrt_xaya = sqrt(xa*xa+ya*ya);
//          k = radius/sqrt_xaya;
//          da = radius/k - radius;
//          dB = fabs(xa*xB+ya*yB-k*radius*radius)/sqrt_xaya;
//              }

                beta = 1+dB/da;
                if(beta>2 || beta<0 || isnan(beta)){ beta = 2; }

            vdem = (*pta)->getVelocity();
                vba = beta*((*ptb)->getVelocity()-vdem);
                vab = -vba;

                (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab*delta_aWab) );
                (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba*delta_bWba) );

                vr_dot = vab*(pta_position-ptb_position);
                if(vr_dot<0){
                mu_ab = smoothLength*vr_dot/(rab*rab+0.01*smoothLength*smoothLength);
                    Gamma_ab = (-alpha_zero*(util::getParam<>("soundSpeed"])*mu_ab)/(rhoa+rhob)*2;
                }
                else{
                 Gamma_ab = 0;
                }

                dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
            demt->addForce((*pta)->getParticleMass()*dva_dt);
            demt->addMoment( (pta_position-demt_curr) % ((*pta)->getParticleMass()*dva_dt) );
//                (*pta)->addVelocityDot(dva_dt);

                dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)*coefficient_b+pa/(rhoa*rhoa)*coefficient_a+Gamma_ab)*delta_bWba;
                (*ptb)->addVelocityDot(dvb_dt);  // the velocities of the ghost particles will not envolve as the same way as others

//                  delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
//                  (*pta)->addVelocityCorrection(delta_a);
                  delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
                  (*ptb)->addVelocityCorrection(delta_b);
        break;

            case 3:  // boundary with free

        // calculate Vab as the method shown in Morris's paper, 1996
                // interact with boundary particles
                da = ptb_position.getZ()-allContainer.getMinCorner().getZ();  // assume with the bottom boundary
        dB = allContainer.getMinCorner().getZ()-pta_position.getX();  // assume with the bottom boundary
        if(pta_position.getX()<allContainer.getMinCorner().getX()){  // with left boundary
            da = ptb_position.getX()-allContainer.getMinCorner().getX();
            dB = allContainer.getMinCorner().getX()-pta_position.getX();
        }
        else if(pta_position.getX()>allContainer.getMaxCorner().getX()){ // with right boundary
            da = allContainer.getMaxCorner().getX()-ptb_position.getX();
            dB = pta_position.getX()-allContainer.getMaxCorner().getX();
        }
        else if(pta_position.getZ()>allContainer.getMaxCorner().getZ()){ // with top boundary
            da = allContainer.getMaxCorner().getZ()-ptb_position.getZ();
            dB = pta_position.getZ()-allContainer.getMaxCorner().getZ();
        }

                beta = 1+dB/da;
                if(beta>2 || beta<0 || isnan(beta)){ beta = 2; }

                vba = beta*(*ptb)->getVelocity();
                vab = -vba;

                (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab*delta_aWab) );
                (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba*delta_bWba) );

                vr_dot = vab*(pta_position-ptb_position);
                if(vr_dot<0){
            mu_ab = smoothLength*vr_dot/(rab*rab+0.01*smoothLength*smoothLength);
                Gamma_ab = (-alpha_zero*(util::getParam<>("soundSpeed"])*mu_ab)/(rhoa+rhob)*2;
                }
                else{
             Gamma_ab = 0;
                }

//                dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
//                (*pta)->addVelocityDot(dva_dt);
                dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)+pa/(rhoa*rhoa)+Gamma_ab)*delta_bWba;
                (*ptb)->addVelocityDot(dvb_dt);  // the velocities of the ghost particles will not envolve as the same way as others

//            delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
//                (*pta)->addVelocityCorrection(delta_a);
                delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
                (*ptb)->addVelocityCorrection(delta_b);
      
//            // apply the boundary forces by Lennard-Jones potential as in Monaghan's paper(1994)
//                if(rab<=spaceInterval){ // ptb is in the smooth kernel
//                dvb_dt = D*(pow(spaceInterval/rab, p1)-pow(spaceInterval/rab, p2))*(ptb_position-pta_position)/(rab*rab);
//                (*ptb)->addVelocityDot(dvb_dt);
//                } // end if

        break;
            default:
            std::cout << "SPH particle type of pta should be 1, 2 or 3!" << std::endl;
            exit(-1);
            } // end switch pta

          break;
        case 2:  // ptb is ghost particle

          if((*pta)->getType()!=1){  // pta is not free sph particles, i.e. pta is ghost or boundary particles
            break;  // do not consider pta, as we treat before
          }
          demt = (*ptb)->getDemParticle();
          demt_curr = demt->getCurrPos();
          pta_local = demt->globalToLocal(pta_position-demt_curr);  // the local position of sph point pta
          ra = demt->getA(); rb = demt->getB(); rc = demt->getC();
          k = 1.0/(sqrt(pta_local.getX()*pta_local.getX()/(ra*ra)+pta_local.getY()*pta_local.getY()/(rb*rb)+pta_local.getZ()*pta_local.getZ()/(rc*rc) ) );
          da = dem::vfabs(pta_local-k*pta_local);  // the distance is the same in rotated coordinates

              // calculate Vab as the method shown in Morris's paper, 1996

              // (1) here I wanna use the distance from the point a/b to the surface of the ellipsoid to simplify the problem
          ptb_local = (*ptb)->getLocalPosition();
          k = 1.0/(sqrt(ptb_local.getX()*ptb_local.getX()/(ra*ra)+ptb_local.getY()*ptb_local.getY()/(rb*rb)+ptb_local.getZ()*ptb_local.getZ()/(rc*rc) ) );
              dB = dem::vfabs(ptb_local-k*ptb_local);

//              // (2) here the Morris's method is used 
//              xa = pta_position.getX(); ya = pta_position.gety();
//              xB = ptb_position.getX(); yB = ptb_position.gety();
//              if(ya==0) {da = xa-radius; dB = radius-xB;}
//              else if(xa==0) {da = ya-radius; dB = radius-yB;}
//              else {
//          sqrt_xaya = sqrt(xa*xa+ya*ya);
//          k = radius/sqrt_xaya;
//          da = radius/k - radius;
//          dB = fabs(xa*xB+ya*yB-k*radius*radius)/sqrt_xaya;
//              }


              beta = 1+dB/da;
              if(beta>2 || beta<0 || isnan(beta)){ beta = 2; }

          vdem = (*ptb)->getVelocity();
              vab = beta*((*pta)->getVelocity()-vdem);
              vba = -vab;

              (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab*delta_aWab) );
              (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba*delta_bWba) );

              vr_dot = vab*(pta_position-ptb_position);
              if(vr_dot<0){
            mu_ab = smoothLength*vr_dot/(rab*rab+0.01*smoothLength*smoothLength);
                Gamma_ab = (-alpha_zero*(util::getParam<>("soundSpeed"])*mu_ab)/(rhoa+rhob)*2;
              }
              else{
             Gamma_ab = 0;
              }

              dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
              (*pta)->addVelocityDot(dva_dt);
              dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)*coefficient_b+pa/(rhoa*rhoa)*coefficient_a+Gamma_ab)*delta_bWba;
          demt->addForce((*ptb)->getParticleMass()*dvb_dt);
          demt->addMoment( (ptb_position-demt_curr) % ((*ptb)->getParticleMass()*dvb_dt) );
//              (*ptb)->addVelocityDot(dvb_dt);  // the velocities of the ghost particles will not envolve as the same way as others

                delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
                (*pta)->addVelocityCorrection(delta_a);
//                delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
//                (*ptb)->addVelocityCorrection(delta_b);
          break;
        case 3:  // ptb is boundary particle

          if((*pta)->getType()!=1){  // pta is not free sph particles, i.e. pta is ghost or boundary particles
            break;  // do not consider pta, as we treat before
          }

          // calculate Vab as the method shown in Morris's paper, 1996
              // interact with boundary particles
              da = pta_position.getZ()-allContainer.getMinCorner().getZ();  // assume with the bottom boundary
          dB = allContainer.getMinCorner().getZ()-ptb_position.getZ();  // assume with the bottom boundary
          if(ptb_position.getX() < allContainer.getMinCorner().getX()){  // with left boundary
        da = pta_position.getX() - allContainer.getMinCorner().getX();
        dB = allContainer.getMinCorner().getX() - ptb_position.getX();
          }
          else if(ptb_position.getX() > allContainer.getMaxCorner().getX()){  // with right boundary
        da = allContainer.getMaxCorner().getX() - pta_position.getX();
        dB = ptb_position.getX() - allContainer.getMaxCorner().getX();
          }
          else if(ptb_position.getZ() > allContainer.getMaxCorner().getZ()){  // with top boundary
        da = allContainer.getMaxCorner().getZ()-pta_position.getZ();
        dB = ptb_position.getZ()-allContainer.getMaxCorner().getZ();
          }

              beta = 1+dB/da;
              if(beta>2 || beta<0 || isnan(beta)){ beta = 2; }

              vab = beta*(*pta)->getVelocity();
              vba = -vab;

              (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab*delta_aWab) );
              (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba*delta_bWba) );

              vr_dot = vab*(pta_position-ptb_position);
              if(vr_dot<0){
        mu_ab = smoothLength*vr_dot/(rab*rab+0.01*smoothLength*smoothLength);
            Gamma_ab = (-alpha_zero*(util::getParam<>("soundSpeed"])*mu_ab)/(rhoa+rhob)*2;
              }
              else{
         Gamma_ab = 0;
              }

              dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
              (*pta)->addVelocityDot(dva_dt);
//              dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)+pa/(rhoa*rhoa)+Gamma_ab)*delta_bWba;
//              (*ptb)->addVelocityDot(dvb_dt);  // the velocities of the ghost particles will not envolve as the same way as others

              delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
              (*pta)->addVelocityCorrection(delta_a);
//              delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
//              (*ptb)->addVelocityCorrection(delta_b);
      
//          // apply the boundary forces by Lennard-Jones potential as in Monaghan's paper(1994)
//              if(rab<=spaceInterval){ // ptb is in the smooth kernel
//            dva_dt = D*(pow(spaceInterval/rab, p1)-pow(spaceInterval/rab, p2))*(pta_position-ptb_position)/(rab*rab);
//            (*pta)->addVelocityDot(dva_dt);
//              } // end if
          break;
        default:
          std::cout << "SPH particle type should be 1, 2 or 3!" << std::endl;
          exit(-1);
    
      } //end swtich type
        } // end if 3h
    } // end for ptb in the same cell

    for(ptb=tmp_particleVec.begin(); ptb!=tmp_particleVec.end(); ptb++){  // all particles in pvec's neighboring cells 
        ptb_position = (*ptb)->getCurrPosition();
        rab = dem::vfabs(pta_position-ptb_position);
        if(rab<=kernelSize){  // ptb is in the smooth kernel
          pb = (*ptb)->getParticlePressure();
//          mub = (*ptb)->getParticleViscosity();
          rhob = (*ptb)->getParticleDensity();

          Wq = kernelFunction(rab/smoothLength);
          phi_4 = pow(Wq/Wqmin,4);
          if(pb>=0){
          Rb = 0.006;
          }
          else{
          Rb = 0.6;
          }
          coefficient_a = 1+Ra*phi_4;
          coefficient_b = 1+Rb*phi_4;

      // add the density dot for pta and ptb
          delta_aWab = gradientKernelFunction(pta_position, ptb_position);  // this is to add SPH pta 
          delta_bWba = - delta_aWab;
          dWab_dra = partialKernelFunction(pta_position, ptb_position);  // this is to add SPH pta
          dWba_drb = dWab_dra;
          Wab = kernelFunction(pta_position, ptb_position);  // this is to add SPH pta
          Wba = Wab; 
      switch((*ptb)->getType()){
        case 1:  // ptb is free SPH particle
          switch((*pta)->getType()){
            case 1:  // free with free
                vab = (*pta)->getVelocity()-(*ptb)->getVelocity();
               vba = -vab;

                (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab*delta_aWab) );
                (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba*delta_bWba) );

                vr_dot = vab*(pta_position-ptb_position);
                if(vr_dot<0){
            mu_ab = smoothLength*vr_dot/(rab*rab+0.01*smoothLength*smoothLength);
                Gamma_ab = (-alpha*(util::getParam<>("soundSpeed"])*mu_ab)/(rhoa+rhob)*2;
                }
                else{
             Gamma_ab = 0;
                }


                dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
                (*pta)->addVelocityDot(dva_dt);
                dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)*coefficient_b+pa/(rhoa*rhoa)*coefficient_a+Gamma_ab)*delta_bWba;
                (*ptb)->addVelocityDot(dvb_dt);

                delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
                (*pta)->addVelocityCorrection(delta_a);
                delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
                (*ptb)->addVelocityCorrection(delta_b);
              break;
            case 2:  // ghost with free

        demt = (*pta)->getDemParticle();
        demt_curr = demt->getCurrPos();
            ptb_local = demt->globalToLocal(ptb_position-demt_curr);  // the local position of sph point ptb
            ra = demt->getA(); rb = demt->getB(); rc = demt->getC();
            k = 1.0/(sqrt(ptb_local.getX()*ptb_local.getX()/(ra*ra)+ptb_local.getY()*ptb_local.getY()/(rb*rb)+ptb_local.getZ()*ptb_local.getZ()/(rc*rc) ) );
            da = dem::vfabs(ptb_local-k*ptb_local);  // the distance is the same in rotated coordinates

                // calculate Vab as the method shown in Morris's paper, 1996

                // (1) here I wanna use the distance from the point a/b to the surface of the ellipsoid to simplify the problem
            pta_local = (*pta)->getLocalPosition();
            k = 1.0/(sqrt(pta_local.getX()*pta_local.getX()/(ra*ra)+pta_local.getY()*pta_local.getY()/(rb*rb)+pta_local.getZ()*pta_local.getZ()/(rc*rc) ) );
                dB = dem::vfabs(pta_local-k*pta_local);

//              // (2) here the Morris's method is used 
//              xa = pta_position.getX(); ya = pta_position.getY();
//              xB = ptb_position.getX(); yB = ptb_position.getY();
//              if(ya==0) {da = xa-radius; dB = radius-xB;}
//              else if(xa==0) {da = ya-radius; dB = radius-yB;}
//              else {
//          sqrt_xaya = sqrt(xa*xa+ya*ya);
//          k = radius/sqrt_xaya;
//          da = radius/k - radius;
//          dB = fabs(xa*xB+ya*yB-k*radius*radius)/sqrt_xaya;
//              }


                beta = 1+dB/da;
                if(beta>2 || beta<0 || isnan(beta)){ beta = 2; }

            vdem = (*pta)->getVelocity();
                vba = beta*((*ptb)->getVelocity()-vdem);
                vab = -vba;

                (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab*delta_aWab) );
                (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba*delta_bWba) );

                vr_dot = vab*(pta_position-ptb_position);
                if(vr_dot<0){
                mu_ab = smoothLength*vr_dot/(rab*rab+0.01*smoothLength*smoothLength);
                    Gamma_ab = (-alpha_zero*(util::getParam<>("soundSpeed"])*mu_ab)/(rhoa+rhob)*2;
                }
                else{
                 Gamma_ab = 0;
                }

                dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
            demt->addForce((*pta)->getParticleMass()*dva_dt);
            demt->addMoment( (pta_position-demt_curr) % ((*pta)->getParticleMass()*dva_dt) );
//                (*pta)->addVelocityDot(dva_dt);

                dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)*coefficient_b+pa/(rhoa*rhoa)*coefficient_a+Gamma_ab)*delta_bWba;
                (*ptb)->addVelocityDot(dvb_dt);  // the velocities of the ghost particles will not envolve as the same way as others

//                  delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
//                  (*pta)->addVelocityCorrection(delta_a);
                  delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
                  (*ptb)->addVelocityCorrection(delta_b);
        break;
            case 3:  // boundary with free

        // calculate Vab as the method shown in Morris's paper, 1996
                // interact with boundary particles
                da = ptb_position.getZ()-allContainer.getMinCorner().getZ();  // assume with the bottom boundary
        dB = allContainer.getMinCorner().getZ()-pta_position.getX();  // assume with the bottom boundary
        if(pta_position.getX()<allContainer.getMinCorner().getX()){  // with left boundary
            da = ptb_position.getX()-allContainer.getMinCorner().getX();
            dB = allContainer.getMinCorner().getX()-pta_position.getX();
        }
        else if(pta_position.getX()>allContainer.getMaxCorner().getX()){ // with right boundary
            da = allContainer.getMaxCorner().getX()-ptb_position.getX();
            dB = pta_position.getX()-allContainer.getMaxCorner().getX();
        }
        else if(pta_position.getZ()>allContainer.getMaxCorner().getZ()){ // with top boundary
            da = allContainer.getMaxCorner().getZ()-ptb_position.getZ();
            dB = pta_position.getZ()-allContainer.getMaxCorner().getZ();
        }

                beta = 1+dB/da;
                if(beta>2 || beta<0 || isnan(beta)){ beta = 2; }

                vba = beta*(*ptb)->getVelocity();
                vab = -vba;

                (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab*delta_aWab) );
                (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba*delta_bWba) );

                vr_dot = vab*(pta_position-ptb_position);
                if(vr_dot<0){
            mu_ab = smoothLength*vr_dot/(rab*rab+0.01*smoothLength*smoothLength);
                Gamma_ab = (-alpha_zero*(util::getParam<>("soundSpeed"])*mu_ab)/(rhoa+rhob)*2;
                }
                else{
             Gamma_ab = 0;
                }

//                dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
//                (*pta)->addVelocityDot(dva_dt);
                dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)+pa/(rhoa*rhoa)+Gamma_ab)*delta_bWba;
                (*ptb)->addVelocityDot(dvb_dt);  // the velocities of the ghost particles will not envolve as the same way as others

//            delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
//                (*pta)->addVelocityCorrection(delta_a);
                delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
                (*ptb)->addVelocityCorrection(delta_b);
      
            // apply the boundary forces by Lennard-Jones potential as in Monaghan's paper(1994)
                if(rab<=spaceInterval){ // ptb is in the smooth kernel
                dvb_dt = D*(pow(spaceInterval/rab, p1)-pow(spaceInterval/rab, p2))*(ptb_position-pta_position)/(rab*rab);
                (*ptb)->addVelocityDot(dvb_dt);
                } // end if

        break;
            default:
            std::cout << "SPH particle type of pta should be 1, 2 or 3!" << std::endl;
            exit(-1);
            } // end switch pta

          break;
        case 2:  // ptb is ghost particle

          if((*pta)->getType()!=1){  // pta is not free sph particles, i.e. pta is ghost or boundary particles
            break;  // do not consider pta, as we treat before
          }

          demt = (*ptb)->getDemParticle();
          demt_curr = demt->getCurrPos();
          pta_local = demt->globalToLocal(pta_position-demt_curr);  // the local position of sph point pta
          ra = demt->getA(); rb = demt->getB(); rc = demt->getC();
          k = 1.0/(sqrt(pta_local.getX()*pta_local.getX()/(ra*ra)+pta_local.getY()*pta_local.getY()/(rb*rb)+pta_local.getZ()*pta_local.getZ()/(rc*rc) ) );
          da = dem::vfabs(pta_local-k*pta_local);  // the distance is the same in rotated coordinates

              // calculate Vab as the method shown in Morris's paper, 1996

              // (1) here I wanna use the distance from the point a/b to the surface of the ellipsoid to simplify the problem
          ptb_local = (*ptb)->getLocalPosition();
          k = 1.0/(sqrt(ptb_local.getX()*ptb_local.getX()/(ra*ra)+ptb_local.getY()*ptb_local.getY()/(rb*rb)+ptb_local.getZ()*ptb_local.getZ()/(rc*rc) ) );
              dB = dem::vfabs(ptb_local-k*ptb_local);

//              // (2) here the Morris's method is used 
//              xa = pta_position.getX(); ya = pta_position.getY();
//              xB = ptb_position.getX(); yB = ptb_position.getY();
//              if(ya==0) {da = xa-radius; dB = radius-xB;}
//              else if(xa==0) {da = ya-radius; dB = radius-yB;}
//              else {
//          sqrt_xaya = sqrt(xa*xa+ya*ya);
//          k = radius/sqrt_xaya;
//          da = radius/k - radius;
//          dB = fabs(xa*xB+ya*yB-k*radius*radius)/sqrt_xaya;
//              }


              beta = 1+dB/da;
              if(beta>2 || beta<0 || isnan(beta)){ beta = 2; }

          vdem = (*ptb)->getVelocity();
              vab = beta*((*pta)->getVelocity()-vdem);
              vba = -vab;

              (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab*delta_aWab) );
              (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba*delta_bWba) );

              vr_dot = vab*(pta_position-ptb_position);
              if(vr_dot<0){
            mu_ab = smoothLength*vr_dot/(rab*rab+0.01*smoothLength*smoothLength);
                Gamma_ab = (-alpha_zero*(util::getParam<>("soundSpeed"])*mu_ab)/(rhoa+rhob)*2;
              }
              else{
             Gamma_ab = 0;
              }

              dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
              (*pta)->addVelocityDot(dva_dt);
              dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)*coefficient_b+pa/(rhoa*rhoa)*coefficient_a+Gamma_ab)*delta_bWba;
          demt->addForce((*ptb)->getParticleMass()*dvb_dt);
          demt->addMoment( (ptb_position-demt_curr) % ((*ptb)->getParticleMass()*dvb_dt) );
//              (*ptb)->addVelocityDot(dvb_dt);  // the velocities of the ghost particles will not envolve as the same way as others

                delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
                (*pta)->addVelocityCorrection(delta_a);
//                delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
//                (*ptb)->addVelocityCorrection(delta_b);
          break;
        case 3:  // ptb is boundary particle

          if((*pta)->getType()!=1){  // pta is not free sph particles, i.e. pta is ghost or boundary particles
            break;  // do not consider pta, as we treat before
          }

          // calculate Vab as the method shown in Morris's paper, 1996
              // interact with boundary particles
              da = pta_position.getZ()-allContainer.getMinCorner().getZ();  // assume with the bottom boundary
          dB = allContainer.getMinCorner().getZ()-ptb_position.getZ();  // assume with the bottom boundary
          if(ptb_position.getX() < allContainer.getMinCorner().getX()){  // with left boundary
        da = pta_position.getX() - allContainer.getMinCorner().getX();
        dB = allContainer.getMinCorner().getX() - ptb_position.getX();
          }
          else if(ptb_position.getX() > allContainer.getMaxCorner().getX()){  // with right boundary
        da = allContainer.getMaxCorner().getX() - pta_position.getX();
        dB = ptb_position.getX() - allContainer.getMaxCorner().getX();
          }
          else if(ptb_position.getZ() > allContainer.getMaxCorner().getZ()){  // with top boundary
        da = allContainer.getMaxCorner().getZ()-pta_position.getZ();
        dB = ptb_position.getZ()-allContainer.getMaxCorner().getZ();
          }

              beta = 1+dB/da;
              if(beta>2 || beta<0 || isnan(beta)){ beta = 2; }

              vab = beta*(*pta)->getVelocity();
              vba = -vab;

              (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab*delta_aWab) );
              (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba*delta_bWba) );

              vr_dot = vab*(pta_position-ptb_position);
              if(vr_dot<0){
        mu_ab = smoothLength*vr_dot/(rab*rab+0.01*smoothLength*smoothLength);
            Gamma_ab = (-alpha_zero*(util::getParam<>("soundSpeed"])*mu_ab)/(rhoa+rhob)*2;
              }
              else{
         Gamma_ab = 0;
              }

              dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
              (*pta)->addVelocityDot(dva_dt);
//              dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)+pa/(rhoa*rhoa)+Gamma_ab)*delta_bWba;
//              (*ptb)->addVelocityDot(dvb_dt);  // the velocities of the ghost particles will not envolve as the same way as others

              delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
              (*pta)->addVelocityCorrection(delta_a);
//              delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
//              (*ptb)->addVelocityCorrection(delta_b);

              // apply the boundary forces by Lennard-Jones potential as in Monaghan's paper(1994)
              if(rab<=spaceInterval){ // ptb is in the smooth kernel
            dva_dt = D*(pow(spaceInterval/rab, p1)-pow(spaceInterval/rab, p2))*(pta_position-ptb_position)/(rab*rab);
            (*pta)->addVelocityDot(dva_dt);
              } // end if
      
          break;
        default:
          std::cout << "SPH particle type should be 1, 2 or 3!" << std::endl;
          exit(-1);
    
      } //end swtich type
        } // end if 3h
    } // end for ptb in neighbor cells
      } // end for pta

      tmp_particleVec.clear();  // clear elements in tmp-vector for particles neighboring cells, it is important

  } // end for pvec, different cells


//      // apply the boundary forces by Lennard-Jones potential as in Monaghan's paper(1994)
//      for(ptb=SPHBoundaryParticleVec.begin(); ptb!=SPHBoundaryParticleVec.end(); ptb++){  
//    ptb_position = (*ptb)->getCurrPosition();
//    rab = dem::vfabs(pta_position-ptb_position);
//    if(rab<=spaceInterval){ // ptb is in the smooth kernel
//        dva_dt = D*(pow(spaceInterval/rab, p1)-pow(spaceInterval/rab, p2))*(pta_position-ptb_position)/(rab*rab);
//        (*pta)->addVelocityDot(dva_dt);
//    } // end if
//      } // end ptb


    } // end calculateSPHDensityDotVelocityDotLinkedList3D()

