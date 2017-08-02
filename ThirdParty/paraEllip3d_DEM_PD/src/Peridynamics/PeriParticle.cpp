// Function Definitions
#include <Peridynamics/PeriParticle.h>
#include <Core/Util/Utility.h>

namespace pd {

using dem::Parameter;

PeriParticle::PeriParticle()
{
  d_id = 0;
  isAlive = true;
  initPosition = dem::Vec(0, 0, 0);
  mass = 0.0;
  volume = 0.0;
  displacement = 0.0;
  velocity = 0.0;
  velocityHalf = 0.0;
  acceleration = 0.0;
  horizonSize = -1.0e16;

  sigma = dem::zeros(3, 3);
  deformationGradient = dem::zeros(3, 3);
  deformationGradientHalf = dem::zeros(3, 3);
  Kinv = dem::zeros(3, 3);
  isv11 = 0;
  if (util::getParam<int>("typeConstitutive") ==
      1) { // 1---implicit, 2---explicit
    isv11 = util::getParam<REAL>("chi");
  } else {
    isv11 = util::getParam<REAL>("c");
  }

  tangentModulus = dem::zeros(6, 6);
  tangentModulus(1, 1) =
    util::getParam<REAL>("tangentModulus11");
  tangentModulus(1, 2) =
    util::getParam<REAL>("tangentModulus12");
  tangentModulus(1, 3) =
    util::getParam<REAL>("tangentModulus13");
  tangentModulus(2, 1) =
    util::getParam<REAL>("tangentModulus21");
  tangentModulus(2, 2) =
    util::getParam<REAL>("tangentModulus22");
  tangentModulus(2, 3) =
    util::getParam<REAL>("tangentModulus23");
  tangentModulus(3, 1) =
    util::getParam<REAL>("tangentModulus31");
  tangentModulus(3, 2) =
    util::getParam<REAL>("tangentModulus32");
  tangentModulus(3, 3) =
    util::getParam<REAL>("tangentModulus33");
  tangentModulus(4, 4) =
    util::getParam<REAL>("tangentModulus44");
  tangentModulus(5, 5) =
    util::getParam<REAL>("tangentModulus55");
  tangentModulus(6, 6) =
    util::getParam<REAL>("tangentModulus66");

  sigma11 = 0;
  sigma12 = 0;
  sigma13 = 0;
  sigma21 = 0;
  sigma22 = 0;
  sigma23 = 0;
  sigma31 = 0;
  sigma32 = 0;
  sigma33 = 0;

  Kinv11 = 0;
  Kinv12 = 0;
  Kinv13 = 0;
  Kinv21 = 0;
  Kinv22 = 0;
  Kinv23 = 0;
  Kinv31 = 0;
  Kinv32 = 0;
  Kinv33 = 0;
  BondedDEMParticleID.clear();
} // end PeriParticle()

PeriParticle::PeriParticle(ParticleID id, REAL x, REAL y, REAL z)
{
  d_id = id;
  isAlive = true;
  initPosition.setX(x);
  initPosition.setY(y);
  initPosition.setZ(z);
  mass = 0.0;
  volume = 0.0;
  displacement = 0.0;
  velocity = 0.0;
  velocityHalf = 0.0;
  acceleration = 0.0;
  horizonSize = -1.0e16;

  sigma = dem::zeros(3, 3);
  deformationGradient = dem::zeros(3, 3);
  deformationGradientHalf = dem::zeros(3, 3);
  Kinv = dem::zeros(3, 3);
  isv11 = 0;
  if (util::getParam<int>("typeConstitutive") ==
      1) { // 1---implicit, 2---explicit
    isv11 = util::getParam<REAL>("chi");
  } else {
    isv11 = util::getParam<REAL>("c");
  }

  tangentModulus = dem::zeros(6, 6);
  tangentModulus(1, 1) =
    util::getParam<REAL>("tangentModulus11");
  tangentModulus(1, 2) =
    util::getParam<REAL>("tangentModulus12");
  tangentModulus(1, 3) =
    util::getParam<REAL>("tangentModulus13");
  tangentModulus(2, 1) =
    util::getParam<REAL>("tangentModulus21");
  tangentModulus(2, 2) =
    util::getParam<REAL>("tangentModulus22");
  tangentModulus(2, 3) =
    util::getParam<REAL>("tangentModulus23");
  tangentModulus(3, 1) =
    util::getParam<REAL>("tangentModulus31");
  tangentModulus(3, 2) =
    util::getParam<REAL>("tangentModulus32");
  tangentModulus(3, 3) =
    util::getParam<REAL>("tangentModulus33");
  tangentModulus(4, 4) =
    util::getParam<REAL>("tangentModulus44");
  tangentModulus(5, 5) =
    util::getParam<REAL>("tangentModulus55");
  tangentModulus(6, 6) =
    util::getParam<REAL>("tangentModulus66");

  sigma11 = 0;
  sigma12 = 0;
  sigma13 = 0;
  sigma21 = 0;
  sigma22 = 0;
  sigma23 = 0;
  sigma31 = 0;
  sigma32 = 0;
  sigma33 = 0;

  Kinv11 = 0;
  Kinv12 = 0;
  Kinv13 = 0;
  Kinv21 = 0;
  Kinv22 = 0;
  Kinv23 = 0;
  Kinv31 = 0;
  Kinv32 = 0;
  Kinv33 = 0;
  BondedDEMParticleID.clear();
} // end PeriParticle()

PeriParticle::PeriParticle(const PeriParticle& pt)
{
  d_id = pt.d_id;
  isAlive = pt.isAlive;
  initPosition = pt.initPosition;
  mass = pt.mass;
  volume = pt.volume;
  displacement = pt.displacement;
  velocity = pt.velocity;
  velocityHalf = pt.velocityHalf;
  acceleration = pt.acceleration;
  sigma = pt.sigma;
  deformationGradient = pt.deformationGradient;
  deformationGradientHalf = pt.deformationGradientHalf;
  Kinv = pt.Kinv;
  isv11 = pt.isv11;
  tangentModulus = pt.tangentModulus;
  horizonSize = pt.horizonSize;
  // in order to keep stress values after gathering
  sigma11 = pt.sigma11;
  sigma12 = pt.sigma12;
  sigma13 = pt.sigma13;

  sigma21 = pt.sigma21;
  sigma22 = pt.sigma22;
  sigma23 = pt.sigma23;

  sigma31 = pt.sigma31;
  sigma32 = pt.sigma32;
  sigma33 = pt.sigma33;

  Kinv11 = pt.Kinv11;
  Kinv12 = pt.Kinv12;
  Kinv13 = pt.Kinv13;

  Kinv21 = pt.Kinv21;
  Kinv22 = pt.Kinv22;
  Kinv23 = pt.Kinv23;

  Kinv31 = pt.Kinv31;
  Kinv32 = pt.Kinv32;
  Kinv33 = pt.Kinv33;
  BondedDEMParticleID = pt.BondedDEMParticleID;
}

PeriParticle::~PeriParticle()
{

  // free the spaces of these pointer vector
  //	    for(std::vector<PeriParticle*>::iterator ip=neighborVec.begin();
  // ip!=neighborVec.end(); ip++){
  //		delete (*ip);
  //	    }
  //	    for(std::vector<PeriBond*>::iterator ib=bondVec.begin();
  // ib!=bondVec.end(); ib++) {
  //		if( (*ib)!=NULL ){
  //		    delete (*ib);
  //		    (*ib)=NULL;
  //	 	}
  //	    }

  //	    neighborVec.clear();
  bondVec.clear();

  sigma.clear();
  deformationGradient.clear();
  deformationGradientHalf.clear();
  Kinv.clear();
  tangentModulus.clear();

} // end PeriParticle()

void
PeriParticle::releaseBondVec()
{

  /*
  for (std::vector<PeriBond *>::iterator ib = bondVec.begin();
       ib != bondVec.end(); ib++) {
    if ((*ib) != NULL) {
      delete (*ib);
      (*ib) = NULL;
    }
  }
  */
  bondVec.clear();
} // releaseBondVec

void
PeriParticle::constructMatrixMember()
{
  sigma.clear();
  deformationGradient.clear();
  deformationGradientHalf.clear();
  Kinv.clear();
  //    	    isv.clear();
  tangentModulus.clear();

  sigma = dem::zeros(3, 3);
  deformationGradient = dem::zeros(3, 3);
  deformationGradientHalf = dem::zeros(3, 3);
  Kinv = dem::zeros(3, 3);
  Kinv(1, 1) = Kinv11;
  Kinv(1, 2) = Kinv12;
  Kinv(1, 3) = Kinv13;
  Kinv(2, 1) = Kinv21;
  Kinv(2, 2) = Kinv22;
  Kinv(2, 3) = Kinv23;
  Kinv(3, 1) = Kinv31;
  Kinv(3, 2) = Kinv32;
  Kinv(3, 3) = Kinv33;
  //	    isv = dem::zeros(1,5);
  //	    if(util::getParam<int>("typeConstitutive") ==
  // 1){
  //// 1---implicit, 2---explicit
  //	    	isv(1,1) = util::getParam<REAL>("chi");
  //	    }
  //	    else{
  //	    	isv(1,1) = util::getParam<REAL>("c");
  //	    }

  tangentModulus = dem::zeros(6, 6);
  tangentModulus(1, 1) =
    util::getParam<REAL>("tangentModulus11");
  tangentModulus(1, 2) =
    util::getParam<REAL>("tangentModulus12");
  tangentModulus(1, 3) =
    util::getParam<REAL>("tangentModulus13");
  tangentModulus(2, 1) =
    util::getParam<REAL>("tangentModulus21");
  tangentModulus(2, 2) =
    util::getParam<REAL>("tangentModulus22");
  tangentModulus(2, 3) =
    util::getParam<REAL>("tangentModulus23");
  tangentModulus(3, 1) =
    util::getParam<REAL>("tangentModulus31");
  tangentModulus(3, 2) =
    util::getParam<REAL>("tangentModulus32");
  tangentModulus(3, 3) =
    util::getParam<REAL>("tangentModulus33");
  tangentModulus(4, 4) =
    util::getParam<REAL>("tangentModulus44");
  tangentModulus(5, 5) =
    util::getParam<REAL>("tangentModulus55");
  tangentModulus(6, 6) =
    util::getParam<REAL>("tangentModulus66");
}

void
PeriParticle::replaceHorizonSizeIfLarger(REAL size)
{
  horizonSize = (horizonSize < size) ? size : horizonSize;
} // end replaceHorizonSizeIfLarger()

void
PeriParticle::calcParticleKinv()
{

  dem::Matrix K(3, 3);
  for (auto& bond : bondVec) {

    // check which pt1 or pt2 in (*bond) is the center, namely (*pt)
    bool is_pt1 = false; // true when (*pt1) is the center
    if (this == bond->getPt1().get()) {
      is_pt1 = true;
    }

    // dem::Vec xi = (*bond)->getXi(is_pt1);
    // K += dyadicProduct(xi,
    // xi)*(*bond)->getVolume(is_pt1)*(*bond)->getWeight();
    K = K + bond->getMicroK(is_pt1);

  } // end bond

  //// for numerical purpose, to be deleted later
  // K = 1.0/(horizonSize*horizonSize)*K;
  //
  //// inverse of matrix K
  // Kinv = K.getInvs()/(horizonSize*horizonSize);

  Kinv = inv(K);
  //std::cout << "Kinv = " << Kinv << "\n";
  assignKinv();

} // end calcParticleKinv()

/*
        void PeriParticle::checkParticleAlive(){

                int num_bonds = 0;	// the number of alive bonds
                for(std::vector<PeriBond*>::iterator bond=bondVec.begin();
bond!=bondVec.end();
bond++){

                    if( (*bond)->getIsAlive() ){
                        REAL bond_length = (*bond)->calcCurrentLength();

                        REAL init_length = (*bond)->getInitLength();
                        REAL stretch = ( bond_length - init_length
)/init_length;

                        if(stretch > stretch_limit || stretch < -2.0 ){
                            (*bond)->setAliveFalse();
                        }
                        else{
                            num_bonds++;
                        }

                    } // if alive
                } // end bond

                // disable a particle
                if(num_bonds < 1){	// as rigid particle
                    isAlive = false;
                    //std::cout << "A particle is disabled due to the lack of
bond" <<
std::endl;
                }

        } // end checkParticleAlive()
*/

void
PeriParticle::checkParticleAlive()
{

  int num_bonds = 0; // the number of alive bonds
  for (auto& bond : bondVec) {
    if (bond->getIsAlive())
      num_bonds++; // if alive
  }                // end bond

  // disable a particle
  if (num_bonds < 1) { // as rigid particle
    isAlive = false;
    //std::cout << "A particle is disabled due to the lack of bond" << std::endl;
  }
} // end checkParticleAlive()

void
PeriParticle::calcParticleStress()
{

  if (!isAlive) { // not alive
    sigma = dem::zeros(3, 3);
  } else {
    // calculate deformation gradient tensor at current and half step
    dem::Matrix N(3, 3);      // matrix N at n+1 step
    dem::Matrix N_half(3, 3); // matrix N at n+1/2 step
    dem::Matrix N_deltaU(3, 3);
    // matrix N, corresponding to \mathbf{u}^{n+1} - \mathbf{u}^{n},
    // used to calculate \nabla (\mathbf{u}^{n+1} - \mathbf{u}^{n})

    for (auto& bond : bondVec) {
      
      // check which pt1 or pt2 in (*bond) is the center, namely (*pt)
      bool is_pt1 = false; // true when (*pt1) is the center
      if (this == bond->getPt1().get()) {
        is_pt1 = true;
      }

      bool bondIsAlive = bond->getIsAlive();

      N = N + bond->getMicroN(is_pt1, bondIsAlive);

      N_half =
        N_half +
        bond->getMicroNHalf(is_pt1, bondIsAlive,
                          util::getParam<REAL>("timeStep"));

      N_deltaU = N_deltaU +
                 bond->getMicroNDeltaU(
                   is_pt1, bondIsAlive,
                   util::getParam<REAL>("timeStep"));

      // if((*bond)->getIsAlive()){

      //	N += (*bond)->getMicroN(is_pt1);

      //	N_half += (*bond)->getMicroNHalf(is_pt1,
      // util::getParam<REAL>("timeStep"));

      //	N_deltaU += (*bond)->getMicroNDeltaU(is_pt1,
      // util::getParam<REAL>("timeStep"));

      //      }

    } // end bond

    deformationGradient = N * Kinv;
    deformationGradientHalf = N_half * Kinv;

    /*
    std::cout << "Particle = " << d_id << " N = " << N << "\n\t Kinv = " << Kinv
              << "\n\t DefGrad = " << deformationGradient << "\n";
    */

    REAL eps = 1.0e-2;
    if (det(deformationGradient) < eps || det(deformationGradientHalf) < eps) {
      // calculate the determinant of deformationGraident and
      // deformationGradientHalf,
      // if the determinants are too small, then this particle is disabled,
      // isAlive = false
      isAlive = false; // disabled particle
      sigma = dem::zeros(3, 3);
      //std::cout << "A particle is disabled because det[F] < 0.0" << std::endl;
    } else {
      if (util::getParam<int>("typeConstitutive") == 1) {
        // Linear Elasticity, for testing purpose
        dem::Matrix identity3x3(3, 3);
        identity3x3(1, 1) = 1;
        identity3x3(2, 2) = 1;
        identity3x3(3, 3) = 1;
        dem::Matrix dudx =
          (deformationGradient - identity3x3) * (inv(deformationGradient));
        dem::Matrix voight_strain(6, 1);
        voight_strain(1, 1) = dudx(1, 1);
        voight_strain(2, 1) = dudx(2, 2);
        voight_strain(3, 1) = dudx(3, 3);
        voight_strain(4, 1) = dudx(2, 3) + dudx(3, 2);
        voight_strain(5, 1) = dudx(1, 3) + dudx(3, 1);
        voight_strain(6, 1) = dudx(1, 2) + dudx(2, 1);
        dem::Matrix voight_sigma = tangentModulus * voight_strain;
        sigma(1, 1) = voight_sigma(1, 1);
        sigma(2, 2) = voight_sigma(2, 1);
        sigma(3, 3) = voight_sigma(3, 1);
        sigma(2, 3) = voight_sigma(4, 1);
        sigma(1, 3) = voight_sigma(5, 1);
        sigma(1, 2) = voight_sigma(6, 1);
        sigma(2, 1) = sigma(1, 2);
        sigma(3, 1) = sigma(1, 3);
        sigma(3, 2) = sigma(2, 3);
      } else if (util::getParam<int>("typeConstitutive") ==
                 2) {
        // calculate G, \nabla \Delta \bf{u}
        dem::Matrix G = N_deltaU * Kinv * inv(deformationGradientHalf);
        dem::Matrix Gsymm = 0.5 * (G + trans(G)); // symmetric part of G
        dem::Matrix Gskew = 0.5 * (G - trans(G)); // skew part of G

        dem::Matrix voight_Gsymm(6, 1);
        voight_Gsymm(1, 1) = Gsymm(1, 1);
        voight_Gsymm(2, 1) = Gsymm(2, 2);
        voight_Gsymm(3, 1) = Gsymm(3, 3);
        voight_Gsymm(4, 1) = Gsymm(2, 3);
        voight_Gsymm(5, 1) = Gsymm(1, 3);
        voight_Gsymm(6, 1) = Gsymm(1, 2);

        dem::Matrix voight_delta_sigma = tangentModulus * voight_Gsymm;

        dem::Matrix delta_sigma(3, 3);
        delta_sigma(1, 1) = voight_delta_sigma(1, 1);
        delta_sigma(2, 2) = voight_delta_sigma(2, 1);
        delta_sigma(3, 3) = voight_delta_sigma(3, 1);
        delta_sigma(2, 3) = voight_delta_sigma(4, 1);
        delta_sigma(1, 3) = voight_delta_sigma(5, 1);
        delta_sigma(1, 2) = voight_delta_sigma(6, 1);
        delta_sigma(2, 1) = delta_sigma(1, 2);
        delta_sigma(3, 1) = delta_sigma(1, 3);
        delta_sigma(3, 2) = delta_sigma(2, 3);

        dem::Matrix identity3x3(3, 3);
        identity3x3(1, 1) = 1;
        identity3x3(2, 2) = 1;
        identity3x3(3, 3) = 1;

        dem::Matrix Q = identity3x3 + inv(identity3x3 - 0.5 * Gskew) * Gskew;
        dem::Matrix trial_sigma = Q * sigma * trans(Q) + delta_sigma;

        // calculate deviatoric trial stress
        dem::Matrix deviatoric_trial_sigma = trial_sigma;
        REAL trace_trial_sigma =
          trial_sigma(1, 1) + trial_sigma(2, 2) + trial_sigma(3, 3);
        deviatoric_trial_sigma(1, 1) =
          trial_sigma(1, 1) - 1.0 / 3.0 * trace_trial_sigma;
        deviatoric_trial_sigma(2, 2) =
          trial_sigma(2, 2) - 1.0 / 3.0 * trace_trial_sigma;
        deviatoric_trial_sigma(3, 3) =
          trial_sigma(3, 3) - 1.0 / 3.0 * trace_trial_sigma;

        REAL L2norm_deviatoric_trial_sigma = 0.0;
        L2norm_deviatoric_trial_sigma =
          deviatoric_trial_sigma(1, 1) * deviatoric_trial_sigma(1, 1) +
          deviatoric_trial_sigma(2, 2) * deviatoric_trial_sigma(2, 2) +
          deviatoric_trial_sigma(3, 3) * deviatoric_trial_sigma(3, 3) +
          2.0 * (deviatoric_trial_sigma(1, 2) * deviatoric_trial_sigma(1, 2)) +
          2.0 * (deviatoric_trial_sigma(1, 3) * deviatoric_trial_sigma(1, 3)) +
          2.0 * (deviatoric_trial_sigma(2, 3) * deviatoric_trial_sigma(2, 3));
        L2norm_deviatoric_trial_sigma = sqrt(L2norm_deviatoric_trial_sigma);

        REAL Aphi = util::getParam<REAL>("Aphi");
        REAL Bphi = util::getParam<REAL>("Bphi");
        REAL cn = isv11; //?
        REAL f_trial = L2norm_deviatoric_trial_sigma -
                       (Aphi * cn - Bphi * 1.0 / 3.0 * trace_trial_sigma);

        if (f_trial < 0) { // elasticity
          sigma = trial_sigma;
          isv11 = cn;
        } else { // plasticity

          REAL Bpsi = util::getParam<REAL>("Bpsi");
          REAL Hc = util::getParam<REAL>("hci");
          REAL KBulk = util::getParam<REAL>("kBulk");
          REAL mu = util::getParam<REAL>("mu");
          REAL delta_gamma =
            f_trial / (2.0 * mu + KBulk * Bphi * Bpsi + Hc * Aphi * Aphi);
          sigma = trial_sigma -
                  delta_gamma * (KBulk * Bpsi * identity3x3 +
                                 2.0 * mu * deviatoric_trial_sigma /
                                   L2norm_deviatoric_trial_sigma);
          isv11 = cn + delta_gamma * Hc * Aphi;
        }
      }

    } // alive particle

  } //

  //std::cout << "ParticleID = " << d_id << " Stress: " << sigma << "\n";

} // end calcParticleStress()

void
PeriParticle::calcParticleAcceleration()
{

  acceleration = 0.0;
  dem::Matrix acceleration_matrix(3, 1);
  dem::Matrix xi_ik_matrix(3, 1);
  dem::Matrix PSi;
  dem::Matrix PSk;
  if (isAlive) {
    for (auto& bond : bondVec) {
      PeriParticle* pti;
      PeriParticle* ptk;
      if (this == bond->getPt1().get()) {
        pti = bond->getPt1().get();
        ptk = bond->getPt2().get();
      } else {
        pti = bond->getPt2().get();
        ptk = bond->getPt1().get();
      }

      /*
      std::cout << "Particle i = " << pti->getId()
                << " def grad = " << pti->deformationGradient << "\n";
      std::cout << "Particle k = " << ptk->getId()
                << " def grad = " << ptk->deformationGradient << "\n";
      */

      // Piola Kirchoff stress of particle i
      PSi = det(pti->deformationGradient) * pti->sigma *
            inv(trans(pti->deformationGradient));
      // Piola Kirchoff stress of particle k
      PSk = det(ptk->deformationGradient) * ptk->sigma *
            inv(trans(ptk->deformationGradient));

      dem::Vec xi_ik = ptk->initPosition - pti->initPosition;
      xi_ik_matrix(1, 1) = xi_ik.x();
      xi_ik_matrix(2, 1) = xi_ik.y();
      xi_ik_matrix(3, 1) = xi_ik.z();

      acceleration_matrix = acceleration_matrix +
                            bond->getWeight() *
                              (PSi * (pti->Kinv) + PSk * (ptk->Kinv)) *
                              xi_ik_matrix * ptk->volume;

    } // end bond
    dem::Matrix grav_vec(3, 1);
    grav_vec(1, 1) = 0;
    grav_vec(2, 1) = 0;
    grav_vec(3, 1) = -util::getParam<REAL>("gravAccel") *
                     (util::getParam<REAL>("gravScale"));

    // actually, if we apply the background damping, then this is acceleration
    // is not the
    // real one. rho0*a_(n+1) = intergral + rho0*b = f_(n+1), the acceleration
    // here is just
    // f_(n+1)/rho0 (we can denote it as acce_here), which is the real
    // acceleration if we do not consider background damping.
    // However, if consider background damping, then rho0*a_(n+1) =
    // f_(n+1)-alpha*rho0*v_(n+1) && from step3 of velocity-verlet integration
    // namely equation (17) in Houfu's note: v_(n+1) = v_(n+1/2)+1/2*a_(n+1)*dt,
    // substitue this into above force equilibrium equation, we
    // get (2+alpha*dt)*v_(n+1) = 2*v_(n+1/2)+f_(n+1)*dt/rho0, where
    // f_(n+1)/rho0 = acce_here. From this equation, we can solve v_(n+1)
    // then go back we can get a_(n+1) = f_(n+1)/rho0-alpha*rho0*v_(n+1).
    // here in this function, we will keep this to calculate acce_here only. the
    // v_(n+1) and a_(n+1) will be calculated in updateVelocity()
    acceleration_matrix =
      acceleration_matrix /
      (util::getParam<REAL>("periDensity") *
       util::getParam<REAL>("massScale")) /*+grav_vec*/;
    acceleration.setX(acceleration_matrix(1, 1));
    acceleration.setY(acceleration_matrix(2, 1));
    acceleration.setZ(acceleration_matrix(3, 1));

  } // alive particle

  //std::cout << "ParticleID = " << d_id << " Acc: " << acceleration << "\n";
} // end calcParticleAcceleration()

void
PeriParticle::updateDisplacement()
{
  REAL deltaT = util::getParam<REAL>("timeStep");
  velocityHalf = velocity + 0.5 * acceleration * deltaT;
  prevDisp = displacement;
  displacement += velocityHalf * deltaT;

  /*
  std::cout << "Disp: P=" << d_id << " delT = " << deltaT 
            << " acc = " << acceleration
            << " v_n = " << velocity << " v_n+1/2 = " << velocityHalf
            << " u_n = " << prevDisp << " u_n+1 = " << displacement << "\n";
  */

} // end updateDisplacement()

void
PeriParticle::updateVelocity()
{
  REAL deltaT = util::getParam<REAL>("timeStep");
  REAL forceDamp = util::getParam<REAL>("forceDamp");

  REAL atf = 2.0 + forceDamp*deltaT;

  // here acceleration is not the real a_(n+1), it is f_(n+1)/rho0
  //	    acceleration = acceleration-dem::DMP_F*velocity;
  velocity = 2.0 * velocityHalf / atf + acceleration * deltaT / atf; 
  acceleration = 2.0 * (velocity - velocityHalf) / deltaT;

  /*
  std::cout << "Vel: P=" << d_id << " delT = " << deltaT << " atf = " << atf
            << " acc = " << acceleration
            << " v_n+1 = " << velocity << " v_n+1/2 = " << velocityHalf << "\n";
  */

} // end updateVelocity()

void
PeriParticle::initial()
{
  displacement = 0.0;
  velocity = 0.0;
  velocityHalf = 0.0;
  acceleration = 0.0;
  sigma = dem::zeros(3, 3);
  deformationGradient = dem::zeros(3, 3);
  deformationGradientHalf = dem::zeros(3, 3);
  tangentModulus = dem::zeros(6, 6);
  isv11 = 0;
  if (util::getParam<int>("typeConstitutive") ==
      1) { // 1---implicit, 2---explicit
    isv11 = util::getParam<REAL>("chi");
  } else {
    isv11 = util::getParam<REAL>("c");
  }
  tangentModulus = dem::zeros(6, 6);
  tangentModulus(1, 1) =
    util::getParam<REAL>("tangentModulus11");
  tangentModulus(1, 2) =
    util::getParam<REAL>("tangentModulus12");
  tangentModulus(1, 3) =
    util::getParam<REAL>("tangentModulus13");
  tangentModulus(2, 1) =
    util::getParam<REAL>("tangentModulus21");
  tangentModulus(2, 2) =
    util::getParam<REAL>("tangentModulus22");
  tangentModulus(2, 3) =
    util::getParam<REAL>("tangentModulus23");
  tangentModulus(3, 1) =
    util::getParam<REAL>("tangentModulus31");
  tangentModulus(3, 2) =
    util::getParam<REAL>("tangentModulus32");
  tangentModulus(3, 3) =
    util::getParam<REAL>("tangentModulus33");
  tangentModulus(4, 4) =
    util::getParam<REAL>("tangentModulus44");
  tangentModulus(5, 5) =
    util::getParam<REAL>("tangentModulus55");
  tangentModulus(6, 6) =
    util::getParam<REAL>("tangentModulus66");
  sigma11 = 0;
  sigma12 = 0;
  sigma13 = 0;
  sigma21 = 0;
  sigma22 = 0;
  sigma23 = 0;
  sigma31 = 0;
  sigma32 = 0;
  sigma33 = 0;
} // end initial()

void
PeriParticle::eraseRecvPeriBonds()
{

  // BB: Feb 3, 2017:
  // Not an efficient operation
  // Better approach may be to use a list if random access of vector
  // members is not needed
  bondVec.erase(std::remove_if(bondVec.begin(), bondVec.end(),
                               [](PeriBondP bond) {
                                 if (bond->getIsRecv() == true) {
                                   return true;
                                 }
                                 return false;
                               }),
                bondVec.end());

  /*
  for (auto bond = bondVec.begin(); bond != bondVec.end();) {
    if ((*bond)->getIsRecv() == true) {
      //		    delete (*bond);
      //		    *bond=NULL;
      bond = bondVec.erase(bond);
    } else
      bond++;
  }
  */
} // eraseRecvPeriBonds()

void
PeriParticle::assignSigma()
{
  sigma11 = sigma(1, 1);
  sigma12 = sigma(1, 2);
  sigma13 = sigma(1, 3);
  sigma21 = sigma(2, 1);
  sigma22 = sigma(2, 2);
  sigma23 = sigma(2, 3);
  sigma31 = sigma(3, 1);
  sigma32 = sigma(3, 2);
  sigma33 = sigma(3, 3);
} // assignSigma

void
PeriParticle::assignKinv()
{
  Kinv11 = Kinv(1, 1);
  Kinv12 = Kinv(1, 2);
  Kinv13 = Kinv(1, 3);
  Kinv21 = Kinv(2, 1);
  Kinv22 = Kinv(2, 2);
  Kinv23 = Kinv(2, 3);
  Kinv31 = Kinv(3, 1);
  Kinv32 = Kinv(3, 2);
  Kinv33 = Kinv(3, 3);
} // assignSigma

} // end pd
