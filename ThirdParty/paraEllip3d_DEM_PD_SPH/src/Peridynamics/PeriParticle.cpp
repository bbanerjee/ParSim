// Function Definitions
#include <Peridynamics/PeriParticle.h>
#include <Core/Util/Utility.h>

namespace pd {

using dem::InputParameter;

PeriParticle::PeriParticle()
{
  d_id = 0;
  d_isAlive = true;
  d_initPosition = dem::Vec(0, 0, 0);
  d_mass = 0.0;
  d_volume = 0.0;
  d_displacement = 0.0;
  d_velocity = 0.0;
  d_velocityHalf = 0.0;
  d_acceleration = 0.0;
  d_horizonSize = -1.0e16;

  d_sigma = dem::zeros(3, 3);
  d_deformationGradient = dem::zeros(3, 3);
  d_deformationGradientHalf = dem::zeros(3, 3);
  d_Kinv = dem::zeros(3, 3);
  d_isv11 = 0;

  // 1---implicit, 2---explicit
  if (util::getParam<int>("typeConstitutive") == 1) { 
    d_isv11 = util::getParam<REAL>("chi");
  } else {
    d_isv11 = util::getParam<REAL>("c");
  }

  d_tangentModulus = dem::zeros(6, 6);
  d_tangentModulus(1, 1) = util::getParam<REAL>("tangentModulus11");
  d_tangentModulus(1, 2) = util::getParam<REAL>("tangentModulus12");
  d_tangentModulus(1, 3) = util::getParam<REAL>("tangentModulus13");
  d_tangentModulus(2, 1) = util::getParam<REAL>("tangentModulus21");
  d_tangentModulus(2, 2) = util::getParam<REAL>("tangentModulus22");
  d_tangentModulus(2, 3) = util::getParam<REAL>("tangentModulus23");
  d_tangentModulus(3, 1) = util::getParam<REAL>("tangentModulus31");
  d_tangentModulus(3, 2) = util::getParam<REAL>("tangentModulus32");
  d_tangentModulus(3, 3) = util::getParam<REAL>("tangentModulus33");
  d_tangentModulus(4, 4) = util::getParam<REAL>("tangentModulus44");
  d_tangentModulus(5, 5) = util::getParam<REAL>("tangentModulus55");
  d_tangentModulus(6, 6) = util::getParam<REAL>("tangentModulus66");

  d_sigma11 = 0;
  d_sigma12 = 0;
  d_sigma13 = 0;
  d_sigma21 = 0;
  d_sigma22 = 0;
  d_sigma23 = 0;
  d_sigma31 = 0;
  d_sigma32 = 0;
  d_sigma33 = 0;

  d_Kinv11 = 0;
  d_Kinv12 = 0;
  d_Kinv13 = 0;
  d_Kinv21 = 0;
  d_Kinv22 = 0;
  d_Kinv23 = 0;
  d_Kinv31 = 0;
  d_Kinv32 = 0;
  d_Kinv33 = 0;
  d_BondedDEMParticleID.clear();
} // end PeriParticle()

PeriParticle::PeriParticle(ParticleID id, REAL x, REAL y, REAL z)
{
  d_id = id;
  d_isAlive = true;
  d_initPosition.setX(x);
  d_initPosition.setY(y);
  d_initPosition.setZ(z);
  d_mass = 0.0;
  d_volume = 0.0;
  d_displacement = 0.0;
  d_velocity = 0.0;
  d_velocityHalf = 0.0;
  d_acceleration = 0.0;
  d_horizonSize = -1.0e16;

  d_sigma = dem::zeros(3, 3);
  d_deformationGradient = dem::zeros(3, 3);
  d_deformationGradientHalf = dem::zeros(3, 3);
  d_Kinv = dem::zeros(3, 3);
  d_isv11 = 0;
  // 1---implicit, 2---explicit
  if (util::getParam<int>("typeConstitutive") == 1) { 
    d_isv11 = util::getParam<REAL>("chi");
  } else {
    d_isv11 = util::getParam<REAL>("c");
  }

  d_tangentModulus = dem::zeros(6, 6);
  d_tangentModulus(1, 1) = util::getParam<REAL>("tangentModulus11");
  d_tangentModulus(1, 2) = util::getParam<REAL>("tangentModulus12");
  d_tangentModulus(1, 3) = util::getParam<REAL>("tangentModulus13");
  d_tangentModulus(2, 1) = util::getParam<REAL>("tangentModulus21");
  d_tangentModulus(2, 2) = util::getParam<REAL>("tangentModulus22");
  d_tangentModulus(2, 3) = util::getParam<REAL>("tangentModulus23");
  d_tangentModulus(3, 1) = util::getParam<REAL>("tangentModulus31");
  d_tangentModulus(3, 2) = util::getParam<REAL>("tangentModulus32");
  d_tangentModulus(3, 3) = util::getParam<REAL>("tangentModulus33");
  d_tangentModulus(4, 4) = util::getParam<REAL>("tangentModulus44");
  d_tangentModulus(5, 5) = util::getParam<REAL>("tangentModulus55");
  d_tangentModulus(6, 6) = util::getParam<REAL>("tangentModulus66");

  d_sigma11 = 0;
  d_sigma12 = 0;
  d_sigma13 = 0;
  d_sigma21 = 0;
  d_sigma22 = 0;
  d_sigma23 = 0;
  d_sigma31 = 0;
  d_sigma32 = 0;
  d_sigma33 = 0;

  d_Kinv11 = 0;
  d_Kinv12 = 0;
  d_Kinv13 = 0;
  d_Kinv21 = 0;
  d_Kinv22 = 0;
  d_Kinv23 = 0;
  d_Kinv31 = 0;
  d_Kinv32 = 0;
  d_Kinv33 = 0;
  d_BondedDEMParticleID.clear();
} // end PeriParticle()

PeriParticle::PeriParticle(const PeriParticle& pt)
{
  d_id = pt.d_id;
  d_isAlive = pt.d_isAlive;
  d_initPosition = pt.d_initPosition;
  d_mass = pt.d_mass;
  d_volume = pt.d_volume;
  d_displacement = pt.d_displacement;
  d_velocity = pt.d_velocity;
  d_velocityHalf = pt.d_velocityHalf;
  d_acceleration = pt.d_acceleration;
  d_sigma = pt.d_sigma;
  d_deformationGradient = pt.d_deformationGradient;
  d_deformationGradientHalf = pt.d_deformationGradientHalf;
  d_Kinv = pt.d_Kinv;
  d_isv11 = pt.d_isv11;
  d_tangentModulus = pt.d_tangentModulus;
  d_horizonSize = pt.d_horizonSize;
  // in order to keep stress values after gathering
  d_sigma11 = pt.d_sigma11;
  d_sigma12 = pt.d_sigma12;
  d_sigma13 = pt.d_sigma13;

  d_sigma21 = pt.d_sigma21;
  d_sigma22 = pt.d_sigma22;
  d_sigma23 = pt.d_sigma23;

  d_sigma31 = pt.d_sigma31;
  d_sigma32 = pt.d_sigma32;
  d_sigma33 = pt.d_sigma33;

  d_Kinv11 = pt.d_Kinv11;
  d_Kinv12 = pt.d_Kinv12;
  d_Kinv13 = pt.d_Kinv13;

  d_Kinv21 = pt.d_Kinv21;
  d_Kinv22 = pt.d_Kinv22;
  d_Kinv23 = pt.d_Kinv23;

  d_Kinv31 = pt.d_Kinv31;
  d_Kinv32 = pt.d_Kinv32;
  d_Kinv33 = pt.d_Kinv33;
  d_BondedDEMParticleID = pt.d_BondedDEMParticleID;
}

PeriParticle::~PeriParticle()
{
  d_bondVec.clear();

  d_sigma.clear();
  d_deformationGradient.clear();
  d_deformationGradientHalf.clear();
  d_Kinv.clear();
  d_tangentModulus.clear();

} // end PeriParticle()

void
PeriParticle::releaseBondVec()
{
  d_bondVec.clear();
} // releaseBondVec

void
PeriParticle::constructMatrixMember()
{
  d_sigma = dem::zeros(3, 3);
  d_deformationGradient = dem::zeros(3, 3);
  d_deformationGradientHalf = dem::zeros(3, 3);
  d_Kinv = dem::zeros(3, 3);
  d_Kinv(1, 1) = d_Kinv11;
  d_Kinv(1, 2) = d_Kinv12;
  d_Kinv(1, 3) = d_Kinv13;
  d_Kinv(2, 1) = d_Kinv21;
  d_Kinv(2, 2) = d_Kinv22;
  d_Kinv(2, 3) = d_Kinv23;
  d_Kinv(3, 1) = d_Kinv31;
  d_Kinv(3, 2) = d_Kinv32;
  d_Kinv(3, 3) = d_Kinv33;
  //	    isv = dem::zeros(1,5);
  //	    if(util::getParam<int>("typeConstitutive") ==
  // 1){
  //// 1---implicit, 2---explicit
  //	    	isv(1,1) = util::getParam<REAL>("chi");
  //	    }
  //	    else{
  //	    	isv(1,1) = util::getParam<REAL>("c");
  //	    }

  d_tangentModulus = dem::zeros(6, 6);
  d_tangentModulus(1, 1) = util::getParam<REAL>("tangentModulus11");
  d_tangentModulus(1, 2) = util::getParam<REAL>("tangentModulus12");
  d_tangentModulus(1, 3) = util::getParam<REAL>("tangentModulus13");
  d_tangentModulus(2, 1) = util::getParam<REAL>("tangentModulus21");
  d_tangentModulus(2, 2) = util::getParam<REAL>("tangentModulus22");
  d_tangentModulus(2, 3) = util::getParam<REAL>("tangentModulus23");
  d_tangentModulus(3, 1) = util::getParam<REAL>("tangentModulus31");
  d_tangentModulus(3, 2) = util::getParam<REAL>("tangentModulus32");
  d_tangentModulus(3, 3) = util::getParam<REAL>("tangentModulus33");
  d_tangentModulus(4, 4) = util::getParam<REAL>("tangentModulus44");
  d_tangentModulus(5, 5) = util::getParam<REAL>("tangentModulus55");
  d_tangentModulus(6, 6) = util::getParam<REAL>("tangentModulus66");
}

void
PeriParticle::replaceHorizonSizeIfLarger(REAL size)
{
  d_horizonSize = (d_horizonSize < size) ? size : d_horizonSize;
} // end replaceHorizonSizeIfLarger()

void
PeriParticle::calcParticleKinv()
{

  dem::Matrix K(3, 3);
  for (auto& bond : d_bondVec) {

    // check which pt1 or pt2 in (*bond) is the center, namely (*pt)
    bool is_pt1 = false; // true when (*pt1) is the center
    if (this == bond->getPt1().get()) {
      is_pt1 = true;
    }

    // dem::Vec xi = (*bond)->getXi(is_pt1);
    // K += dyadicProduct(xi,
    // xi)*(*bond)->volume(is_pt1)*(*bond)->getWeight();
    K = K + bond->getMicroK(is_pt1);

  } // end bond

  //// for numerical purpose, to be deleted later
  // K = 1.0/(horizonSize*horizonSize)*K;
  //
  //// inverse of matrix K
  // Kinv = K.getInvs()/(horizonSize*horizonSize);

  d_Kinv = inv(K);

  /*
  std::ostringstream out;
  out << "ParticleID = " << d_id  << " Kinv = " << Kinv << "\n";
  std::cout << out.str();
  */

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
  for (auto& bond : d_bondVec) {
    if (bond->getIsAlive())
      num_bonds++; // if alive
  }                // end bond

  // disable a particle
  if (num_bonds < 1) { // as rigid particle
    d_isAlive = false;
    //std::cout << "A particle is disabled due to the lack of bond" << std::endl;
  }
} // end checkParticleAlive()

void
PeriParticle::calcParticleStress()
{

  if (!d_isAlive) { // not alive
    d_sigma = dem::zeros(3, 3);
  } else {
    // calculate deformation gradient tensor at current and half step
    dem::Matrix N(3, 3);      // matrix N at n+1 step
    dem::Matrix N_half(3, 3); // matrix N at n+1/2 step
    dem::Matrix N_deltaU(3, 3);
    // matrix N, corresponding to \mathbf{u}^{n+1} - \mathbf{u}^{n},
    // used to calculate \nabla (\mathbf{u}^{n+1} - \mathbf{u}^{n})

    for (auto& bond : d_bondVec) {
      
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

    d_deformationGradient = N * d_Kinv;
    d_deformationGradientHalf = N_half * d_Kinv;

    /*
    std::ostringstream out;
    out << "Particle = " << d_id << " N = " << N << " Kinv = " << Kinv
              << " DefGrad = " << deformationGradient << "\n";
    std::cout << out.str();
    */

    REAL eps = 1.0e-2;
    if (det(d_deformationGradient) < eps || det(d_deformationGradientHalf) < eps) {
      // calculate the determinant of deformationGraident and
      // deformationGradientHalf,
      // if the determinants are too small, then this particle is disabled,
      // isAlive = false
      d_isAlive = false; // disabled particle
      d_sigma = dem::zeros(3, 3);
      //std::cout << "A particle is disabled because det[F] < 0.0" << std::endl;
    } else {
      if (util::getParam<int>("typeConstitutive") == 1) {
        // Linear Elasticity, for testing purpose
        dem::Matrix identity3x3(3, 3);
        identity3x3(1, 1) = 1;
        identity3x3(2, 2) = 1;
        identity3x3(3, 3) = 1;
        dem::Matrix dudx =
          (d_deformationGradient - identity3x3) * (inv(d_deformationGradient));
        dem::Matrix voight_strain(6, 1);
        voight_strain(1, 1) = dudx(1, 1);
        voight_strain(2, 1) = dudx(2, 2);
        voight_strain(3, 1) = dudx(3, 3);
        voight_strain(4, 1) = dudx(2, 3) + dudx(3, 2);
        voight_strain(5, 1) = dudx(1, 3) + dudx(3, 1);
        voight_strain(6, 1) = dudx(1, 2) + dudx(2, 1);
        dem::Matrix voight_sigma = d_tangentModulus * voight_strain;
        d_sigma(1, 1) = voight_sigma(1, 1);
        d_sigma(2, 2) = voight_sigma(2, 1);
        d_sigma(3, 3) = voight_sigma(3, 1);
        d_sigma(2, 3) = voight_sigma(4, 1);
        d_sigma(1, 3) = voight_sigma(5, 1);
        d_sigma(1, 2) = voight_sigma(6, 1);
        d_sigma(2, 1) = d_sigma(1, 2);
        d_sigma(3, 1) = d_sigma(1, 3);
        d_sigma(3, 2) = d_sigma(2, 3);
      } else if (util::getParam<int>("typeConstitutive") ==
                 2) {
        // calculate G, \nabla \Delta \bf{u}
        dem::Matrix G = N_deltaU * d_Kinv * inv(d_deformationGradientHalf);
        dem::Matrix Gsymm = 0.5 * (G + trans(G)); // symmetric part of G
        dem::Matrix Gskew = 0.5 * (G - trans(G)); // skew part of G

        dem::Matrix voight_Gsymm(6, 1);
        voight_Gsymm(1, 1) = Gsymm(1, 1);
        voight_Gsymm(2, 1) = Gsymm(2, 2);
        voight_Gsymm(3, 1) = Gsymm(3, 3);
        voight_Gsymm(4, 1) = Gsymm(2, 3);
        voight_Gsymm(5, 1) = Gsymm(1, 3);
        voight_Gsymm(6, 1) = Gsymm(1, 2);

        dem::Matrix voight_delta_sigma = d_tangentModulus * voight_Gsymm;

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
        dem::Matrix trial_sigma = Q * d_sigma * trans(Q) + delta_sigma;

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
        REAL cn = d_isv11; //?
        REAL f_trial = L2norm_deviatoric_trial_sigma -
                       (Aphi * cn - Bphi * 1.0 / 3.0 * trace_trial_sigma);

        if (f_trial < 0) { // elasticity
          d_sigma = trial_sigma;
          d_isv11 = cn;
        } else { // plasticity

          REAL Bpsi = util::getParam<REAL>("Bpsi");
          REAL Hc = util::getParam<REAL>("hci");
          REAL KBulk = util::getParam<REAL>("kBulk");
          REAL mu = util::getParam<REAL>("mu");
          REAL delta_gamma =
            f_trial / (2.0 * mu + KBulk * Bphi * Bpsi + Hc * Aphi * Aphi);
          d_sigma = trial_sigma -
                  delta_gamma * (KBulk * Bpsi * identity3x3 +
                                 2.0 * mu * deviatoric_trial_sigma /
                                   L2norm_deviatoric_trial_sigma);
          d_isv11 = cn + delta_gamma * Hc * Aphi;
        }
      }

    } // alive particle

  } //

  //std::cout << "ParticleID = " << d_id << " Stress: " << sigma << "\n";

} // end calcParticleStress()

void
PeriParticle::calcParticleAcceleration()
{

  d_acceleration = 0.0;
  dem::Matrix acceleration_matrix(3, 1);
  dem::Matrix xi_ik_matrix(3, 1);
  dem::Matrix PSi;
  dem::Matrix PSk;
  if (d_isAlive) {
    for (auto& bond : d_bondVec) {
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
      PSi = det(pti->d_deformationGradient) * pti->d_sigma *
            inv(trans(pti->d_deformationGradient));
      // Piola Kirchoff stress of particle k
      PSk = det(ptk->d_deformationGradient) * ptk->d_sigma *
            inv(trans(ptk->d_deformationGradient));

      dem::Vec xi_ik = ptk->d_initPosition - pti->d_initPosition;
      xi_ik_matrix(1, 1) = xi_ik.x();
      xi_ik_matrix(2, 1) = xi_ik.y();
      xi_ik_matrix(3, 1) = xi_ik.z();

      acceleration_matrix = acceleration_matrix +
                            bond->getWeight() *
                              (PSi * (pti->d_Kinv) + PSk * (ptk->d_Kinv)) *
                              xi_ik_matrix * ptk->d_volume;

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
    d_acceleration.setX(acceleration_matrix(1, 1));
    d_acceleration.setY(acceleration_matrix(2, 1));
    d_acceleration.setZ(acceleration_matrix(3, 1));

  } // alive particle

  /*
  std::ostringstream out;
  out << "ParticleID = " << d_id << " Acc: " << acceleration << "\n";
  std::cout << out.str();
  */

} // end calcParticleAcceleration()

void
PeriParticle::updateDisplacement()
{
  REAL deltaT = util::getParam<REAL>("timeStep");
  d_velocityHalf = d_velocity + 0.5 * d_acceleration * deltaT;
  d_prevDisp = d_displacement;
  d_displacement += d_velocityHalf * deltaT;

  /*
  std::ostringstream out;
  out << "Disp: P=" << d_id << " delT = " << deltaT 
            << " acc = " << acceleration
            << " v_n = " << velocity << " v_n+1/2 = " << velocityHalf
            << " u_n = " << prevDisp << " u_n+1 = " << displacement << "\n";
  std::cout << out.str();
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
  d_velocity = 2.0 * d_velocityHalf / atf + d_acceleration * deltaT / atf; 
  d_acceleration = 2.0 * (d_velocity - d_velocityHalf) / deltaT;

  /*
  std::ostringstream out;
  out << "Vel: P=" << d_id << " delT = " << deltaT << " atf = " << atf
            << " acc = " << acceleration
            << " v_n+1 = " << velocity << " v_n+1/2 = " << velocityHalf << "\n";
  std::cout << out.str();
  */

} // end updateVelocity()

void
PeriParticle::initial()
{
  d_displacement = 0.0;
  d_velocity = 0.0;
  d_velocityHalf = 0.0;
  d_acceleration = 0.0;
  d_sigma = dem::zeros(3, 3);
  d_deformationGradient = dem::zeros(3, 3);
  d_deformationGradientHalf = dem::zeros(3, 3);
  d_tangentModulus = dem::zeros(6, 6);
  d_isv11 = 0;
  // 1---implicit, 2---explicit
  if (util::getParam<int>("typeConstitutive") == 1) { 
    d_isv11 = util::getParam<REAL>("chi");
  } else {
    d_isv11 = util::getParam<REAL>("c");
  }
  d_tangentModulus = dem::zeros(6, 6);
  d_tangentModulus(1, 1) = util::getParam<REAL>("tangentModulus11");
  d_tangentModulus(1, 2) = util::getParam<REAL>("tangentModulus12");
  d_tangentModulus(1, 3) = util::getParam<REAL>("tangentModulus13");
  d_tangentModulus(2, 1) = util::getParam<REAL>("tangentModulus21");
  d_tangentModulus(2, 2) = util::getParam<REAL>("tangentModulus22");
  d_tangentModulus(2, 3) = util::getParam<REAL>("tangentModulus23");
  d_tangentModulus(3, 1) = util::getParam<REAL>("tangentModulus31");
  d_tangentModulus(3, 2) = util::getParam<REAL>("tangentModulus32");
  d_tangentModulus(3, 3) = util::getParam<REAL>("tangentModulus33");
  d_tangentModulus(4, 4) = util::getParam<REAL>("tangentModulus44");
  d_tangentModulus(5, 5) = util::getParam<REAL>("tangentModulus55");
  d_tangentModulus(6, 6) = util::getParam<REAL>("tangentModulus66");
  d_sigma11 = 0;
  d_sigma12 = 0;
  d_sigma13 = 0;
  d_sigma21 = 0;
  d_sigma22 = 0;
  d_sigma23 = 0;
  d_sigma31 = 0;
  d_sigma32 = 0;
  d_sigma33 = 0;
} // end initial()

void
PeriParticle::eraseRecvPeriBonds()
{

  // BB: Feb 3, 2017:
  // Not an efficient operation
  // Better approach may be to use a list if random access of vector
  // members is not needed
  d_bondVec.erase(std::remove_if(d_bondVec.begin(), d_bondVec.end(),
                               [](PeriBondP bond) {
                                 if (bond->getIsRecv() == true) {
                                   return true;
                                 }
                                 return false;
                               }),
                d_bondVec.end());

} // eraseRecvPeriBonds()

void
PeriParticle::assignSigma()
{
  d_sigma11 = d_sigma(1, 1);
  d_sigma12 = d_sigma(1, 2);
  d_sigma13 = d_sigma(1, 3);
  d_sigma21 = d_sigma(2, 1);
  d_sigma22 = d_sigma(2, 2);
  d_sigma23 = d_sigma(2, 3);
  d_sigma31 = d_sigma(3, 1);
  d_sigma32 = d_sigma(3, 2);
  d_sigma33 = d_sigma(3, 3);
} // assignSigma

void
PeriParticle::assignKinv()
{
  d_Kinv11 = d_Kinv(1, 1);
  d_Kinv12 = d_Kinv(1, 2);
  d_Kinv13 = d_Kinv(1, 3);
  d_Kinv21 = d_Kinv(2, 1);
  d_Kinv22 = d_Kinv(2, 2);
  d_Kinv23 = d_Kinv(2, 3);
  d_Kinv31 = d_Kinv(3, 1);
  d_Kinv32 = d_Kinv(3, 2);
  d_Kinv33 = d_Kinv(3, 3);
} // assignSigma

} // end pd
