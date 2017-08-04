#include "Particle.h"
#include "const.h"
#include "ran.h"
#include "root6.h"
#include <iostream>

//#define MOMENT
#ifdef MOMENT
const std::size_t START = 10000;  // at which time step to apply moment? for moment rotation test only.
#define SLIP  // if defined, stick and slip; otherwise slide.
#endif
// for moment case: timeStep = 5.0e-07; totalSteps = 12000

//#define MINDLIN_ASSUMED

namespace dem {

  Particle::Particle()
    :id(0), type(0), a(0), b(0), c(0),
     young(0), poisson(0),
     currPos(0), prevPos(0),
     currDirecA(0), currDirecB(0), currDirecC(0),
     prevDirecA(0), prevDirecB(0), prevDirecC(0),
     currVeloc(0), prevVeloc(0),
     currOmga(0), prevOmga(0),
     force(0), prevForce(0),
     moment(0), prevMoment(0),
     constForce(0), constMoment(0),
     density(0), mass(0), volume(0), momentJ(0),
     kinetEnergy(0), contactNum(0), inContact(false) {
    for (std::size_t i = 0; i < 10; ++i)
      coef[i] = 0;
  }
  

  void Particle::init () {
    // generate orientation of axle a/b/c using Euler angles

    REAL angle1, angle2, angle3; // angle1=[0,Pi], angle2=[0,2*Pi), angle3=[0,2*Pi]
    angle1 = ran(&idum)*Pi;
    angle2 = ran(&idum)*Pi*2;
    angle3 = ran(&idum)*Pi*2;

    REAL c1, c2, c3, s1, s2, s3;
    c1 = cos(angle1);
    c2 = cos(angle2);
    c3 = cos(angle3);
    s1 = sin(angle1);
    s2 = sin(angle2);
    s3 = sin(angle3);

    REAL l1, m1, n1, l2, m2, n2, l3, m3, n3;
    l1 = c2*c3-c1*s2*s3;  m1 = s2*c3+c1*c2*s3;  n1 = s1*s3;
    l2 = -c2*s3-c1*s2*c3; m2 = -s2*s3+c1*c2*c3; n2 = s1*c3;
    l3 = s1*s2;           m3 = -s1*c2;          n3 = c1;

    currDirecA = Vec(0.000000000000e+00, 1.570796326794e+00, 1.570796326794e+00);
    currDirecB = Vec(1.570796326794e+00, 0.000000000000e+00, 1.570796326794e+00);
    currDirecC = Vec(1.570796326794e+00, 1.570796326794e+00, 0.000000000000e+00);
    //currDirecC = vacos(normalize(vcos(currDirecA) * vcos(currDirecB)));

    prevPos = currPos;
    prevDirecA = currDirecA;
    prevDirecB = currDirecB;
    prevDirecC = currDirecC;
    prevVeloc = currVeloc = 0;
    prevOmga = currOmga = 0;
    force = prevForce = 0;
    moment = prevMoment = 0;
    constForce = constMoment = 0;
    density = dem::Parameter::getSingleton().parameter["specificG"] * 1.0e+3;
    volume = 4.0/3.0*Pi*a*b*c;
    mass = density*volume;
    momentJ = Vec(mass/5*(b*b+c*c), mass/5*(a*a+c*c), mass/5*(a*a+b*b));
    contactNum = 0;
    inContact = false;
    globalCoef();
  }


  Particle::Particle(std::size_t n, std::size_t tp, Vec center, REAL r, REAL yng, REAL poi)
    :id(n), type(tp), a(r), b(r), c(r), young(yng), poisson(poi), currPos(center) {
    init();
  }


  Particle::Particle(std::size_t n, std::size_t tp, Vec center, REAL ra, REAL rb, REAL rc, REAL yng, REAL poi)
    :id(n), type(tp), a(ra), b(rb), c(rc), young(yng), poisson(poi), currPos(center) {
    init();
  }


  Particle::Particle(std::size_t n, std::size_t tp, Vec center, Gradation &grad, REAL yng, REAL poi)
    :id(n), type(tp), young(yng), poisson(poi), currPos(center)  {
    // generate particle size in terms of gradation distribution
    REAL sievenum = grad.getSieveNum();
    REAL randnum = ran(&idum);
    for (std::size_t k = 0; k < sievenum; ++k) {
      if ( randnum <= grad.getPercent()[sievenum-1-k]) {
	a = grad.getSize()[sievenum-1-k]; // use a for sieving (where a >= b >= c) 
	break;
      }
    }
  
#ifdef RANDOM_SHAPE
    grad.setPtclRatioBA(ran(&idum));
    grad.setPtclRatioCA(ran(&idum));
#endif
  
    b = 1;
    c = a * grad.getPtclRatioCA();
  
    init();
  }
  

  Particle::Particle(std::size_t n, std::size_t tp, Vec dim, Vec position, Vec dirca, Vec dircb, Vec dircc, REAL yng, REAL poi)
    :id(n), type(tp), young(yng), poisson(poi) {
    a = dim.getX();
    b = dim.getY();
    c = dim.getZ();
    currPos = prevPos = position;
    currDirecA = prevDirecA = dirca;
    currDirecB = prevDirecB = dircb;
    currDirecC = prevDirecC = dircc;
    currVeloc = prevVeloc = 0;
    currOmga = prevOmga = 0;
    force = prevForce = 0;
    moment = prevMoment = 0;
    constForce = constMoment = 0;
    contactNum = 0;
    density = dem::Parameter::getSingleton().parameter["specificG"] * 1.0e3;
    volume = 4.0/3.0*Pi*a*b*c;
    mass = density*volume;
    momentJ = Vec(mass/5*(b*b + c*c), mass/5*(a*a + c*c), mass/5*(a*a + b*b));
    inContact = false;
    globalCoef();
  }


  Vec Particle::globalToLocal(Vec input) const {
    Vec lmn, local;
    lmn = vcos(getCurrDirecA()); local.setX(lmn * input); // l1,m1,n1
    lmn = vcos(getCurrDirecB()); local.setY(lmn * input); // l2,m2,n2
    lmn = vcos(getCurrDirecC()); local.setZ(lmn * input); // l3,m3,n3
    return local;
  }


  Vec Particle::localToGlobal(Vec input) const {
    Vec lmn, global;
    lmn = vcos( Vec(currDirecA.getX(),currDirecB.getX(),currDirecC.getX()) ); global.setX(lmn * input); // l1,l2,l3
    lmn = vcos( Vec(currDirecA.getY(),currDirecB.getY(),currDirecC.getY()) ); global.setY(lmn * input); // m1,m2,n3
    lmn = vcos( Vec(currDirecA.getZ(),currDirecB.getZ(),currDirecC.getZ()) ); global.setZ(lmn * input); // n1,n2,n3
    return global;
  }

  
  // 1: rotational energy is 1/2(I1*w1^2+I2*w2^2+I3*w3^2), where each term is expressed in local frame.
  // 2. angular velocities in global frame needs to be converted to those in local frame.
  REAL Particle::getTransEnergy() const {
    return mass*pow(vfabs(currVeloc),2)/2;
  }
  
  
  REAL Particle::getRotatEnergy() const {
    Vec currLocalOmga = globalToLocal(currOmga);
  
    return momentJ.getX()*pow(currLocalOmga.getX(),2)/2 
      + momentJ.getY()*pow(currLocalOmga.getY(),2)/2
      + momentJ.getZ()*pow(currLocalOmga.getZ(),2)/2;
  }
  
  
  REAL Particle::getKinetEnergy() const {
    return getTransEnergy() + getRotatEnergy();
  }


  REAL Particle::getPotenEnergy(REAL ref) const {
    return dem::Parameter::getSingleton().parameter["gravAccel"]*mass*(currPos.getZ() - ref);
  }

  
  void Particle::getGlobalCoef(REAL coef[]) const {
    for (std::size_t i = 0; i < 10; ++i)
      coef[i] = this->coef[i];
  }


  REAL Particle::surfaceError(Vec pt) const {
    REAL x = pt.getX();
    REAL y = pt.getY();
    REAL z = pt.getZ();
    return coef[0]*x*x + coef[1]*y*y + coef[2]*z*z + coef[3]*x*y + coef[4]*y*z 
      + coef[5]*z*x + coef[6]*x + coef[7]*y + coef[8]*z + coef[9];
  }


  void Particle::globalCoef() {
    // coef[0]-x^2, coef[1]-y^2, coef[2]-z^2, coef[3]-xy, coef[4]-yz, coef[5]-zx
    // coef[6]-x, coef[7]-y, coef[8]-z, coef[9]-const
    if(a == b && b == c){
      coef[0] = 1;
      coef[1] = 1;
      coef[2] = 1;
      coef[3] = 0;
      coef[4] = 0;
      coef[5] = 0;
      coef[6] = -2*currPos.getX();
      coef[7] = -2*currPos.getY();
      coef[8] = -2*currPos.getZ();
      coef[9] = pow(vfabs(currPos),2)-a*a;
      return;
    }
    Vec  v1 = vcos(currDirecA);
    Vec  v2 = vcos(currDirecB);
    Vec  v3 = vcos(currDirecC);
    REAL X0 = currPos.getX();
    REAL Y0 = currPos.getY();
    REAL Z0 = currPos.getZ();
    REAL l1 = v1.getX();
    REAL m1 = v1.getY();
    REAL n1 = v1.getZ();
    REAL l2 = v2.getX();
    REAL m2 = v2.getY();
    REAL n2 = v2.getZ();
    REAL l3 = v3.getX();
    REAL m3 = v3.getY();
    REAL n3 = v3.getZ();
    coef[0] = l1*l1/a/a + l2*l2/b/b + l3*l3/c/c;
    coef[1] = m1*m1/a/a + m2*m2/b/b + m3*m3/c/c;
    coef[2] = n1*n1/a/a + n2*n2/b/b + n3*n3/c/c;
    coef[3] = (2*l1*m1)/a/a + (2*l2*m2)/b/b + (2*l3*m3)/c/c;
    coef[4] = (2*m1*n1)/a/a + (2*m2*n2)/b/b + (2*m3*n3)/c/c;
    coef[5] = (2*l1*n1)/a/a + (2*l2*n2)/b/b + (2*l3*n3)/c/c;
    coef[6] = 
      -2*l1*m1*Y0*pow(a,-2) - 2*l1*n1*Z0*pow(a,-2) - 
      2*l2*m2*Y0*pow(b,-2) - 2*l2*n2*Z0*pow(b,-2) - 
      2*l3*m3*Y0*pow(c,-2) - 2*l3*n3*Z0*pow(c,-2) - 
      2*X0*pow(a,-2)*pow(l1,2) - 2*X0*pow(b,-2)*pow(l2,2) - 
      2*X0*pow(c,-2)*pow(l3,2);
    coef[7] = 
      (-2*l1*m1*X0)/a/a - (2*l2*m2*X0)/b/b - 
      (2*l3*m3*X0)/c/c - (2*m1*m1*Y0)/a/a - 
      (2*m2*m2*Y0)/b/b - (2*m3*m3*Y0)/c/c - 
      (2*m1*n1*Z0)/a/a - (2*m2*n2*Z0)/b/b - 
      (2*m3*n3*Z0)/c/c;
    coef[8] = 
      (-2*l1*n1*X0)/a/a - (2*l2*n2*X0)/b/b - 
      (2*l3*n3*X0)/c/c - (2*m1*n1*Y0)/a/a - 
      (2*m2*n2*Y0)/b/b - (2*m3*n3*Y0)/c/c - 
      (2*n1*n1*Z0)/a/a - (2*n2*n2*Z0)/b/b - 
      (2*n3*n3*Z0)/c/c;
    coef[9] = 
      -1 + 2*l1*m1*X0*Y0*pow(a,-2) + 2*l1*n1*X0*Z0*pow(a,-2) + 
      2*m1*n1*Y0*Z0*pow(a,-2) + 2*l2*m2*X0*Y0*pow(b,-2) + 
      2*l2*n2*X0*Z0*pow(b,-2) + 2*m2*n2*Y0*Z0*pow(b,-2) + 
      2*l3*m3*X0*Y0*pow(c,-2) + 2*l3*n3*X0*Z0*pow(c,-2) + 
      2*m3*n3*Y0*Z0*pow(c,-2) + 
      pow(a,-2)*pow(l1,2)*pow(X0,2) + 
      pow(b,-2)*pow(l2,2)*pow(X0,2) + 
      pow(c,-2)*pow(l3,2)*pow(X0,2) + 
      pow(a,-2)*pow(m1,2)*pow(Y0,2) + 
      pow(b,-2)*pow(m2,2)*pow(Y0,2) + 
      pow(c,-2)*pow(m3,2)*pow(Y0,2) + 
      pow(a,-2)*pow(n1,2)*pow(Z0,2) + 
      pow(b,-2)*pow(n2,2)*pow(Z0,2) + 
      pow(c,-2)*pow(n3,2)*pow(Z0,2);
    REAL divd = coef[0];
    for (std::size_t i = 0; i < 10; ++i)  // when a particle is initialized or updated, coef[0] is set as 1.0
      coef[i] /= divd;
  }
  
  
  bool Particle::intersectWithLine(Vec v, Vec dirc, Vec rt[]) const {
    REAL x0 = v.getX();
    REAL y0 = v.getY();
    REAL z0 = v.getZ();
    REAL p = dirc.getX();
    REAL q = dirc.getY();
    REAL r = dirc.getZ();
    REAL a = coef[0];
    REAL b = coef[1];
    REAL c = coef[2];
    REAL d = coef[3];
    REAL e = coef[4];
    REAL f = coef[5];
    REAL g = coef[6];
    REAL h = coef[7];
    REAL i = coef[8];
    REAL j = coef[9];
  
    REAL A = a*p*p + b*q*q + c*r*r + d*p*q + e*q*r + f*r*p;
    REAL B = 2*a*p*x0 + 2*b*q*y0 + 2*c*r*z0
      + d*p*y0 + d*q*x0 + e*q*z0 + e*r*y0 + f*p*z0 + f*r*x0
      + g*p + h*q + i*r;
    REAL C = a*x0*x0 + b*y0*y0 + c*z0*z0 + d*x0*y0 + e*y0*z0 +f*z0*x0
      + g*x0 + h*y0 + i*z0 + j;
  
    REAL delta = B*B - 4*A*C;
    if (delta < 0) {
      debugInf << "Particle.cpp: iter=" << iteration
		<< " delta < 0 in intersectWithLine()" << std::endl;
      return false;
    }
    else{
      REAL t1 = (-B + sqrt(delta)) / (2*A);
      REAL t2 = (-B - sqrt(delta)) / (2*A);
    
      rt[0].setX(t1*p + x0);
      rt[0].setY(t1*q + y0);
      rt[0].setZ(t1*r + z0);
      rt[1].setX(t2*p + x0);
      rt[1].setY(t2*q + y0);
      rt[1].setZ(t2*r + z0);   
      return true;
    }
  }
  
  
  // 1. This member function is coded based on Mathematical equations in local frame, 
  //    x^2/a^2 + y^2/b^2 + z^2/c^2 = 1, seeking appropriate osculating circle among 
  //    an infinite number of  osculating circles passing through the contact point.
  // 2. r = 2*r1*r2/(r1+r2)
  // 3. It is important to eliminate float exceptions in computations, that is, 
  //    when dz/dx == infinite, coordinate x & z are switched to use dx/dz == 0.
  // 4. When a point is close to the equator, for example, fabs(z) == 0,
  //    float exception is prone to occurring, then a switch is needed as above.
  REAL Particle::getRadius(Vec v) const {
    if(a == b && b == c)
      return a;
  
    REAL per = 1.0e-4; // define when a point is close to equator
    REAL ra = a;       // semi-axles of ellipsoid
    REAL rb = b;
    REAL rc = c;
  
    // get the local coodinates of vector v, the point on the particle's surface
    Vec  v1 = vcos(currDirecA);
    Vec  v2 = vcos(currDirecB);
    Vec  v3 = vcos(currDirecC);
    REAL X0 = currPos.getX();
    REAL Y0 = currPos.getY();
    REAL Z0 = currPos.getZ();
    REAL x1 = v.getX() - X0;
    REAL y1 = v.getY() - Y0;
    REAL z1 = v.getZ() - Z0;
    REAL l1 = v1.getX();
    REAL m1 = v1.getY();
    REAL n1 = v1.getZ();
    REAL l2 = v2.getX();
    REAL m2 = v2.getY();
    REAL n2 = v2.getZ();
    REAL l3 = v3.getX();
    REAL m3 = v3.getY();
    REAL n3 = v3.getZ();
    REAL x = l1*x1 + m1*y1 + n1*z1;
    REAL y = l2*x1 + m2*y1 + n2*z1;
    REAL z = l3*x1 + m3*y1 + n3*z1;
  
    REAL tmp;
    if (fabs(z) <= c*per) {   // switch x & z, use 0 instead of infinity
      tmp = ra; ra = rc; rc = tmp;
      tmp = x; x = z; z = tmp; 
      if (fabs(z) <= a*per) { // switch y & z, use 0 instead of infinity
	tmp = ra; ra = rb; rb = tmp;
	tmp = y; y = z; z = tmp; 
      }     
    }
  
    REAL p = -rc*rc/ra/ra*x/z;
    REAL q = -rc*rc/rb/rb*y/z;
    REAL r = -rc*rc/ra/ra*(1/z + rc*rc/ra/ra*x*x/pow(z,3));
    REAL t = -rc*rc/rb/rb*(1/z + rc*rc/rb/rb*y*y/pow(z,3));
    REAL s = -pow(rc,4)/ra/ra/rb/rb*x*y/pow(z,3);
    REAL n = sqrt(1 + p*p + q*q);
  
    REAL A,B,C;
    A = r*t - s*s;
    B = n*(2*p*q*s - (1 + p*p)*t - (1 + q*q)*r);
    C = n*n*n*n;
  
    // if delta < 0, then it is usually -1.0e-20, caused by computational precision.
    /*
      if (B*B-4*A*C < 0){
      debugInf<< "Particle.cpp: iter=" << iteration
      << " delta < 0 in getRadius()"
      << " delta=" << B*B-4*A*C
      << " -C/B=" << -C/B
      << std::endl;
      }
    */
    return fabs(-C/B*2.0); // 2*r1*r2/(r1 + r2)
  }
  
  
  void Particle::clearContactForce() {

    REAL gravAccel = dem::Parameter::getSingleton().parameter["gravAccel"];
    REAL gravScale = dem::Parameter::getSingleton().parameter["gravScale"];

    force = constForce;
    moment = constMoment;
  
    force += Vec(0, 0, -gravAccel * mass * gravScale); // Unit is Newton, gravScale is for amplification.
    inContact = false;
  
    if (getType() == 3) // ellipsoidal pile
      force -= Vec(0, 0, -gravAccel * mass * gravScale); 
  
#ifdef MOMENT
    REAL m[20] = { 1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
		   80, 70, 60, 50, 40, 30, 20, 10, 0};
#ifdef SLIP
    for (std::size_t i = 0;i < 20; ++i) m[i] *= 1.0e-8;
#else
    for (std::size_t i = 0;i < 20; ++i) m[i] *= 2.0e-8; 
#endif
    std::size_t s[20];
    for (std::size_t i = 0;i < 20; ++i)
      s[i] = START + i*100;
  
    for (std::size_t i = 0;i < 19; ++i)
      if (iteration >= s[i] && iteration < s[i+1] )
	moment += Vec(0,m[i],0);
    if (iteration >= s[19] )
      moment += Vec(0,m[19],0);
#endif
  }
  
  
  void Particle::dragForce() {
    REAL Cd  = dem::Parameter::getSingleton().parameter["Cd"];
    REAL rho = dem::Parameter::getSingleton().parameter["fluidDensity"];

    REAL ux = currVeloc.getX(); 
    REAL uy = currVeloc.getY(); 
    REAL uz = currVeloc.getZ();
    Vec globalDelta = Vec(fabs(ux)*ux, fabs(uy)*uy, fabs(uz)*uz);
    Vec localDelta  = globalToLocal(globalDelta);
    Vec localForce, globalForce;
    // localDelta needs to project in local frame in order to calculate local drag forces
    localForce.setX(-0.5*rho*localDelta.getX()*Cd*Pi*b*c);
    localForce.setY(-0.5*rho*localDelta.getY()*Cd*Pi*c*a);
    localForce.setZ(-0.5*rho*localDelta.getZ()*Cd*Pi*a*b);
    globalForce = localToGlobal(localForce);
    addForce(globalForce);
  }


  // central difference integration method
  void Particle::update() {

    REAL forceDamp = dem::Parameter::getSingleton().parameter["forceDamp"];
    REAL momentDamp = dem::Parameter::getSingleton().parameter["momentDamp"];
    REAL massScale = dem::Parameter::getSingleton().parameter["massScale"];
    REAL mntScale = dem::Parameter::getSingleton().parameter["mntScale"];
    REAL pileRate = dem::Parameter::getSingleton().parameter["pileRate"];
  
    if (getType() == 0 || getType() == 5) { // 0-free, 1-fixed, 5-free bounary particle
      // It is important to distinguish global frame from local frame!
      Vec prevLocalOmga;
      Vec currLocalOmga;
      Vec localMoment;
      REAL atf = forceDamp*2; 
      REAL atm = momentDamp*2; 
    
// do not allow move in y direction, only rotate along y axis
force.setY(0);


      // force: translational kinetics equations are in global frame
      currVeloc = prevVeloc * (2-atf) / (2+atf) + force / (mass * massScale) * timeStep * 2 / (2+atf);
      currPos = prevPos + currVeloc * timeStep;
    
      // moment: angular kinetics (rotational) equations are in local frame,
      // so global values need to be converted to those in local frame when applying equations
      localMoment = globalToLocal(moment);
      prevLocalOmga = globalToLocal(prevOmga);
    
      currLocalOmga.setX( prevLocalOmga.getX() * (2-atm) / (2+atm) + localMoment.getX() / (momentJ.getX() * mntScale) * timeStep * 2 / (2+atm) ); 
      currLocalOmga.setY( prevLocalOmga.getY() * (2-atm) / (2+atm) + localMoment.getY() / (momentJ.getY() * mntScale) * timeStep * 2 / (2+atm) );
      currLocalOmga.setZ( prevLocalOmga.getZ() * (2-atm) / (2+atm) + localMoment.getZ() / (momentJ.getZ() * mntScale) * timeStep * 2 / (2+atm) );
    
      // convert local angular velocities to those in global frame in order to rotate a particle in global space
      currOmga = localToGlobal(currLocalOmga);
    
currOmga.setX(0);
currOmga.setZ(0);

      currDirecA = vacos(normalize(rotateVec(vcos(prevDirecA),currOmga * timeStep)));
      currDirecB = vacos(normalize(rotateVec(vcos(prevDirecB),currOmga * timeStep)));
      currDirecC = vacos(normalize(rotateVec(vcos(prevDirecC),currOmga * timeStep)));
    }
#ifdef MOMENT
    else if (getType() == 2) { //special case 2 (moment): translate first, then rotate
      Vec prevLocalOmga;
      Vec currLocalOmga;
      Vec localMoment;
      REAL atf = forceDamp*2; 
      REAL atm = momentDamp*2; 
      currVeloc = prevVeloc * (2-atf) / (2+atf) + force / (mass * massScale) * timeStep * 2 / (2+atf);
      if (iteration < START)
	currPos = prevPos + currVeloc * timeStep;
	
      localMoment = globalToLocal(moment);
      prevLocalOmga = globalToLocal(prevOmga);
    
      currLocalOmga.setX( prevLocalOmga.getX() * (2-atm) / (2+atm) + localMoment.getX() / (momentJ.getX() * mntScale) * timeStep * 2 / (2+atm) ); 
      currLocalOmga.setY( prevLocalOmga.getY() * (2-atm) / (2+atm) + localMoment.getY() / (momentJ.getY() * mntScale) * timeStep * 2 / (2+atm) );
      currLocalOmga.setZ( prevLocalOmga.getZ() * (2-atm) / (2+atm) + localMoment.getZ() / (momentJ.getZ() * mntScale) * timeStep * 2 / (2+atm) );
    
      if (iteration >= START) {	
	currOmga = localToGlobal(currLocalOmga);
      
	currDirecA = vacos(normalize(rotateVec(vcos(prevDirecA),currOmga * timeStep)));
	currDirecB = vacos(normalize(rotateVec(vcos(prevDirecB),currOmga * timeStep)));
	currDirecC = vacos(normalize(rotateVec(vcos(prevDirecC),currOmga * timeStep)));
      }
    }
#endif
    else if (getType() == 3) { // special case 3 (displacemental ellipsoidal pile): translate in vertical direction only
      currVeloc.setX(0);	
      currVeloc.setY(0);
      currVeloc.setZ(-pileRate);
      currPos = prevPos + currVeloc * timeStep;
    }
    else if (getType() == 4) { // special case 4 (impacting ellipsoidal penetrator): impact with inital velocity in vertical direction only 
      REAL atf = forceDamp*2; 
      currVeloc = prevVeloc * (2-atf) / (2+atf) + force / (mass * massScale) * timeStep * 2 / (2+atf);
      currVeloc.setX(0);	
      currVeloc.setY(0);
      currPos = prevPos + currVeloc * timeStep;
    }
    else if (getType() == 6) { // translation only, no rotation
      REAL atf = forceDamp*2; 
      currVeloc = prevVeloc * (2-atf) / (2+atf) + force / (mass * massScale) * timeStep * 2 / (2+atf);
      currPos = prevPos + currVeloc * timeStep;
    }
  
    // Below is needed for all cases
    // ensure three axles perpendicular to each other, and being unit vector
    if(currDirecA == 0)
      currDirecA = vacos(normalize(vcos(currDirecB) % vcos(currDirecC)));
    if(currDirecB == 0)
      currDirecB = vacos(normalize(vcos(currDirecC) % vcos(currDirecA)));
    if(currDirecC == 0)
      currDirecC = vacos(normalize(vcos(currDirecA) % vcos(currDirecB)));
  
    prevPos = currPos;
    prevDirecA = currDirecA;
    prevDirecB = currDirecB;
    prevDirecC = currDirecC;
    prevVeloc = currVeloc;
    prevOmga = currOmga;
    prevForce = force; 
    prevMoment = moment;
  
    contactNum = 0;
    globalCoef(); // every time the particle is updated, the algebra expression is also updated
  }
  

  bool Particle::nearestPTOnPlane(REAL p, REAL q, REAL r, REAL s, Vec &ptnp) const {
    if(a == b && b == c) {
      Vec tnm = Vec(p,q,r) / sqrt(p * p + q * q + r * r);
      // signed distance from particle center to plane
      REAL l_nm = (currPos.getX() * p + currPos.getY() * q + currPos.getZ() * r + s) / sqrt(p * p + q * q + r * r); 
      ptnp = currPos - l_nm * tnm;
      if( (a - fabs(l_nm)) / (2.0*a) > dem::Parameter::getSingleton().parameter["minRelaOverlap"]) // intersect
	return true;
      else  // no intersect
	return false;
    }
  
    REAL a = coef[0];
    REAL b = coef[1];
    REAL c = coef[2];
    REAL d = coef[3];
    REAL e = coef[4];
    REAL f = coef[5];
    REAL g = coef[6];
    REAL h = coef[7];
    REAL i = coef[8];
    REAL j = coef[9];

    REAL domi = 
      e*e*p*p + 4*c*d*p*q - 4*a*c*q*q + 
      f*f*q*q - 2*d*f*q*r + 
      d*d*r*r - 
      2*e*(f*p*q + d*p*r - 2*a*q*r) - 
      4*b*(c*p*p + r*(-f*p + a*r));

    REAL x = 
      (-(f*i*q*q) - 2*b*i*p*r + f*h*q*r + d*i*q*r + 
       2*b*g*r*r - d*h*r*r - e*e*p*s - 
       2*b*f*r*s - 2*c*(h*p*q - g*q*q - 2*b*p*s + 
			d*q*s) + e*(i*p*q + h*p*r - 2*g*q*r + f*q*s + 
				    d*r*s))/domi;
    REAL y = 
      (f*i*p*q - 2*f*h*p*r + d*i*p*r + f*g*q*r - 
       2*a*i*q*r - d*g*r*r + 2*a*h*r*r - 
       f*f*q*s + d*f*r*s + 
       2*c*(h*p*p - g*p*q - d*p*s + 2*a*q*s) + 
       e*(-i*p*p + g*p*r + f*p*s - 2*a*r*s))/domi;

    REAL z = 
      (f*h*p*q - 2*d*i*p*q - f*g*q*q + 
       2*a*i*q*q + d*h*p*r + d*g*q*r - 2*a*h*q*r + 
       d*f*q*s - d*d*r*s + 
       e*(-h*p*p + g*p*q + d*p*s - 2*a*q*s) + 
       2*b*(i*p*p - g*p*r - f*p*s + 2*a*r*s))/domi;

    ptnp = Vec(x, y, z);
  
    REAL val = a*x*x + b*y*y + c*z*z + d*x*y + e*y*z + f*x*z + g*x + h*y + i*z + j;
  
    if (val >= 0) // not intersect
      return false;
    else // intersect
      return true;
  }
  

  void Particle::planeRBForce(planeBoundary *plane,
			      std::map<std::size_t, std::vector<BoundaryTgt> > &BdryTgtMap,
			      std::vector<BoundaryTgt> &vtmp) {
    // (p, q, r) are in the same direction as the outward normal vector,
    // hence it is not necessary to provide information about which side the particle is about the plane.
    REAL p,q,r,s;
    Vec dirc = normalize(plane->getDirec());
    p = dirc.getX();
    q = dirc.getY();
    r = dirc.getZ();
    s = -dirc * plane->getPoint(); // plane equation: p(x-x0) + q(y-y0) + r(z-z0) = 0, that is, px + qy + rz + s = 0

    Vec pt1;
    if (!nearestPTOnPlane(p, q, r, s, pt1)) // the particle and the plane does not intersect
      return;
  
    // if particle and plane intersect:
    ++contactNum;
    inContact = true;
    Vec rt[2];
    if (!intersectWithLine(pt1, dirc, rt)) // the line and ellipsoid surface does not intersect
      return;
  
    Vec pt2;
    ///* universal, allow for large overlap
    if (p*rt[0].getX()+q*rt[0].getY()+r*rt[0].getZ()+s > 0)
      pt2 = rt[0];
    else
      pt2 = rt[1];
    //*/
    /* not universal, only allow for small overlap
    if (vfabs(rt[0]-pt1) < vfabs(rt[1]-pt1) )
      pt2 = rt[0];
    else
      pt2 = rt[1];
    */
  
    // obtain normal force
    REAL penetr = vfabs(pt1 - pt2);
    if (penetr / (2.0*getRadius(pt2) ) <= dem::Parameter::getSingleton().parameter["minRelaOverlap"])
      return;
  
    REAL R0 = getRadius(pt2);
    REAL E0 = young/(1-poisson*poisson); // rigid wall has infinite young's modulus
    REAL allowedOverlap = 2.0 * R0 * dem::Parameter::getSingleton().parameter["maxRelaOverlap"];
    if (penetr > allowedOverlap) {
#ifndef NDEBUG
      std::stringstream inf;
      inf.setf(std::ios::scientific, std::ios::floatfield);
      inf << "Particle.cpp: iter=" << std::setw(8) << iteration 
	  << "  ptcl=" << std::setw(8) << getId()
	  << "  bdry=" << std::setw(8) << plane->getId()
	  << " penetr=" << std::setw(OWID) << penetr 
	  << " allow="  << std::setw(OWID) << allowedOverlap
	  << std::endl;
      MPI_Status status;
      int length = OWID*2 + 8*3 + 19 + 7*3 + 8 + 1;
      MPI_File_write_shared(overlapInf, const_cast<char*> (inf.str().c_str()), length, MPI_CHAR, &status);
#endif
      penetr = allowedOverlap;
    }
  
    REAL measureOverlap = dem::Parameter::getSingleton().parameter["measureOverlap"];  
    penetr = nearbyint (penetr/measureOverlap) * measureOverlap;
    REAL contactRadius = sqrt(penetr*R0);
    Vec normalDirc = -dirc;
    // pow(penetr,1.5), a serious bug
    Vec normalForce = sqrt(penetr * penetr * penetr) * sqrt(R0) * 4 * E0/3 * normalDirc;
  
    /*
      debugInf << ' ' << iteration
      << ' ' << getId()
      << ' ' << plane->boundaryId
      << ' ' << pt1.getX()
      << ' ' << pt1.getY()
      << ' ' << pt1.getZ()
      << ' ' << rt[0].getX()
      << ' ' << rt[0].getY()
      << ' ' << rt[0].getZ()
      << ' ' << rt[1].getX()
      << ' ' << rt[1].getY()
      << ' ' << rt[1].getZ()
      << ' ' << vfabs(rt[0]-pt1)
      << ' ' << vfabs(rt[1]-pt1)
      << ' ' << penetr
      << std::endl;
    */
  
    // apply normal force
    addForce(normalForce);
    addMoment(((pt1 + pt2)/2 - currPos) % normalForce);
  
    // obtain normal damping force
    Vec veloc2 = getCurrVeloc() + getCurrOmga() % ((pt1 + pt2)/2 - getCurrPos());
    REAL kn = pow(6 * vfabs(normalForce) * R0 * pow(E0,2), 1.0/3.0);
    REAL dampCritical = 2 * sqrt(getMass() * kn); // critical damping
    Vec cntDampingForce = dem::Parameter::getSingleton().parameter["contactDamp"] * dampCritical * ((-veloc2) * normalDirc) * normalDirc;
  
    // apply normal damping force
    addForce(cntDampingForce);
    addMoment(((pt1 + pt2)/2 - currPos) % cntDampingForce);
  
    Vec tgtForce = 0;
    if (dem::Parameter::getSingleton().parameter["boundaryFric"] != 0) {
      // checkin previous tangential force and displacement
      Vec  prevTgtForce;
      Vec  prevTgtDisp;
      bool prevTgtLoading = true;
      Vec  tgtDispStart;
      REAL tgtPeak = 0;
    
      bool tgtLoading = true;
      std::vector<BoundaryTgt>::iterator it;
      for (it = BdryTgtMap[plane->getId()].begin(); it != BdryTgtMap[plane->getId()].end(); ++it){
	if (id == it->particleId) {
	  prevTgtForce   = it->tgtForce;
	  prevTgtDisp    = it->tgtDisp;
	  prevTgtLoading = it->tgtLoading;
	  tgtDispStart   = it->tgtDispStart;
	  tgtPeak        = it->tgtPeak;
	  break;
	}
      }
    
      // obtain tangtential force
      REAL G0 = young/2/(1+poisson);
      // Vr = Vb + w (crossdot) r, each item needs to be in either global or local frame; 
      //      here global frame is used for better convenience.
      Vec relaDispInc = (currVeloc + currOmga % ((pt1+pt2)/2 - currPos)) * timeStep;
      Vec tgtDispInc  = relaDispInc - (relaDispInc * normalDirc) * normalDirc;
      Vec tgtDisp     = prevTgtDisp + tgtDispInc; // prevTgtDisp read by checkin
      Vec TgtDirc;
    
      if (vfabs(tgtDisp) == 0)
	TgtDirc = 0;
      else
	TgtDirc = normalize(-tgtDisp); // TgtDirc points along tangential forces exerted on particle 1
    
      /////////////////////////////////////////////////////////////////////////////////////////////////////////
      // linear friction model
      REAL fP  = dem::Parameter::getSingleton().parameter["boundaryFric"] * vfabs(normalForce);
      REAL ks  = 4 * G0 * contactRadius / (2 - poisson);
      tgtForce = prevTgtForce + ks * (-tgtDispInc); // prevTgtForce read by checkin
    
      Vec fricDampingForce = 0;
      if (vfabs(tgtForce) > fP) // slide case
	tgtForce = fP * TgtDirc;
      else { // adhered/slip case
      
	// obtain tangential damping force
	Vec relaVel = currVeloc + currOmga % ((pt1 + pt2)/2 - currPos);  
	Vec TgtVel  = relaVel - (relaVel * normalDirc) * normalDirc;
	REAL dampCritical = 2 * sqrt(getMass() * ks); // critical damping
	fricDampingForce = 1.0 * dampCritical * (-TgtVel);
      }
    
      /////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Mindlin's model (loading/unloading condition assumed)
      // This model is not recommended as it is impossible to strictly determine loading/unloading condition
      // unless load is known (the case of pure moment rotation).
#ifdef MINDLIN_ASSUMED
      REAL val = 0;
      fP = contactFric * vfabs(normalForce);
      tgtLoading = (prevTgtDisp * tgtDispInc >= 0); 
    
      if (tgtLoading) {       // loading
	if (!prevTgtLoading) { // pre-step is unloading
	  val = 8 * G0 * contactRadius * vfabs(tgtDispInc)/(3 * (2 - poisson) * fP);
	  tgtDispStart = prevTgtDisp;
	}
	else                  // pre-step is loading
	  val = 8 * G0 * contactRadius * vfabs(tgtDisp - tgtDispStart)/(3 * (2 - poisson) * fP);
      
	if (val > 1.0)              
	  tgtForce = fP * TgtDirc;
	else {
	  ks = 4 * G0 * contactRadius / (2 - poisson) * sqrt(1 - val);
	  // incremental method
	  tgtForce = prevTgtForce + ks * (-tgtDispInc); // tgtDispInc determines signs
	  // total value method: tgtForce = fP*(1-pow(1-val, 1.5))*TgtDirc;
	}
      }
      else {                 // unloading
	if (prevTgtLoading) { // pre-step is loading
	  val = 8 * G0 * contactRadius * vfabs(tgtDisp - tgtDispStart)/(3 * (2 - poisson) * fP);
	  tgtPeak = vfabs(prevTgtForce);
	}
	else                 // pre-step is unloading
	  val = 8 * G0 * contactRadius * vfabs(tgtDisp - tgtDispStart)/(3 * (2 - poisson) * fP);
      
	if (val > 1.0 || tgtPeak > fP)  
	  tgtForce = fP * TgtDirc;
	else {
	  ks = 2 * sqrt(2) * G0 * contactRadius/(2 - poisson) * sqrt(1 + pow(1 - tgtPeak/fP, 2.0/3.0) + val);
	  // incremental method
	  tgtForce = prevTgtForce + ks*(-tgtDispInc); // tgtDispInc determines signs
	  // total value method: tgtForce = (tgtPeak-2*fP*(1-sqrt(2)/4*pow(1+ pow(1-tgtPeak/fP,2.0/3.0) + val,1.5)))*TgtDirc;
	}
      }
    
      if (vfabs(tgtForce) > fP) // slice case
	tgtForce = fP * TgtDirc;
      else { // adhered/slip case
      
	// obtain tangential damping force
	Vec relaVel = currVeloc + currOmga * ((pt1 + pt2)/2 - currPos);  
	Vec TgtVel  = relaVel - (relaVel * normalDirc) * normalDirc;
	REAL dampCritical = 2 * sqrt(getMass() * ks); // critical damping
	fricDampingForce = 1.0 * dampCritical * (-TgtVel);
      }
    
#endif
      /////////////////////////////////////////////////////////////////////////////////////////////////////////	
    
      /*
	if (iteration % 100 == 0)  
	debugInf << "Particle.cpp, iter=" << iteration
	<< " normalForce=" << vfabs(normalForce)
	<< " cntDampingForce= " << vfabs(cntDampingForce)
	<< " kn=" << kn
	<< " tgtForce=" << vfabs(tgtForce)
	<< " fricDampingForce=" << vfabs(fricDampingForce)
	<< " ks=" << ks
	<< std::endl;
      */
      /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
      // apply tangential force
      addForce(tgtForce);
      addMoment(((pt1 + pt2)/2 - currPos) % tgtForce); 
    
      // apply tangential damping force for adhered/slip case
      addForce(fricDampingForce);
    
      // update current tangential force and displacement, don't checkout.
      // checkout in rigidBF() ensures BdryTgtMap update after each particles
      // contacting this boundary is processed.
      vtmp.push_back(BoundaryTgt(id, tgtForce, tgtDisp, tgtLoading, tgtDispStart, tgtPeak));
    
    }
  
    plane->getContactInfo().push_back(BdryContact(this, pt1, -normalForce, -tgtForce, penetr));
    // update forces acting on boundary in class Boundary, not here
  }
  
  
  Vec Particle::cylinderRBForce(std::size_t boundaryId, const Cylinder &S, int side) {
    // side == -1, the particles are inside the cylinder
    // side == +1, the particles are outside the cylinder
    REAL x0 = S.getCenter().getX();
    REAL y0 = S.getCenter().getY();
    REAL r = S.getRadius();
    REAL coef2[10] = {1, 1, 0, 0, 0, 0, -2*x0, -2*y0, 0, x0*x0 + y0*y0 - r*r};
    Vec pt1;
    if (!root6(coef, coef2, pt1)) // on the cylinder and within the particle
      return 0; // no contact
    ++contactNum;
    Vec rt[2];
    Vec cz = Vec(S.getCenter().getX(), S.getCenter().getY(), pt1.getZ());
    Vec tmp = pt1 - cz;
    intersectWithLine(pt1, normalize(tmp), rt);
    Vec pt2;
  
    if ((rt[0] - pt1) * tmp*side < 0)
      pt2 = rt[0];
    else
      pt2 = rt[1];
    // Vec pt2 = vfabs(rt[0]-cz)>vfabs(rt[1]-cz)?rt[0]:rt[1];
    REAL radius = getRadius(pt2);
    REAL E0 = 0.5*young/(1 - poisson*poisson);
    REAL R0 = (r*radius)/(r + radius);
    REAL rou = vfabs(pt1 - pt2);
    Vec normalDirc = normalize(pt1 - pt2);
    REAL nfc = sqrt(rou*rou*rou) * sqrt(R0) * 4 * E0/3; // pow(rou,1.5), a serious bug
    Vec normalForce = nfc * normalDirc;
  
    addForce(normalForce);
    addMoment(((pt1 + pt2)/2 - getCurrPos()) % normalForce);	    
  
    return normalForce;
  }

  void Particle::clearFluidGrid() {
    fluidGrid.clear();
  }

  void Particle::recordFluidGrid(std::size_t i, std::size_t j, std::size_t k) {
    std::vector<REAL> vec;
    vec.push_back(static_cast<REAL> (i));
    vec.push_back(static_cast<REAL> (j));
    vec.push_back(static_cast<REAL> (k));
    fluidGrid.push_back(vec);
  }

  void Particle::setDemParticleInSPHParticle(){
    for(std::vector<sph::SPHParticle*>::iterator it=SPHGhostParticleVec.begin(); it!=SPHGhostParticleVec.end(); it++){
	(*it)->setDemParticle(this);	// this is the self pointer
    }
  }

  void Particle::setNULLDemParticleInSPHParticle(){
    for(std::vector<sph::SPHParticle*>::iterator it=SPHGhostParticleVec.begin(); it!=SPHGhostParticleVec.end(); it++){
	(*it)->setNULLDemParticle();	// this is the self pointer
    }
  }

    void Particle::updateSPHGhostParticle(){

	dem::Vec sph_local, sph_curr;
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHGhostParticleVec.begin(); pt!=SPHGhostParticleVec.end(); pt++){
	    sph_local = (*pt)->getLocalPosition();	
	    sph_curr = localToGlobal(sph_local)+currPos;
	    sph_curr.setY(0);
	    (*pt)->setCurrPosition(sph_curr);
	    (*pt)->setCurrVelocity(currVeloc+currOmga%(sph_curr-currPos));
	    (*pt)->updateParticleDensity();
	}
    }

} // namespace dem ends
