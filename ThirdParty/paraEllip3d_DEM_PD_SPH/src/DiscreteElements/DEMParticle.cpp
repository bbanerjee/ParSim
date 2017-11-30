#include <Core/Const/Constants.h>
#include <Core/Math/ran.h>
#include <Core/Math/root6.h>
#include <Core/Util/Utility.h>
#include <DiscreteElements/DEMParticle.h>
#include <iostream>

//#define MOMENT
#ifdef MOMENT
const std::size_t START =
  10000; // at which time step to apply moment? for moment rotation test only.
#define SLIP // if defined, stick and slip; otherwise slide.
#endif
// for moment case: timeStep = 5.0e-07; totalSteps = 12000

//#define MINDLIN_ASSUMED

using namespace dem;

std::string 
DEMParticle::getDEMParticleShape(DEMParticleShape shape)
{
  if (shape == DEMParticleShape::ELLIPSOID)
    return "ellipsoid";
  else if (shape == DEMParticleShape::SPHERE)
    return "sphere";
  else if (shape == DEMParticleShape::POLYELLIPSOID)
    return "polyellipsoid";
  else if (shape == DEMParticleShape::CUBE)
    return "cube";

  return "none";
}

DEMParticle::DEMParticleShape
DEMParticle::getDEMParticleShape(const std::string& shape)
{
  if (shape == "ellipsoid") 
    return DEMParticleShape::ELLIPSOID;
  else if (shape == "sphere") 
    return DEMParticleShape::SPHERE;
  else if (shape == "polyellipsoid") 
    return DEMParticleShape::POLYELLIPSOID;
  else if (shape == "cube") 
    return DEMParticleShape::CUBE;

  return DEMParticleShape::NONE;
}

std::string 
DEMParticle::getDEMParticleType(DEMParticle::DEMParticleType type)
{
  std::string str = "";
  switch (type) {
    case DEMParticleType::FREE :
      str = "free";
      break;
    case DEMParticleType::FIXED :
      str = "fixed";
      break;
    case DEMParticleType::ROTATE_ONLY :
      str = "rotate_only";
      break;
    case DEMParticleType::TRANSLATE_Z_ONLY :
      str = "translate_z_only";
      break;
    case DEMParticleType::IMPACT_Z_ONLY :
      str = "impact_z_only";
      break;
    case DEMParticleType::BOUNDARY_FREE :
      str = "boundary_free";
      break;
    case DEMParticleType::TRANSLATE_ONLY :
      str = "translate_only";
      break;
    case DEMParticleType::GHOST :
      str = "ghost";
      break;
    default :
      str = "free";
      break;
  }
  return str;
}

DEMParticle::DEMParticleType 
DEMParticle::getDEMParticleType(const std::string& type)
{
  DEMParticle::DEMParticleType val = DEMParticleType::FREE;
  if (type == "free")
    val = DEMParticleType::FREE;
  else if (type == "fixed")
    val = DEMParticleType::FIXED;
  else if (type == "rotate_only")
    val = DEMParticleType::ROTATE_ONLY;
  else if (type == "translate_z_only")
    val = DEMParticleType::TRANSLATE_Z_ONLY;
  else if (type == "impact_z_only")
    val = DEMParticleType::IMPACT_Z_ONLY;
  else if (type == "boundary_free")
    val = DEMParticleType::BOUNDARY_FREE;
  else if (type == "translate_only")
    val = DEMParticleType::TRANSLATE_ONLY;
  else if (type == "ghost")
    val = DEMParticleType::GHOST;
  else
    val = DEMParticleType::FREE;

  return val;
}

DEMParticle::DEMParticle()
  : d_id(0)
  , d_shape(DEMParticleShape::ELLIPSOID)
  , d_type(DEMParticleType::FREE)
  , d_a(0)
  , d_b(0)
  , d_c(0)
  , d_young(0)
  , d_poisson(0)
  , d_currPos(0)
  , d_prevPos(0)
  , d_currDirecA(0)
  , d_currDirecB(0)
  , d_currDirecC(0)
  , d_prevDirecA(0)
  , d_prevDirecB(0)
  , d_prevDirecC(0)
  , d_currentVelocity(0)
  , d_previousVelocity(0)
  , d_currOmga(0)
  , d_prevOmga(0)
  , d_force(0)
  , d_moment(0)
  , d_prevForce(0)
  , d_prevMoment(0)
  , d_bodyForce(0)
  , d_constForce(0)
  , d_constMoment(0)
  , d_density(0)
  , d_mass(0)
  , d_volume(0)
  , d_momentJ(0)
  , d_kineticEnergy(0)
  , d_contactNum(0)
  , d_inContact(false)
{
  for (double& coef : d_coef) coef = 0;
}

void
DEMParticle::init(bool randomOrientation, bool randomVelocity)
{
  // generate orientation of axle a/b/c using Euler angles

  // angle1=[0,Pi], angle2=[0,2*Pi), angle3=[0,2*Pi]
  REAL angle1 = Pi, angle2 = 0.99999 * Pi * 2, angle3 = Pi * 2; 
  if (randomOrientation) {
    angle1 *= ran(&idum);
    angle2 *= ran(&idum);
    angle3 *= ran(&idum);
  } 

  REAL c1, c2, c3, s1, s2, s3;
  c1 = cos(angle1);
  c2 = cos(angle2);
  c3 = cos(angle3);
  s1 = sin(angle1);
  s2 = sin(angle2);
  s3 = sin(angle3);

  REAL l1, m1, n1, l2, m2, n2, l3, m3, n3;
  l1 = c2 * c3 - c1 * s2 * s3;
  m1 = s2 * c3 + c1 * c2 * s3;
  n1 = s1 * s3;
  l2 = -c2 * s3 - c1 * s2 * c3;
  m2 = -s2 * s3 + c1 * c2 * c3;
  n2 = s1 * c3;
  l3 = s1 * s2;
  m3 = -s1 * c2;
  n3 = c1;

  d_currDirecA = Vec(l1, m1, n1);
  d_currDirecB = Vec(l2, m2, n2);
  d_currDirecC = Vec(l3, m3, n3);
  // d_currDirecC = normalize(cross(d_currDirecA, d_currDirecB));

  d_prevPos = d_currPos;
  d_prevDirecA = d_currDirecA;
  d_prevDirecB = d_currDirecB;
  d_prevDirecC = d_currDirecC;

  if (randomVelocity) {
    d_currentVelocity.setX(2 * ran(&idum) - 1);
    d_currentVelocity.setY(2 * ran(&idum) - 1);
    d_currentVelocity.setZ(-ran(&idum));
    d_previousVelocity = d_currentVelocity;
  } else {
    d_previousVelocity = d_currentVelocity = 0;
  }

  d_prevOmga = d_currOmga = 0;
  d_force = d_prevForce = 0;
  d_moment = d_prevMoment = 0;
  d_constForce = d_constMoment = 0;
  d_density = util::getParam<REAL>("specificG") * 1.0e+3;
  d_volume = 4 / 3.0 * Pi * d_a * d_b * d_c;
  d_mass = d_density * d_volume;
  d_momentJ = Vec(d_mass / 5 * (d_b * d_b + d_c * d_c),
                  d_mass / 5 * (d_a * d_a + d_c * d_c),
                  d_mass / 5 * (d_a * d_a + d_b * d_b));
  d_contactNum = 0;
  d_inContact = false;

  // Assumes body force is in the z- direction **TODO** make more general
  computeBodyForce(Vec(0, 0, -1));

  computeAndSetGlobalCoef();
}

DEMParticle::DEMParticle(std::size_t n, 
                         DEMParticleShape shape, DEMParticleType type,
                         Vec center, REAL r, REAL yng, REAL poi)
  : d_id(n)
  , d_shape(shape)
  , d_type(type)
  , d_a(r)
  , d_b(r)
  , d_c(r)
  , d_young(yng)
  , d_poisson(poi)
  , d_currPos(std::move(center))
{
  init();
}

DEMParticle::DEMParticle(std::size_t n, 
                         DEMParticleShape shape, DEMParticleType type,
                         Vec center, REAL ra, REAL rb,
                         REAL rc, REAL yng, REAL poi)
  : d_id(n)
  , d_shape(shape)
  , d_type(type)
  , d_a(ra)
  , d_b(rb)
  , d_c(rc)
  , d_young(yng)
  , d_poisson(poi)
  , d_currPos(std::move(center))
{
  init();
}

DEMParticle::DEMParticle(std::size_t n, 
                         DEMParticleShape shape, DEMParticleType type,
                         Vec center, Gradation& grad,
                         REAL yng, REAL poi,
                         bool randomOrientation, bool randomRadiusRatio,
                         bool randomVelocity)
  : d_id(n)
  , d_shape(shape)
  , d_type(type)
  , d_young(yng)
  , d_poisson(poi)
  , d_currPos(std::move(center))
{
  // generate particle size in terms of gradation distribution
  REAL sievenum = grad.getSieveNum();
  REAL randnum = ran(&idum);
  for (std::size_t k = 0; k < sievenum; ++k) {
    if (randnum <= grad.getPercent()[sievenum - 1 - k]) {
      // use a for sieving (where a >= b >= c)
      d_a = grad.getSize()[sievenum - 1 - k];
      break;
    }
  }

  if (randomRadiusRatio) {
    grad.setPtclRatioBA(ran(&idum));
    grad.setPtclRatioCA(ran(&idum));
  }

  d_b = d_a * grad.getPtclRatioBA();
  d_c = d_a * grad.getPtclRatioCA();

  init(randomOrientation, randomVelocity);
}

DEMParticle::DEMParticle(std::size_t n, 
                         DEMParticleShape shape, DEMParticleType type,
                         Vec dim, Vec position,
                         Vec dirca, Vec dircb, Vec dircc, REAL yng, REAL poi)
  : d_id(n)
  , d_shape(shape)
  , d_type(type)
  , d_young(yng)
  , d_poisson(poi)
{
  d_a = dim.x();
  d_b = dim.y();
  d_c = dim.z();
  d_currPos = d_prevPos = position;
  d_currDirecA = d_prevDirecA = vcos(dirca);
  d_currDirecB = d_prevDirecB = vcos(dircb);
  d_currDirecC = d_prevDirecC = vcos(dircc);
  d_currentVelocity = d_previousVelocity = 0;
  d_currOmga = d_prevOmga = 0;
  d_force = d_prevForce = 0;
  d_moment = d_prevMoment = 0;
  d_constForce = d_constMoment = 0;
  d_contactNum = 0;
  d_density = util::getParam<REAL>("specificG") * 1.0e3;
  d_volume = 4 / 3.0 * Pi * d_a * d_b * d_c;
  d_mass = d_density * d_volume;
  d_momentJ = Vec(d_mass / 5 * (d_b * d_b + d_c * d_c),
                  d_mass / 5 * (d_a * d_a + d_c * d_c),
                  d_mass / 5 * (d_a * d_a + d_b * d_b));
  d_inContact = false;

  // Assumes body force is in the z- direction **TODO** make more general
  computeBodyForce(Vec(0, 0, -1));

  computeAndSetGlobalCoef();
}

Vec
DEMParticle::globalToLocal(Vec input) const
{
  Vec local;
  Vec lmn = currentAxisA();
  local.setX(dot(lmn, input)); // l1,m1,n1
  lmn = currentAxisB();
  local.setY(dot(lmn, input)); // l2,m2,n2
  lmn = currentAxisC();
  local.setZ(dot(lmn, input)); // l3,m3,n3
  return local;
}

Vec
DEMParticle::localToGlobal(Vec input) const
{
  Vec global;
  Vec lmn = Vec(d_currDirecA.x(), d_currDirecB.x(), d_currDirecC.x());
  global.setX(dot(lmn, input)); // l1,l2,l3
  lmn = Vec(d_currDirecA.y(), d_currDirecB.y(), d_currDirecC.y());
  global.setY(dot(lmn, input)); // m1,m2,n3
  lmn = Vec(d_currDirecA.z(), d_currDirecB.z(), d_currDirecC.z());
  global.setZ(dot(lmn, input)); // n1,n2,n3
  return global;
}

Vec
DEMParticle::globalToLocalPrev(Vec input) const
{
  Vec local;
  Vec lmn = previousAxisA();
  local.setX(dot(lmn, input)); // l1,m1,n1
  lmn = previousAxisB();
  local.setY(dot(lmn, input)); // l2,m2,n2
  lmn = previousAxisC();
  local.setZ(dot(lmn, input)); // l3,m3,n3
  return local;
}

Vec
DEMParticle::localToGlobalPrev(Vec input) const
{
  Vec global;
  Vec lmn = Vec(d_prevDirecA.x(), d_prevDirecB.x(), d_prevDirecC.x());
  global.setX(dot(lmn, input)); // l1,l2,l3
  lmn = Vec(d_prevDirecA.y(), d_prevDirecB.y(), d_prevDirecC.y());
  global.setY(dot(lmn, input)); // m1,m2,n3
  lmn = Vec(d_prevDirecA.z(), d_prevDirecB.z(), d_prevDirecC.z());
  global.setZ(dot(lmn, input)); // n1,n2,n3
  return global;
}

// 1: rotational energy is 1/2(I1*w1^2+I2*w2^2+I3*w3^2), where each term is
// expressed in local frame.
// 2. angular velocities in global frame needs to be converted to those in local
// frame.
REAL
DEMParticle::computeTranslationalEnergy() const
{
  return d_mass * pow(vnormL2(d_currentVelocity), 2) / 2;
}

REAL
DEMParticle::computeRotationalEnergy() const
{
  Vec currLocalOmga = globalToLocal(d_currOmga);

  return d_momentJ.x() * pow(currLocalOmga.x(), 2) / 2 +
         d_momentJ.y() * pow(currLocalOmga.y(), 2) / 2 +
         d_momentJ.z() * pow(currLocalOmga.z(), 2) / 2;
}

REAL
DEMParticle::computeKineticEnergy() const
{
  return computeTranslationalEnergy() + computeRotationalEnergy();
}

REAL
DEMParticle::computePotentialEnergy(REAL ref) const
{
  auto gravity = util::getParam<REAL>("gravAccel");
  //std::cout << "g = " << gravity << " m = " << d_mass 
  //          << " (z - z_ref) = " << d_currPos.z() << "-" << ref << std::endl;
  return gravity * d_mass * (d_currPos.z() - ref);
}

void
DEMParticle::getGlobalCoef(REAL coef[]) const
{
  for (std::size_t i = 0; i < 10; ++i)
    coef[i] = this->d_coef[i];
}

REAL
DEMParticle::surfaceError(Vec pt) const
{
  REAL x = pt.x();
  REAL y = pt.y();
  REAL z = pt.z();
  return d_coef[0] * x * x + d_coef[1] * y * y + d_coef[2] * z * z +
         d_coef[3] * x * y + d_coef[4] * y * z + d_coef[5] * z * x +
         d_coef[6] * x + d_coef[7] * y + d_coef[8] * z + d_coef[9];
}

void
DEMParticle::computeAndSetGlobalCoef()
{
  // d_coef[0]-x^2, d_coef[1]-y^2, d_coef[2]-z^2, d_coef[3]-xy, d_coef[4]-yz,
  // d_coef[5]-zx
  // d_coef[6]-x, d_coef[7]-y, d_coef[8]-z, d_coef[9]-const
  if (d_a == d_b && d_b == d_c) {
    d_coef[0] = 1;
    d_coef[1] = 1;
    d_coef[2] = 1;
    d_coef[3] = 0;
    d_coef[4] = 0;
    d_coef[5] = 0;
    d_coef[6] = -2 * d_currPos.x();
    d_coef[7] = -2 * d_currPos.y();
    d_coef[8] = -2 * d_currPos.z();
    d_coef[9] = pow(vnormL2(d_currPos), 2) - d_a * d_a;
    return;
  }

  Vec v1 = d_currDirecA;
  Vec v2 = d_currDirecB;
  Vec v3 = d_currDirecC;

  REAL X0 = d_currPos.x();
  REAL Y0 = d_currPos.y();
  REAL Z0 = d_currPos.z();

  REAL l1 = v1.x() / d_a;
  REAL m1 = v1.y() / d_a;
  REAL n1 = v1.z() / d_a;
  REAL l2 = v2.x() / d_b;
  REAL m2 = v2.y() / d_b;
  REAL n2 = v2.z() / d_b;
  REAL l3 = v3.x() / d_c;
  REAL m3 = v3.y() / d_c;
  REAL n3 = v3.z() / d_c;

  auto xsq = l1 * l1 + l2 * l2 + l3 * l3 ;
  auto ysq = m1 * m1 + m2 * m2 + m3 * m3 ;
  auto zsq = n1 * n1 + n2 * n2 + n3 * n3 ;
  auto xy = 2 * (l1 * m1 +  l2 * m2 + l3 * m3) ;
  auto yz = 2 * (m1 * n1 +  m2 * n2 + m3 * n3) ;
  auto zx = 2 * (l1 * n1 +  l2 * n2 + l3 * n3) ;
  auto x = 2 * ( -l1 * m1 * Y0 - l1 * n1 * Z0 -
                  l2 * m2 * Y0 - l2 * n2 * Z0 - 
                  l3 * m3 * Y0 - l3 * n3 * Z0 -
                  X0 * (l1 * l1) - X0 * (l2 * l2) - X0 * (l3 * l3) );
  auto y = 2 * ( (-l1 * m1 * X0) - ( l2 * m2 * X0) -
                  ( l3 * m3 * X0) - ( m1 * m1 * Y0) -
                  ( m2 * m2 * Y0) - ( m3 * m3 * Y0) -
                  ( m1 * n1 * Z0) - ( m2 * n2 * Z0) -
                  ( m3 * n3 * Z0) );
  auto z = 2 * ( (-l1 * n1 * X0) - ( l2 * n2 * X0) -
                  ( l3 * n3 * X0) - ( m1 * n1 * Y0) -
                  ( m2 * n2 * Y0) - ( m3 * n3 * Y0) -
                  ( n1 * n1 * Z0) - ( n2 * n2 * Z0) -
                  ( n3 * n3 * Z0) );
  auto c = -1 + 2 * (l1 * m1 * X0 * Y0 +
                      l1 * n1 * X0 * Z0 +
                      m1 * n1 * Y0 * Z0 +
                      l2 * m2 * X0 * Y0 +
                      l2 * n2 * X0 * Z0 +
                      m2 * n2 * Y0 * Z0 +
                      l3 * m3 * X0 * Y0 +
                      l3 * n3 * X0 * Z0 +
                      m3 * n3 * Y0 * Z0) +
                    (l1 * l1) * (X0 * X0) +
                    (l2 * l2) * (X0 * X0) +
                    (l3 * l3) * (X0 * X0) +
                    (m1 * m1) * (Y0 * Y0) +
                    (m2 * m2) * (Y0 * Y0) +
                    (m3 * m3) * (Y0 * Y0) +
                    (n1 * n1) * (Z0 * Z0) +
                    (n2 * n2) * (Z0 * Z0) +
                    (n3 * n3) * (Z0 * Z0);

  REAL fac = 1.0/xsq;
  d_coef[0] = 1;
  d_coef[1] = ysq * fac;
  d_coef[2] = zsq * fac;
  d_coef[3] = xy * fac;
  d_coef[4] = yz * fac;
  d_coef[5] = zx * fac;
  d_coef[6] = x * fac;
  d_coef[7] = y * fac;
  d_coef[8] = z * fac;
  d_coef[9] = c * fac;
}

bool
DEMParticle::intersectWithLine(Vec v, Vec dirc, Vec rt[]) const
{
  REAL x0 = v.x();
  REAL y0 = v.y();
  REAL z0 = v.z();
  REAL p = dirc.x();
  REAL q = dirc.y();
  REAL r = dirc.z();
  REAL a = d_coef[0];
  REAL b = d_coef[1];
  REAL c = d_coef[2];
  REAL d = d_coef[3];
  REAL e = d_coef[4];
  REAL f = d_coef[5];
  REAL g = d_coef[6];
  REAL h = d_coef[7];
  REAL i = d_coef[8];
  REAL j = d_coef[9];

  REAL A =
    a * p * p + b * q * q + c * r * r + d * p * q + e * q * r + f * r * p;
  REAL B = 2 * a * p * x0 + 2 * b * q * y0 + 2 * c * r * z0 + d * p * y0 +
           d * q * x0 + e * q * z0 + e * r * y0 + f * p * z0 + f * r * x0 +
           g * p + h * q + i * r;
  REAL C = a * x0 * x0 + b * y0 * y0 + c * z0 * z0 + d * x0 * y0 + e * y0 * z0 +
           f * z0 * x0 + g * x0 + h * y0 + i * z0 + j;

  REAL delta = B * B - 4 * A * C;
  if (delta < 0) {
    std::cerr << "DEMParticle.cpp: iter=" << g_iteration
              << " delta < 0 in intersectWithLine()" << std::endl;
    return false;
  } else {
    REAL t1 = (-B + sqrt(delta)) / (2 * A);
    REAL t2 = (-B - sqrt(delta)) / (2 * A);

    rt[0].setX(t1 * p + x0);
    rt[0].setY(t1 * q + y0);
    rt[0].setZ(t1 * r + z0);
    rt[1].setX(t2 * p + x0);
    rt[1].setY(t2 * q + y0);
    rt[1].setZ(t2 * r + z0);
    return true;
  }
}

// 1. This member function is coded based on Mathematical equations in local
// frame,
//    x^2/a^2 + y^2/b^2 + z^2/c^2 = 1, seeking appropriate osculating circle
// among
//    an infinite number of  osculating circles passing through the contact
// point.
// 2. r = 2*r1*r2/(r1+r2)
// 3. It is important to eliminate float exceptions in computations, that is,
//    when dz/dx == infinite, coordinate x & z are switched to use dx/dz == 0.
// 4. When a point is close to the equator, for example, fabs(z) == 0,
//    float exception is prone to occurring, then a switch is needed as above.
REAL
DEMParticle::computeRadius(Vec v) const
{
  if (d_a == d_b && d_b == d_c)
    return d_a;

  REAL per = 1.0e-4; // define when a point is close to equator
  REAL ra = d_a;     // semi-axles of ellipsoid
  REAL rb = d_b;
  REAL rc = d_c;

  // get the local coordinates of vector v, the point on the particle's surface
  Vec v1 = d_currDirecA;
  Vec v2 = d_currDirecB;
  Vec v3 = d_currDirecC;

  Vec xx = v - d_currPos;
  REAL x = dot(v1, xx);
  REAL y = dot(v2, xx);
  REAL z = dot(v3, xx);

  REAL tmp;
  if (fabs(z) <= d_c * per) { // switch x & z, use 0 instead of infinity
    tmp = ra; ra = rc; rc = tmp;
    tmp = x; x = z; z = tmp;
    if (fabs(z) <= d_a * per) { // switch y & z, use 0 instead of infinity
      tmp = ra; ra = rb; rb = tmp;
      tmp = y; y = z; z = tmp;
    }
  }

  auto raSq = ra * ra;
  auto rbSq = rb * rb;
  auto rcSq = rc * rc;
  auto rcSq_raSq = rcSq/raSq;
  auto rcSq_rbSq = rcSq/rbSq;
  auto xSq = x * x;
  auto ySq = y * y;
  auto zSq = z * z;

  auto p = -(rcSq_raSq * x) / z;
  auto q = -(rcSq_rbSq * y) / z;
  auto r = -rcSq_raSq * (1 + rcSq_raSq * xSq / zSq) / z;
  auto t = -rcSq_rbSq * (1 + rcSq_rbSq * ySq / zSq) / z;
  auto s = -(rcSq * rcSq * x * y) / (raSq * rbSq * zSq * z);

  auto pSq = p * p; 
  auto qSq = q * q; 

  auto n = std::sqrt( 1 + pSq + qSq);

  auto B = n * ( 2 * p * q * s - (1 + pSq) * t - (1 + qSq) * r);
  auto C = n * n * n * n;

  /*
  // if delta < 0, then it is usually -1.0e-20, caused by computational
  // precision.
    REAL A = r * t - s * s;
    if (B*B-4*A*C < 0){
    //std::cout<< "DEMParticle.cpp: iter=" << iteration
    << " delta < 0 in radius()"
    << " delta=" << B*B-4*A*C
    << " -C/B=" << -C/B
    << std::endl;
    }
  */

  return std::abs(-C / B * 2.0); // 2*r1*r2/(r1 + r2)
}


void
DEMParticle::clearContactForce()
{
  d_forceIDMap.clear();
  d_momentIDMap.clear();

  d_inContact = false;

#ifdef MOMENT
  REAL m[20] = { 1,   10, 20, 30, 40, 50, 60, 70, 80, 90,
                 100, 80, 70, 60, 50, 40, 30, 20, 10, 0 };
#ifdef SLIP
  for (std::size_t i = 0; i < 20; ++i)
    m[i] *= 1.0e-8;
#else
  for (std::size_t i = 0; i < 20; ++i)
    m[i] *= 2.0e-8;
#endif
  std::size_t s[20];
  for (std::size_t i = 0; i < 20; ++i)
    s[i] = START + i * 100;

  for (std::size_t i = 0; i < 19; ++i)
    if (iteration >= s[i] && iteration < s[i + 1])
      d_moment += Vec(0, m[i], 0);
  if (iteration >= s[19])
    d_moment += Vec(0, m[19], 0);
#endif
}

// central difference integration method
void
DEMParticle::update(REAL timeStep)
{

  REAL forceDamp = util::getParam<REAL>("forceDamp");
  REAL momentDamp = util::getParam<REAL>("momentDamp");
  REAL massScale = util::getParam<REAL>("massScale");
  REAL mntScale = util::getParam<REAL>("mntScale");
  REAL pileRate = util::getParam<REAL>("pileRate");

  // First update the force and moment by summing quantities computed during
  // contact (initialize with current d_force & d_moment)
  Vec resultantForce = d_force;;
  auto forceMap = forceIDMap();
  for (auto it = forceMap.cbegin(); it != forceMap.cend(); ++it ) {
    resultantForce += it->second;
  }
  setForce(resultantForce);
  Vec resultantMoment = d_moment;
  auto momentMap = momentIDMap();
  for (auto it = momentMap.cbegin(); it != momentMap.cend(); ++it ) {
    resultantMoment += it->second;
  }
  setMoment(resultantMoment);

  //if (getId() == 284) {
  //  std::cout << "DEMParticle::update:: id = " << d_id 
  //            << " force = " << resultantForce
  //            << " moment = " << resultantMoment << "\n";
  //}

  /*
  if (getId() == 2) {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    //std::cerr << std::setprecision(16) 
              << "Before update: MPI rank " << world_rank 
              << " id = " << getId() << "\n\t"
              << " mass = " << d_mass
              << " massScale = " << massScale 
              << " timeStep = " << timeStep << "\n";
    //std::cout << std::setprecision(16) 
              << "\t id:force = [";
    for (auto it = forceIDMap.cbegin(); it != forceIDMap.cend(); ++it ) {
      //std::cout << "\t" << it->first << ":" << it->second << "\n\t ";
    }
    //std::cout << "\t]\n";
    //std::cout << std::setprecision(16) 
              << " force = " << d_force << "\n"
              << " resultant = " << resultantForce << "\n";
    //std::cout << std::setprecision(16) 
              << "\t currPos = " << d_currPos << "\n\t"
              << " currVeloc = " << d_currentVelocity << "\n\t"
              << " prevVelocity = " << d_previousVelocity << " atf = " << forceDamp*2
              << " currDirecA = " << d_currDirecA << "\n";
   
  }
  */

  // 0-free, 1-fixed, 5-free bounary particle
  // It is important to distinguish global frame from local frame!
  if (getType() == DEMParticle::DEMParticleType::FREE || 
      getType() == DEMParticle::DEMParticleType::BOUNDARY_FREE) {
    Vec prevLocalOmga;
    Vec currLocalOmga;
    Vec localMoment;
    REAL atf = forceDamp * 2;
    REAL atm = momentDamp * 2;

    // force: translational kinetics equations are in global frame
    d_currentVelocity = d_previousVelocity * (2 - atf) / (2 + atf) +
                  d_force / (d_mass * massScale) * timeStep * 2 / (2 + atf);
    d_currPos = d_prevPos + d_currentVelocity * timeStep;

    // moment: angular kinetics (rotational) equations are in local frame,
    // so global values need to be converted to those in local frame when
    // applying equations
    localMoment = globalToLocal(d_moment);
    prevLocalOmga = globalToLocal(d_prevOmga);

    currLocalOmga.setX(prevLocalOmga.x() * (2 - atm) / (2 + atm) +
                       localMoment.x() / (d_momentJ.x() * mntScale) * timeStep *
                         2 / (2 + atm));
    currLocalOmga.setY(prevLocalOmga.y() * (2 - atm) / (2 + atm) +
                       localMoment.y() / (d_momentJ.y() * mntScale) * timeStep *
                         2 / (2 + atm));
    currLocalOmga.setZ(prevLocalOmga.z() * (2 - atm) / (2 + atm) +
                       localMoment.z() / (d_momentJ.z() * mntScale) * timeStep *
                         2 / (2 + atm));

    // convert local angular velocities to those in global frame in order to
    // rotate a particle in global space
    d_currOmga = localToGlobal(currLocalOmga);

    d_currDirecA =
      normalize(rotateVec(d_prevDirecA, d_currOmga * timeStep));
    d_currDirecB =
      normalize(rotateVec(d_prevDirecB, d_currOmga * timeStep));
    d_currDirecC =
      normalize(rotateVec(d_prevDirecC, d_currOmga * timeStep));

    //if (getId() == 284) {
      //int world_rank;
      //MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
      //auto forceMap = forceIDMap();
      //std::cout << std::setprecision(16) 
      //          << "After update: MPI rank " << world_rank 
      //          << " id = " << getId() << "\n\t"
      //          << " force = " << d_force << " mass = " << d_mass
      //          << " massScale = " << massScale 
      //          << " timeStep = " << timeStep << "\n";
      //std::cout << std::setprecision(16) 
      //          << "\t id:force = [";
      //for (auto it = forceMap.cbegin(); it != forceMap.cend(); ++it ) {
      //  std::cout << it->first << ":" << it->second << "\n\t ";
      //}
      //std::cout << "\t]\n";
    //  std::cout << std::setprecision(16) 
    //            << "id = " << d_id
    //            << " currPos = " << d_currPos << " "
    //            << " currVeloc = " << d_currentVelocity << " "
    //            << " prevVelocity = " << d_previousVelocity << " atf = " << atf
    //            << " currDirecA = " << d_currDirecA << "\n";
    //}
  }
#ifdef MOMENT
  // special case 2 (moment): translate first, then rotate
  else if (getType() == DEMParticle::DEMParticleType::ROTATE_ONLY) {
    Vec prevLocalOmga;
    Vec currLocalOmga;
    Vec localMoment;
    REAL atf = forceDamp * 2;
    REAL atm = momentDamp * 2;
    d_currentVelocity = d_previousVelocity * (2 - atf) / (2 + atf) +
                  d_force / (d_mass * massScale) * timeStep * 2 / (2 + atf);
    if (iteration < START)
      d_currPos = d_prevPos + d_currentVelocity * timeStep;

    localMoment = globalToLocal(d_moment);
    prevLocalOmga = globalToLocal(d_prevOmga);

    currLocalOmga.setX(prevLocalOmga.x() * (2 - atm) / (2 + atm) +
                       localMoment.x() / (d_momentJ.x() * mntScale) * timeStep *
                         2 / (2 + atm));
    currLocalOmga.setY(prevLocalOmga.y() * (2 - atm) / (2 + atm) +
                       localMoment.y() / (d_momentJ.y() * mntScale) * timeStep *
                         2 / (2 + atm));
    currLocalOmga.setZ(prevLocalOmga.z() * (2 - atm) / (2 + atm) +
                       localMoment.z() / (d_momentJ.z() * mntScale) * timeStep *
                         2 / (2 + atm));

    if (iteration >= START) {
      d_currOmga = localToGlobal(currLocalOmga);

      d_currDirecA =
        normalize(rotateVec(d_prevDirecA, d_currOmga * timeStep));
      d_currDirecB =
        normalize(rotateVec(d_prevDirecB, d_currOmga * timeStep));
      d_currDirecC =
        normalize(rotateVec(d_prevDirecC, d_currOmga * timeStep));
    }
  }
#endif
  // special case 3 (displacemental ellipsoidal
  // pile): translate in vertical direction only
  else if (getType() == DEMParticle::DEMParticleType::TRANSLATE_Z_ONLY) {
    d_currentVelocity.setX(0);
    d_currentVelocity.setY(0);
    d_currentVelocity.setZ(-pileRate);
    d_currPos = d_prevPos + d_currentVelocity * timeStep;
  }
  // special case 4 (impacting ellipsoidal penetrator): impact
  // with inital velocity in vertical direction only
  else if (getType() == DEMParticle::DEMParticleType::IMPACT_Z_ONLY) {
    REAL atf = forceDamp * 2;
    d_currentVelocity = d_previousVelocity * (2 - atf) / (2 + atf) +
                  d_force / (d_mass * massScale) * timeStep * 2 / (2 + atf);
    d_currentVelocity.setX(0);
    d_currentVelocity.setY(0);
    d_currPos = d_prevPos + d_currentVelocity * timeStep;
  }
  // translation only, no rotation
  else if (getType() == DEMParticle::DEMParticleType::TRANSLATE_ONLY) {
    REAL atf = forceDamp * 2;
    d_currentVelocity = d_previousVelocity * (2 - atf) / (2 + atf) +
                  d_force / (d_mass * massScale) * timeStep * 2 / (2 + atf);
    d_currPos = d_prevPos + d_currentVelocity * timeStep;
  }
  // special case 10: pull out a DEM particle in
  // peri-domain, prescribed constant velocity
  else if (getType() == DEMParticle::DEMParticleType::GHOST) {
    d_currPos = d_prevPos + d_currentVelocity * timeStep;
  }

  // Below is needed for all cases
  // ensure three axles perpendicular to each other, and being unit vector
  if (d_currDirecA == 0)
    d_currDirecA = normalize(cross(d_currDirecB, d_currDirecC));
  if (d_currDirecB == 0)
    d_currDirecB = normalize(cross(d_currDirecC, d_currDirecA));
  if (d_currDirecC == 0)
    d_currDirecC = normalize(cross(d_currDirecA, d_currDirecB));

  d_prevPos = d_currPos;
  d_prevDirecA = d_currDirecA;
  d_prevDirecB = d_currDirecB;
  d_prevDirecC = d_currDirecC;
  d_previousVelocity = d_currentVelocity;
  d_prevOmga = d_currOmga;
  d_prevForce = d_force;
  d_prevMoment = d_moment;

  d_contactNum = 0;
  computeAndSetGlobalCoef(); // every time the particle is updated, the algebra expression is
                // also updated
}

bool
DEMParticle::nearestPTOnPlane(REAL p, REAL q, REAL r, REAL s, Vec& ptnp) const
{
  if (d_a == d_b && d_b == d_c) {
    Vec tnm = Vec(p, q, r) / sqrt(p * p + q * q + r * r);
    // signed distance from particle center to plane
    REAL l_nm =
      (d_currPos.x() * p + d_currPos.y() * q + d_currPos.z() * r + s) /
      sqrt(p * p + q * q + r * r);
    ptnp = d_currPos - l_nm * tnm;
    if ((d_a - fabs(l_nm)) / (2.0 * d_a) >
        util::getParam<REAL>("minAllowableRelativeOverlap")) // intersect
      return true;
    else // no intersect
      return false;
  }

  REAL a = d_coef[0];
  REAL b = d_coef[1];
  REAL c = d_coef[2];
  REAL d = d_coef[3];
  REAL e = d_coef[4];
  REAL f = d_coef[5];
  REAL g = d_coef[6];
  REAL h = d_coef[7];
  REAL i = d_coef[8];
  REAL j = d_coef[9];

  REAL domi = e * e * p * p + 4 * c * d * p * q - 4 * a * c * q * q +
              f * f * q * q - 2 * d * f * q * r + d * d * r * r -
              2 * e * (f * p * q + d * p * r - 2 * a * q * r) -
              4 * b * (c * p * p + r * (-f * p + a * r));

  REAL x =
    (-(f * i * q * q) - 2 * b * i * p * r + f * h * q * r + d * i * q * r +
     2 * b * g * r * r - d * h * r * r - e * e * p * s - 2 * b * f * r * s -
     2 * c * (h * p * q - g * q * q - 2 * b * p * s + d * q * s) +
     e * (i * p * q + h * p * r - 2 * g * q * r + f * q * s + d * r * s)) /
    domi;
  REAL y = (f * i * p * q - 2 * f * h * p * r + d * i * p * r + f * g * q * r -
            2 * a * i * q * r - d * g * r * r + 2 * a * h * r * r -
            f * f * q * s + d * f * r * s +
            2 * c * (h * p * p - g * p * q - d * p * s + 2 * a * q * s) +
            e * (-i * p * p + g * p * r + f * p * s - 2 * a * r * s)) /
           domi;

  REAL z =
    (f * h * p * q - 2 * d * i * p * q - f * g * q * q + 2 * a * i * q * q +
     d * h * p * r + d * g * q * r - 2 * a * h * q * r + d * f * q * s -
     d * d * r * s + e * (-h * p * p + g * p * q + d * p * s - 2 * a * q * s) +
     2 * b * (i * p * p - g * p * r - f * p * s + 2 * a * r * s)) /
    domi;

  ptnp = Vec(x, y, z);

  REAL val = a * x * x + b * y * y + c * z * z + d * x * y + e * y * z +
             f * x * z + g * x + h * y + i * z + j;

  if (val >= 0) // not intersect
    return false;
  else // intersect
    return true;
}

void
DEMParticle::planeRBForce(PlaneBoundary* plane,
                          BoundaryTangentArrayMap& BdryTangentMap,
                          BoundaryTangentArray& vtmp,
                          REAL timeStep,
                          REAL minOverlapFactor, REAL maxOverlapFactor,
                          std::size_t iteration)
{
  // (p, q, r) are in the same direction as the outward normal vector,
  // hence it is not necessary to provide information about which side the
  // particle is about the plane.
  REAL p, q, r, s;
  Vec dirc = normalize(plane->getDirection());
  p = dirc.x();
  q = dirc.y();
  r = dirc.z();

  // plane equation: p(x-x0) + q(y-y0) + r(z-z0)  = 0, 
  // that is, px + qy + rz + s = 0
  s = dot(-dirc, plane->getPosition()); 

  Vec pt1;
  // the particle and the plane does not intersect
  if (!nearestPTOnPlane(p, q, r, s, pt1)) 
    return;

  // if particle and plane intersect:
  ++d_contactNum;
  d_inContact = true;
  Vec rt[2];
  // the line and ellipsoid surface does not intersect
  if (!intersectWithLine( pt1, dirc, rt)) 
    return;

  Vec pt2;
  ///* universal, allow for large overlap
  if (p * rt[0].x() + q * rt[0].y() + r * rt[0].z() + s > 0)
    pt2 = rt[0];
  else
    pt2 = rt[1];
  //*/
  /* not universal, only allow for small overlap
  if (vnormL2(rt[0]-pt1) < vnormL2(rt[1]-pt1) )
    pt2 = rt[0];
  else
    pt2 = rt[1];
  */

  // obtain normal force
  REAL R0 = computeRadius(pt2);
  REAL allowedOverlap = 2.0 * R0 * maxOverlapFactor;
  REAL d_penetration = vnormL2(pt1 - pt2);
  if (d_penetration <= allowedOverlap)
    return;

  // rigid wall has infinite young's modulus
  REAL E0 = d_young / (1 - d_poisson * d_poisson); 
  if (d_penetration > allowedOverlap) {
    std::stringstream inf;
    inf.setf(std::ios::scientific, std::ios::floatfield);
    inf << "DEMParticle.cpp: iter=" << std::setw(8) << iteration
        << "  ptcl=" << std::setw(8) << getId() << "  bdry=" << std::setw(8)
        << static_cast<int>(plane->getId()) << " d_penetration=" << std::setw(OWID) << d_penetration
        << " allow=" << std::setw(OWID) << allowedOverlap << "\n";
    MPI_Status status;
    //std::cout << "inf.str().length = " << inf.str().length()
    //          << " length = " << length << "\n";
    MPI_File_write_shared(overlapInf, const_cast<char*>(inf.str().c_str()),
                          inf.str().length(), MPI_CHAR, &status);

    //d_penetration = allowedOverlap;
  }

  //REAL minMeasurableOverlap = util::getParam<REAL>("minMeasurableOverlap");
  //d_penetration = nearbyint(d_penetration / minMeasurableOverlap) * minMeasurableOverlap;
  REAL contactRadius = sqrt(d_penetration * R0);
  Vec normalDirc = -dirc;
  // pow(d_penetration,1.5), a serious bug
  Vec normalForce =
    sqrt(d_penetration * d_penetration * d_penetration) * sqrt(R0) * 4 * E0 / 3 * normalDirc;

  /*
    //std::cout << ' ' << iteration
    << ' ' << getId()
    << ' ' << plane->BoundaryID
    << ' ' << pt1.x()
    << ' ' << pt1.y()
    << ' ' << pt1.z()
    << ' ' << rt[0].x()
    << ' ' << rt[0].y()
    << ' ' << rt[0].z()
    << ' ' << rt[1].x()
    << ' ' << rt[1].y()
    << ' ' << rt[1].z()
    << ' ' << vnormL2(rt[0]-pt1)
    << ' ' << vnormL2(rt[1]-pt1)
    << ' ' << d_penetration
    << std::endl;
  */

  // apply normal force
  addForce(normalForce);
  addMoment(cross((pt1 + pt2) / 2 - d_currPos, normalForce));

  // obtain normal damping force
  auto dampingCoeff = util::getParam<REAL>("contactDamp");
  Vec veloc2 = currentVelocity() + cross(currentAngularVelocity(), ((pt1 + pt2) / 2 - currentPosition()));
  REAL kn = pow(6 * vnormL2(normalForce) * R0 * pow(E0, 2), 1.0 / 3.0);
  REAL dampCritical = 2 * sqrt(mass() * kn); // critical damping
  Vec contactDampingForce = dampingCoeff * dampCritical *
                        dot(-veloc2, normalDirc) * normalDirc;

  // apply normal damping force
  addForce(contactDampingForce);
  addMoment(cross(((pt1 + pt2) / 2 - d_currPos), contactDampingForce));

  std::cout << "Boundary-Particle:" << d_force << " f_n = " << normalForce
            << " f_d = " << contactDampingForce << "\n";
  Vec tangentForce = 0;
  if (util::getParam<REAL>("boundaryFric") != 0) {
    // checkin previous tangential force and displacement
    Vec prevTangentForce;
    Vec prevTangentDisp;
    // bool prevTangentLoading = true;
    Vec tangentDispStart;
    REAL tangentPeak = 0;

    bool tangentLoading = true;
    auto boundaryID = static_cast<size_t>(plane->getId());
    for (const auto& tangent : BdryTangentMap[boundaryID]) {
      if (d_id == tangent.particleId) {
        prevTangentForce = tangent.tangentForce;
        prevTangentDisp = tangent.tangentDisp;
        // prevTangentLoading = it->tangentLoading;
        tangentDispStart = tangent.tangentDispStart;
        tangentPeak = tangent.tangentPeak;
        break;
      }
    }

    // obtain tangtential force
    REAL G0 = d_young / 2 / (1 + d_poisson);
    // Vr = Vb + w (crossdot) r, each item needs to be in either global or local
    // frame; here global frame is used for better convenience.
    Vec relaDispInc =
      (d_currentVelocity + cross(d_currOmga, ((pt1 + pt2) / 2 - d_currPos))) * timeStep;
    Vec tangentDispInc = relaDispInc - dot(relaDispInc, normalDirc) * normalDirc;
    Vec tangentDisp = prevTangentDisp + tangentDispInc; // prevTangentDisp read by checkin
    Vec TangentDirc;

    if (vnormL2(tangentDisp) == 0)
      TangentDirc = 0;
    else
      TangentDirc = normalize(-tangentDisp); // TangentDirc points along tangential forces
                                     // exerted on particle 1

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // linear friction model
    REAL fP = util::getParam<REAL>("boundaryFric") * vnormL2(normalForce);
    REAL ks = 4 * G0 * contactRadius / (2 - d_poisson);
    tangentForce =
      prevTangentForce + ks * (-tangentDispInc); // prevTangentForce read by checkin

    Vec fricDampingForce = 0;
    if (vnormL2(tangentForce) > fP) // slide case
      tangentForce = fP * TangentDirc;
    else { // adhered/slip case

      // obtain tangential damping force
      Vec relaVel = d_currentVelocity + cross(d_currOmga, ((pt1 + pt2) / 2 - d_currPos));
      Vec TangentVel = relaVel - dot(relaVel, normalDirc) * normalDirc;
      REAL dampCritical = 2 * sqrt(mass() * ks); // critical damping
      fricDampingForce = 1.0 * dampCritical * (-TangentVel);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mindlin's model (loading/unloading condition assumed)
// This model is not recommended as it is impossible to strictly determine
// loading/unloading condition
// unless load is known (the case of pure moment rotation).
#ifdef MINDLIN_ASSUMED
    REAL val = 0;
    fP = contactFric * vnormL2(normalForce);
    tangentLoading = (prevTangentDisp * tangentDispInc >= 0);

    if (tangentLoading) {        // loading
      if (!prevTangentLoading) { // pre-step is unloading
        val = 8 * G0 * contactRadius * vnormL2(tangentDispInc) /
              (3 * (2 - d_poisson) * fP);
        tangentDispStart = prevTangentDisp;
      } else // pre-step is loading
        val = 8 * G0 * contactRadius * vnormL2(tangentDisp - tangentDispStart) /
              (3 * (2 - d_poisson) * fP);

      if (val > 1.0)
        tangentForce = fP * TangentDirc;
      else {
        ks = 4 * G0 * contactRadius / (2 - d_poisson) * sqrt(1 - val);
        // incremental method
        tangentForce =
          prevTangentForce + ks * (-tangentDispInc); // tangentDispInc determines signs
        // total value method: tangentForce = fP*(1-pow(1-val, 1.5))*TangentDirc;
      }
    } else {                // unloading
      if (prevTangentLoading) { // pre-step is loading
        val = 8 * G0 * contactRadius * vnormL2(tangentDisp - tangentDispStart) /
              (3 * (2 - d_poisson) * fP);
        tangentPeak = vnormL2(prevTangentForce);
      } else // pre-step is unloading
        val = 8 * G0 * contactRadius * vnormL2(tangentDisp - tangentDispStart) /
              (3 * (2 - d_poisson) * fP);

      if (val > 1.0 || tangentPeak > fP)
        tangentForce = fP * TangentDirc;
      else {
        ks = 2 * sqrt(2) * G0 * contactRadius / (2 - d_poisson) *
             sqrt(1 + pow(1 - tangentPeak / fP, 2.0 / 3.0) + val);
        // incremental method
        tangentForce =
          prevTangentForce + ks * (-tangentDispInc); // tangentDispInc determines signs
        // total value method: tangentForce = (tangentPeak-2*fP*(1-sqrt(2)/4*pow(1+
        // pow(1-tangentPeak/fP,2.0/3.0) + val,1.5)))*TangentDirc;
      }
    }

    if (vnormL2(tangentForce) > fP) // slice case
      tangentForce = fP * TangentDirc;
    else { // adhered/slip case

      // obtain tangential damping force
      Vec relaVel = d_currentVelocity + d_currOmga * ((pt1 + pt2) / 2 - d_currPos);
      Vec TangentVel = relaVel - (relaVel * normalDirc) * normalDirc;
      REAL dampCritical = 2 * sqrt(mass() * ks); // critical damping
      fricDampingForce = 1.0 * dampCritical * (-TangentVel);
    }

#endif
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    /*
        if (iteration % 100 == 0)
        //std::cout << "DEMParticle.cpp, iter=" << iteration
        << " normalForce=" << vnormL2(normalForce)
        << " contactDampingForce= " << vnormL2(contactDampingForce)
        << " kn=" << kn
        << " tangentForce=" << vnormL2(tangentForce)
        << " fricDampingForce=" << vnormL2(fricDampingForce)
        << " ks=" << ks
        << std::endl;
      */
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    // apply tangential force
    addForce(tangentForce);
    addMoment(cross(((pt1 + pt2) / 2 - d_currPos), tangentForce));

    // apply tangential damping force for adhered/slip case
    addForce(fricDampingForce);

    // update current tangential force and displacement, don't checkout.
    // checkout in rigidBF() ensures BdryTangentMap update after each particles
    // contacting this boundary is processed.
    vtmp.push_back(BoundaryTangent(d_id, tangentForce, tangentDisp, tangentLoading,
                                   tangentDispStart, tangentPeak));
  }

  plane->getBoundaryContacts().push_back(
    BoundaryContact(this, pt1, -normalForce, -tangentForce, d_penetration));
  // update forces acting on boundary in class Boundary, not here
}

Vec
DEMParticle::cylinderRBForce(std::size_t BoundaryID, const Cylinder& S, int side)
{
  // side == -1, the particles are inside the cylinder
  // side == +1, the particles are outside the cylinder
  REAL x0 = S.center().x();
  REAL y0 = S.center().y();
  REAL r = S.radius();
  REAL d_coef2[10] = { 1, 1,       0,       0, 0,
                       0, -2 * x0, -2 * y0, 0, x0 * x0 + y0 * y0 - r * r };
  Vec pt1;
  if (!root6(d_coef, d_coef2, pt1)) // on the cylinder and within the particle
    return 0;                       // no contact
  ++d_contactNum;
  Vec rt[2];
  Vec cz = Vec(S.center().x(), S.center().y(), pt1.z());
  Vec tmp = pt1 - cz;
  intersectWithLine(pt1, normalize(tmp), rt);
  Vec pt2;

  if (dot((rt[0] - pt1), tmp) * side < 0)
    pt2 = rt[0];
  else
    pt2 = rt[1];
  // Vec pt2 = vnormL2(rt[0]-cz)>vnormL2(rt[1]-cz)?rt[0]:rt[1];
  REAL radius = computeRadius(pt2);
  REAL E0 = 0.5 * d_young / (1 - d_poisson * d_poisson);
  REAL R0 = (r * radius) / (r + radius);
  REAL rou = vnormL2(pt1 - pt2);
  Vec normalDirc = normalize(pt1 - pt2);
  REAL nfc = sqrt(rou * rou * rou) * sqrt(R0) * 4 * E0 /
             3; // pow(rou,1.5), a serious bug
  Vec normalForce = nfc * normalDirc;

  addForce(normalForce);
  addMoment(cross(((pt1 + pt2) / 2 - currentPosition()), normalForce));

  return normalForce;
}

void
DEMParticle::clearFluidGrid()
{
  d_fluidGrid.clear();
}

void
DEMParticle::recordFluidGrid(std::size_t i, std::size_t j, std::size_t k,
                          REAL volFrac)
{
  std::vector<REAL> vec;
  vec.push_back(static_cast<REAL>(i));
  vec.push_back(static_cast<REAL>(j));
  vec.push_back(static_cast<REAL>(k));
  vec.push_back(volFrac);
  d_fluidGrid.push_back(vec);
}

bool 
DEMParticle::containsPoint(const dem::Vec& point,
                           const dem::Vec& dem_pos,
                           const REAL& bufferLength,
                           dem::Vec& localCoord,
                           bool& insideGhostLayer) const
{
  REAL radius_a = d_a;
  REAL radius_b = d_b;
  REAL radius_c = d_c;
  REAL inner_a = d_a - (bufferLength < 0 ? 0 : bufferLength); 
  REAL inner_b = d_b - (bufferLength < 0 ? 0 : bufferLength); 
  REAL inner_c = d_c - (bufferLength < 0 ? 0 : bufferLength);

  localCoord = globalToLocal(point - dem_pos);
  REAL outer_x = localCoord.x()/radius_a;
  REAL outer_y = localCoord.y()/radius_b;
  REAL outer_z = localCoord.z()/radius_c;
  REAL inner_x = localCoord.x()/(inner_a + 1.0e-20);
  REAL inner_y = localCoord.y()/(inner_b + 1.0e-20);
  REAL inner_z = localCoord.z()/(inner_c + 1.0e-20);

  REAL outer_rad_sq = outer_x*outer_x + outer_y*outer_y + outer_z*outer_z;
  REAL inner_rad_sq = inner_x*inner_x + inner_y*inner_y + inner_z*inner_z;

  insideGhostLayer = false;
  if (outer_rad_sq <= 1) {
    if (inner_rad_sq > 1 || inner_a <= 0 || inner_b <= 0 || inner_c <= 0) {
      insideGhostLayer = true;
    }
    return true;
  }
  return false;
}

/**
 *  Approximate shortest distance by projection of the point on to the ellipsoid
 *  (see Zimmermann and Svoboda, "Approximation of Euclidean Distance of Point from Ellipse")
 *
 *  d = || (p - p/||T.p||) ||
 *    where p is the point in the local coord system of the ellipsoid
 *    and   T is the transformation matrix from the ellipsoid to the unit sphere
 *          T = [[1/a 0 0],[0 1/b 0];[0 0 1/c]]
 */
REAL 
DEMParticle::shortestDistToBoundary(const Vec& point) const
{
  Vec radii(d_a, d_b, d_c);
  auto centroid = d_currPos;
  auto pp = globalToLocal(point - centroid);
  auto norm_T_dot_p = (pp/radii).length();
  auto proj = pp - pp/norm_T_dot_p;
  return proj.length();
}

void
DEMParticle::dragForce()
{
  REAL Cd = util::getParam<REAL>("Cd");
  REAL rho = util::getParam<REAL>("fluidDensity");

  REAL ux = d_currentVelocity.x();
  REAL uy = d_currentVelocity.y();
  REAL uz = d_currentVelocity.z();
  Vec globalDelta(fabs(ux) * ux, fabs(uy) * uy, fabs(uz) * uz);
  Vec localDelta = globalToLocal(globalDelta);
  Vec localForce(0);
  // localDelta needs to project in local frame in order to calculate local drag
  // forces
  localForce.setX(-0.5 * rho * localDelta.x() * Cd * Pi * d_b * d_c);
  localForce.setY(-0.5 * rho * localDelta.y() * Cd * Pi * d_c * d_a);
  localForce.setZ(-0.5 * rho * localDelta.z() * Cd * Pi * d_a * d_b);
  Vec globalForce = localToGlobal(localForce);
  addForce(globalForce);
}

