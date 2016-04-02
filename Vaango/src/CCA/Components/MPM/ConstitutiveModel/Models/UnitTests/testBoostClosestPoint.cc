#include <iostream>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/assign.hpp>

#include <boost/numeric/conversion/bounds.hpp>
#include <boost/foreach.hpp>

typedef boost::geometry::model::d2::point_xy<double> point_type;
typedef boost::geometry::model::polygon<point_type> polygon_type;

std::vector<double> linspace(double start, double end, int num)
{
  double delta = (end - start) / (double)num;

  std::vector<double> linspaced;
  for(int i=0; i < num+1; ++i)
    {
      linspaced.push_back(start + delta * (double) i);
    }
  return linspaced;
}

std::vector<point_type> getYieldSurfacePoints(double capX, double zeta)
{
  double PEAKI1 = 1.0e3;
  double FSLOPE = 0.453;
  double STREN = 1.0e7;
  double YSLOPE = 0.31;
  double BETA = 1.0;
  double CR = 0.5;
  //double T1 = 0.0;
  //double T2 = 0.0;

  //double capX = -1000.0;
  //double zeta = 0.0;

  double sqrtKG = 0.500805;

  // Compute a1, a2, a3, a4
  double a1 = STREN;
  double a2 = (FSLOPE-YSLOPE)/(STREN-YSLOPE*PEAKI1);
  double a3 = (STREN-YSLOPE*PEAKI1)*std::exp(-a2*PEAKI1);
  double a4 = YSLOPE;

  // Compute kappa
  double kappa = PEAKI1 - CR*(PEAKI1 - capX);

  // Set up I1 values
  int num_points = 10;
  std::vector<double> I1_vec = linspace(0.99999*capX+zeta, 0.99999*PEAKI1+zeta, num_points);
  std::vector<double> J2_vec;
  for (auto I1 : I1_vec) {

    // Compute F_f
    double I1_minus_zeta = I1 - zeta;
    double Ff = a1 - a3*std::exp(a2*I1_minus_zeta) - a4*(I1_minus_zeta);
    double Ff_sq = Ff*Ff;

    // Compute Fc
    double Fc_sq = 1.0;
    if ((I1_minus_zeta < kappa) && (capX < I1_minus_zeta)) {
      double ratio = (kappa - I1_minus_zeta)/(kappa - capX);
      Fc_sq = 1.0 - ratio*ratio;
    }

    // Compute J2
    double J2 = Ff_sq*Fc_sq;
    J2_vec.push_back(J2);
  }

  // Convert I1 vs J2 to r' vs. z
  std::vector<double> z_vec;
  for (auto I1 : I1_vec) {
    z_vec.push_back(I1/std::sqrt(3.0));
  }

  std::vector<double> rprime_vec;
  for (auto J2 : J2_vec) {
    rprime_vec.push_back(BETA*std::sqrt(2.0*J2)*sqrtKG);
  }

  // Create a point_type vector
  std::vector<point_type> polyline;
  auto z_iter = z_vec.begin();
  auto r_iter = rprime_vec.begin();
  while (z_iter != z_vec.end() || r_iter != rprime_vec.end()) {
    double z = *z_iter;
    double r = *r_iter;
    std::cout << "(" << z << "," << r << ")";
    polyline.push_back(point_type(z, r));
    ++z_iter; ++r_iter;
  }
  std::cout << std::endl;

  auto rev_z_iter = z_vec.rbegin();
  auto rev_r_iter = rprime_vec.rbegin();
  while (rev_z_iter != z_vec.rend() || rev_r_iter != rprime_vec.rend()) {
    double z = *rev_z_iter;
    double r = *rev_r_iter;
    std::cout << "(" << z << "," << -r << ")";
    polyline.push_back(point_type(z, -r));
    ++rev_z_iter; ++rev_r_iter;
  }
  polyline.push_back(point_type(*(z_vec.begin()), *(rprime_vec.begin())));
  std::cout << std::endl;

  return polyline;

}

// From https://www.mathworks.com/matlabcentral/fileexchange/19398-distance-from-a-point-to-polygon/content/p_poly_dist.m
point_type findClosestPoint(const point_type& p, const std::vector<point_type>& poly)
{
  // Get point coordinates
  double xx = boost::geometry::get<0>(p);
  double yy = boost::geometry::get<1>(p);

  std::vector<double> A, B, C, AB, VV;
  std::vector<point_type> XP;

  auto iterStart = poly.begin();
  auto iterEnd   = poly.end();
  auto iterNext = iterStart;
  ++iterNext;
  for ( ; iterNext != iterEnd; ++iterStart, ++iterNext) {
    double xstart = boost::geometry::get<0>(*iterStart);
    double ystart = boost::geometry::get<1>(*iterStart);
    double xnext = boost::geometry::get<0>(*iterNext);
    double ynext = boost::geometry::get<1>(*iterNext);

    // segments that connect the vertices
    double xab = xnext - xstart;
    double yab = ynext - ystart;

    // segment length (squared)
    double abSq = xab*xab + yab*yab;

    // find the projection of point p = (x,y) on each segment
    double xpa = xx - xstart;
    double ypa = yy - ystart;

    // find t = (p - a)/(b - a);
    double pa_dot_ab = xpa*xab + ypa*yab;
    double tt = pa_dot_ab/abSq;

    std::cout << " tt = " << tt << " xp = " <<  xstart + tt*xab << " yp = " << ystart + tt*yab
              << std::endl;

    // Find projction point
    if (!(tt < 0.0 || tt > 1.0)) {
      double xp = xstart + tt*xab;
      double yp = ystart + tt*yab;
      XP.push_back(point_type(xp, yp));
    }
  }

  if (!(XP.size() > 0)) {
    std::cout << "No closest point" << std::endl;
   
    point_type min_p;
    double min_d = boost::numeric::bounds<double>::highest();
    for (auto pa : poly) {
      double xa = boost::geometry::get<0>(pa);
      double ya = boost::geometry::get<1>(pa);
      double dist = (xx - xa)*(xx - xa) + (yy - ya)*(yy - ya);
      if (dist < min_d) {
        min_d = dist;
        min_p = pa; 
      }
    }
    return min_p;
  }

  point_type min_p;
  double min_d = boost::numeric::bounds<double>::highest();
  BOOST_FOREACH(point_type const& xp, XP)
    {
        double d = boost::geometry::comparable_distance(p, xp);
        if (d < min_d)
        {
            min_d = d;
            min_p = xp;
        }
    }
    
    std::cout 
        << "Closest: " << boost::geometry::dsv(min_p) << std::endl
        << "At: " << boost::geometry::distance(p, min_p) << std::endl;
  return min_p;
}

int main()
{

    point_type p(-728.967, 0.0);

    double capX = -1000.0;
    double zeta = 0.0;
    std::vector<point_type> poly_points = getYieldSurfacePoints(capX, zeta);

    polygon_type yield_surface;
    boost::geometry::assign_points(yield_surface, poly_points);
    std::cout << "Polygon = " << boost::geometry::dsv(yield_surface) << std::endl;

    double d = boost::geometry::distance(p, yield_surface);
    std::cout << "Closest distance: " << d << std::endl;

    point_type xp = findClosestPoint(p, poly_points);
    std::cout << " x = " << boost::geometry::get<0>(xp)
              << " y = " << boost::geometry::get<1>(xp) << std::endl;

    point_type p1(-242.201, 1371.81);
    point_type xp1 = findClosestPoint(p1, poly_points);
    std::cout << " x = " << boost::geometry::get<0>(xp1)
              << " y = " << boost::geometry::get<1>(xp1) << std::endl;

    capX = -1000.0;
    zeta = -28985.9;
    std::vector<point_type> poly_points1 = getYieldSurfacePoints(capX, zeta);

    point_type p2(-6290.1, 658.137);
    point_type xp2 = findClosestPoint(p2, poly_points1);
    std::cout << " x = " << boost::geometry::get<0>(xp2)
              << " y = " << boost::geometry::get<1>(xp2) << std::endl;

    return 0;
}
