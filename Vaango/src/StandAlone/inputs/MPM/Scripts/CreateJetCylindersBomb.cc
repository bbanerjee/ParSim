/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>




// define INSERT_P to generate a cylinder distribution that is compatible with
// the insert particles feature

// define FULL_LENGTH to generate the traditional long string of cylinders

//#define FULL_LENGTH
#undef FULL_LENGTH
#define INSERT_P
//#undef INSERT_P




int main()
{
  string outfile_name = "cylinders.xml";
  ofstream dest(outfile_name.c_str());
  if(!dest){
    std::cerr <<  "File " << outfile_name << " can't be opened." << std::endl;
  }

  outfile_name = "InsertParticles.dat";
  ofstream dest_IP(outfile_name.c_str());
  if(!dest_IP){
    std::cerr <<  "File " << outfile_name << " can't be opened." << std::endl;
  }

  dest << "<?xml version='1.0' encoding='ISO-8859-1' ?>" << std::endl;
  dest << "<Uintah_Include>" << std::endl;

  // Density fits
  double ae = -27739.2e1;  // Tail region of jet
  double be = 8248.14e3;   // Tail region of jet
  double a = -6200.44e1;   // Tip region of jet
  double b = 5466.92e3;    // Tip region of jet

  // Radius fits
  double c = -0.0105897;  // Tail region of jet
  double d =  0.0029218e2;  // Tail region of jet
  double q = -0.00253423; // Tip region of jet
  double r =  0.001756e2;   // Tip region of jet
  double rho_W = 17600.0e3;
  double v_tip=.7830;
  double L0=30.;
  double delT=0.25;
  double total_mass=0.;

  double x=12.;          // Position at which to monitor densities

  dest << "<!--\n" << std::endl;
  dest << "x = " << x << std::endl;
  dest << "a = " << a << std::endl;
  dest << "b = " << b << std::endl;
  dest << "c = " << c << std::endl;
  dest << "d = " << d << std::endl;
  dest << "q = " << q << std::endl;
  dest << "r = " << r << std::endl;
  dest << "rho_W = " << rho_W << std::endl;
  dest << "delT = " << delT << std::endl;
  dest << "-->\n" << std::endl;

  double T0=(x-L0)/v_tip; // Time at which tip is at x
  double t=T0; 
  double X=100.;  // Initialize ref position to large value to get started
  int n=0;
  while(X>5.){  // Continue until the tail of the reference jet is reached
    double F=1.0+v_tip*(t/L0);
    double rho;
    X=x/F;
    double rad = 0.;
    if(X <= 30. && X > 13.){
      rad = q*X + r;
      rho=(a*(x/F) + b)/F;
    }
    else if(X <=13.){
      rad = c*X + d;
      rho=(ae*(x/F) + be)/F;
    }
    else{
     std::cout << "SHOULDN'T GET HERE!  X = " << X << std::endl;
    }

    double vel = v_tip*X/L0;
    // L_int = Length of this element of the jet
    double L_int=vel*delT; 
    // L_seg = Length of the segment representing this element of the jet
    double L_seg=L_int*(rho/rho_W);

    total_mass+=(M_PI*rad*rad)*L_seg*rho_W;
    // Make the radius of the cylinder slightly smaller than the radius
    // from the curve fit.  This is done for a few reasons:
    // 1.  Accounts for a non-flat density profile.
    // 2.  Reduces interference with the tunnel walls, a problem that wouldn't
    //     occur in a particulated jet
    double volume = (M_PI*rad*rad)*L_seg;
    if(X<25.){
      rad*=.9;
      L_seg = volume/(M_PI*rad*rad);
    }
    double elapT = t-T0;
    double init_pos = vel*(-elapT) - 0.1;
    dest_IP << t-T0 << " " << n << " " << " 0.0  1.0 0.0  0.0 " << vel << " 0.0\n";
    // adjust delT as to keep cylinders from shrinking too much
    if(X>25.){
      delT=delT*pow(v_tip/vel,.05);
    } else if(X>20.){
      delT=delT*pow(v_tip/vel,.025);
    } else {
      delT=delT*pow(v_tip/vel,.010);
    } 
    t+=delT;

    // Write cylinder to the cylinders.xml file
    dest << "  <geom_object>\n";
    dest << "    <cylinder label = \"" << n << "\">\n";
#ifdef FULL_LENGTH 
    dest << "       <top>[0.0, " << init_pos << ", 0.0]</top>\n";
    dest << "       <bottom>[0.0, " << init_pos-L_seg << ", 0.0]</bottom>\n";
#endif
#ifdef INSERT_P
    dest << "       <top>[0.0, " << -1. << ", 0.0]</top>\n";
    dest << "       <bottom>[0.0, " << -1. - L_seg << ", 0.0]</bottom>\n";
#endif
    dest << "       <radius>" << rad << "</radius>\n";
    dest << "    </cylinder>\n";
    dest << "    <res>[4,4,1]</res>\n";
#ifdef FULL_LENGTH 
    dest << "    <velocity>[0.0, " << vel << ", 0.0]</velocity>\n";
#endif
#ifdef INSERT_P
    dest << "    <velocity>[0.0, 0.0, 0.0]</velocity>\n";
#endif
    dest << "    <temperature>294.0</temperature>\n";
    dest << "    <color>" << n++ << "</color>\n";
    dest << "  </geom_object>\n\n";
  }
  std::cout << "Total Mass = " << total_mass << std::endl;

  dest << "</Uintah_Include>" << std::endl;
}
