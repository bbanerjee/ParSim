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

// To compile:

// > g++ -O3 -o resizeCylinders resizeCylinders.cc

// To run:

// > resizeCylinders.cc

// The code will create a PositionRadius.RS.txt file, and also an xml file
// called Test2D.RS.xml which is compatible for inclusion in a Uintah Problem
// Specification (ups) file.




bool isCylInsideRVE(double partDia, double RVEsize,
                       double xCent, double yCent);

bool isCylCenterInsideRVE(double RVEsize,
                          double xCent, double yCent);

bool doesCylIntersectOthers(double partDia, vector<double> diaLocs,
                            double &gap,
                            double &xCent, double &yCent,
                            std::vector<double> xLocs,
                            std::vector<double> yLocs, int &i_this);

void printCylLocs(vector<vector<double> > xLocs,
                     std::vector<vector<double> > yLocs,
                     std::vector<vector<double> > diaLocs,
                     int n_bins, const double RVEsize, double diam_max);

int main()
{

  // Parameters for user to change - BEGIN
  double RVEsize = 0.1;
  // If you just want to specify a minimum and maximum grain size and
  // have a distribution in between, you only need to change diam_min, diam_max
  // and n_sizes.
  // diam_min is the minimum grain diameter
  double diam_min = 0.000625;
  // diam_max is the maximum grain diameter
  double diam_max = 0.011250;
  // How many sizes of grains do you want to us
  int n_sizes = 10;
  // Parameters for user to change - END IF AN EQUAL SIZE DISTRIBUTION IS O.K.
  //  If you want to specify the distribution more carefully, see below.

  // Part of optimizing the search for intersections
  int n_bins = RVEsize/diam_max;
  std::cout << "n_bins = " << n_bins << std::endl;

  double bin_width = RVEsize/((double) n_bins);

  double RVE_area   =  (RVEsize*RVEsize);

  //Store the locations in n_bins separate vectors, so we have a smaller region
  //to search for intersections
  std::vector<double> xLocs;
  std::vector<double> yLocs;
  std::vector<double> diaLocs;

  //Open file to receive cylinder descriptions
  string infile_name = "Position_Radius.txt";
  ifstream source(infile_name.c_str());
  if(!source){
    std::cerr <<  "File " << infile_name << " can't be opened." << std::endl;
  }
  double x,y,r,d;

  int outOfRVE=0;
  while(source >> x >> y >> r){
    if(isCylCenterInsideRVE(RVEsize,x,y)){
     xLocs.push_back(x);
     yLocs.push_back(y);
     diam_max=max(diam_max,2.*r);
     diaLocs.push_back(2.*r);
    } else{
      outOfRVE++;
    }
  }

  double total_cyl_area_orig = 0.0;
  double total_cyl_area_new  = 0.0;
  std::cout << xLocs.size() << std::endl;

  int numInts=0;
  for(int i = 0;i<xLocs.size();i++){
    d = diaLocs[i];
    total_cyl_area_orig+= 0.25*M_PI*(d*d);;
    int i_this = i;
    double gap=9.e99;

    bool cylsIntersect = doesCylIntersectOthers(d, diaLocs, gap,
                                                xLocs[i],yLocs[i],
                                                xLocs,yLocs,i_this);

    if(cylsIntersect){
      numInts++;
    }
    if(gap>9.e90){
      gap = 0.0;
    }
    diaLocs[i] += 0.95*gap;
    d = diaLocs[i];
    diam_max=max(diam_max,d);
    total_cyl_area_new+= 0.25*M_PI*(d*d);;
  }

  std::cout << "numInts = " << numInts << std::endl;
  std::cout << "Spheres out of RVE = " << outOfRVE << std::endl;
  std::cout << "Total cylinder area orig = " << total_cyl_area_orig << std::endl;
  std::cout << "Total cylinder area new  = " << total_cyl_area_new  << std::endl;
  std::cout << "New Maximum Diameter = " << diam_max << std::endl;

  std::vector<vector<double> > xbinLocs(n_bins);
  std::vector<vector<double> > ybinLocs(n_bins);
  std::vector<vector<double> > dbinLocs(n_bins);

  for(int i = 0; i<xLocs.size(); i++){
    int index = (xLocs[i]/RVEsize)*((double) n_bins);
    xbinLocs[index].push_back(xLocs[i]);
    ybinLocs[index].push_back(yLocs[i]);
    dbinLocs[index].push_back(diaLocs[i]);
  }

  printCylLocs(xbinLocs,ybinLocs,dbinLocs,n_bins,RVEsize,diam_max);

}

bool isCylInsideRVE(double partDia, double RVEsize,
                       double xCent, double yCent)
{

    // Find if the particle fits in the box
    double rad = 0.5*partDia;
    double xMinPartBox = xCent-rad;
    double xMaxPartBox = xCent+rad;
    double yMinPartBox = yCent-rad;
    double yMaxPartBox = yCent+rad;
    if (xMinPartBox >= 0.0 && xMaxPartBox <= RVEsize &&
        yMinPartBox >= 0.0 && yMaxPartBox <= RVEsize) {
      return true;
    }
    return false;
}

bool isCylCenterInsideRVE(double RVEsize,
                          double xCent, double yCent)
{

    // Find if the particle center is in the box
    if (xCent >= 0.0 && xCent <= RVEsize &&
        yCent >= 0.0 && yCent <= RVEsize) {
      return true;
    }
    return false;
}

bool doesCylIntersectOthers(double partDia, vector<double> diaLocs,
                            double &gap,
                            double &xCent, double &yCent,
                            std::vector<double> xLocs,
                            std::vector<double> yLocs, int &i_this)
{
  for(unsigned int i = 0; i<xLocs.size(); i++){
   if(i!=i_this){
    // Compute distance between centers
    double distCent = sqrt((xCent-xLocs[i])*(xCent-xLocs[i]) +
                           (yCent-yLocs[i])*(yCent-yLocs[i]));

    double sumRad = 0.5*(partDia + diaLocs[i]);

    double space = distCent - sumRad;
    gap = min(gap, space);

    if(space < 0.0){
      return true;
    }
   }
  }

  // None of the cyls intersected
  return false;
}

void printCylLocs(vector<vector<double> > xLocs, vector<vector<double> > yLocs,
                  std::vector<vector<double> > diaLocs,
                  int n_bins, const double RVEsize, double diam_max)
{
  //Open file to receive cyl descriptions
  string outfile_name = "Test2D.RS.xml";
  ofstream dest(outfile_name.c_str());
  if(!dest){
    std::cerr <<  "File " << outfile_name << " can't be opened." << std::endl;
  }

  dest << "<?xml version='1.0' encoding='ISO-8859-1' ?>" << std::endl;
  dest << "<Uintah_Include>" << std::endl;
  dest << "<union>\n\n";

  int cylcount = 0;
  for(int k=0;k<n_bins;k++){
    for(unsigned int i = 0; i<xLocs[k].size(); i++){
         dest << "    <cylinder label = \"" << cylcount++ << "\">\n";
         dest << "       <top>[" << xLocs[k][i] << ", " << yLocs[k][i] << ", " << 10000 << "]</top>\n";
         dest << "       <bottom>[" << xLocs[k][i] << ", " << yLocs[k][i] << ", " << -10000.0 << "]</bottom>\n";
         dest << "       <radius>" << 0.5*diaLocs[k][i] << "</radius>\n";
         dest << "    </cylinder>\n";
    }
  }

  dest << "</union>\n\n";
  dest << "</Uintah_Include>" << std::endl;

  string outfile_name2 = "Position_Radius.RS.txt";
  ofstream dest2(outfile_name2.c_str());
  if(!dest2){
    std::cerr <<  "File " << outfile_name << " can't be opened." << std::endl;
  }

  dest2.precision(15);

  for(int k=0;k<n_bins;k++){
   for(unsigned int i = 0; i<xLocs[k].size(); i++){
       dest2 <<  xLocs[k][i] << " " << yLocs[k][i] << " " << 0.5*diaLocs[k][i] << "\n";
   }
  }
}
