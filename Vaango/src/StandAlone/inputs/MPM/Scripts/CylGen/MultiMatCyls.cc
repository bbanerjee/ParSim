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

// Distribute cylinders into various materials according to a user specified
// size distribution.  A set of cylinders is provided in a file called
// Position_Radius.txt.  These are then printed into files in both xml format
// and in plain text

// To compile:

//                  g++ -O3 -o MultiMatCyls MultiMatCyls.cc

// To run:

//                  MultiMatCyls


#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <vector>



void printCylLocs(vector<vector<double> > xLocs, vector<vector<double> > yLocs,
                  std::vector<vector<double> > diaLocs,
                  int n_bins, double RVEsize, double diam_max, int matl_num);

int main()
{
  // Parameters for user to change - BEGIN
  double RVEsize = 0.1;
  double diam_max = 0.01375;
  int n_bins  = 10;
  int n_matls = 10;

  //Open file to receive cylinder descriptions
  string infile_name = "Position_Radius.txt";
  ifstream source(infile_name.c_str());
  if(!source){
    std::cerr <<  "File " << infile_name << " can't be opened." << std::endl;
  }

  double x,y,r;

  std::vector<double> xpos;
  std::vector<double> ypos;
  std::vector<double> radius;

  // Read in the cylinders
  int num_cyls=0;
  while(source >> x >> y >> r){
     num_cyls++;
     xpos.push_back(x);
     ypos.push_back(y);
     radius.push_back(r);
  }

  std::vector<vector<double> > xmatlLocs(n_matls);
  std::vector<vector<double> > ymatlLocs(n_matls);
  std::vector<vector<double> > rmatlLocs(n_matls);
  
  // Distribute the cylinders to the appropriate materials based on
  for(int n = 0; n<num_cyls; n++){
    int i = n%n_matls;
    xmatlLocs[i].push_back(xpos[n]);
    ymatlLocs[i].push_back(ypos[n]);
    rmatlLocs[i].push_back(radius[n]);
  }

  // Bin the cylinders according to their x-position.  No real good reason
  // to do this other than compatibility with an existing print function
  for(int i = 0; i<n_matls; i++){
    std::vector<vector<double> > xbinLocs(n_bins);
    std::vector<vector<double> > ybinLocs(n_bins);
    std::vector<vector<double> > rbinLocs(n_bins);

    for(int k = 0; k<xmatlLocs[i].size(); k++){
      int index = (xmatlLocs[i][k]/RVEsize)*((double) n_bins);
      xbinLocs[index].push_back(xmatlLocs[i][k]);
      ybinLocs[index].push_back(ymatlLocs[i][k]);
      rbinLocs[index].push_back(rmatlLocs[i][k]);
    }

    printCylLocs(xbinLocs,ybinLocs,rbinLocs,n_bins,RVEsize,diam_max,i);
  }

}

void printCylLocs(vector<vector<double> > xLocs,
                  std::vector<vector<double> > yLocs,
                  std::vector<vector<double> > radLocs, int n_bins,
                  double RVEsize, double diam_max, int matl_num)
{
     std::stringstream out;
    string s;

    out << matl_num;
    s = out.str();
   
    //Open file to receive cylinder descriptions
    string outfile_name = "Test2D." + s + ".xml";
    ofstream dest(outfile_name.c_str());
    if(!dest){
      std::cerr <<  "File " << outfile_name << " can't be opened." << std::endl;
    }

    dest << "<?xml version='1.0' encoding='ISO-8859-1' ?>" << std::endl;
    dest << "<Uintah_Include>" << std::endl;
    dest << "<intersection>\n";
    dest << "  <box>\n";
    dest << "    <min>[0.0, 0.0, -10000.0]</min>" << std::endl;
    dest << "    <max>[" << RVEsize << ", " << RVEsize << ",  10000.0]</max>" << std::endl;
    dest << "  </box>\n\n";

    dest << "  <union>\n";
    for(int k=0;k<n_bins;k++){
      if(xLocs[k].size()>0){
        for(unsigned int i = 0; i<xLocs[k].size(); i++){
             dest << "    <cylinder>\n";
             dest << "       <top>[" << xLocs[k][i] << ", " << yLocs[k][i] << ", " << 10000 << "]</top>\n";
             dest << "       <bottom>[" << xLocs[k][i] << ", " << yLocs[k][i] << ", " << -10000.0 << "]</bottom>\n";
             dest << "       <radius>" << radLocs[k][i] << "</radius>\n";
             dest << "    </cylinder>\n";
        }
      }
    }
    dest << "  </union>\n\n";
    dest << " </intersection>\n\n";

    dest << "</Uintah_Include>" << std::endl;

    string outfile_name2 = "Position_RadiusNew." + s + ".txt";
    ofstream dest2(outfile_name2.c_str());
    if(!dest2){
      std::cerr <<  "File " << outfile_name << " can't be opened." << std::endl;
    }

    for(int k=0;k<n_bins;k++){
      if(xLocs[k].size()>0){
        for(unsigned int i = 0; i<xLocs[k].size(); i++){
         dest2 <<  xLocs[k][i] << " " << yLocs[k][i] << " " << radLocs[k][i] << "\n";
        }
      }
    }
}
