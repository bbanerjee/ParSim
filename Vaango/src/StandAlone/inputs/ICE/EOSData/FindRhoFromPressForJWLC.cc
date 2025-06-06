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

  // Pointwise computation of microscopic density
  // P=P(rho,T) is not invertable to get rho=rho(P,T)
  // so I'm using Newton's method to find the rhoM
  // such that 
  // press - (A*exp(-R1*rho0/rhoM) +
  //          B*exp(-R2*rho0/rhoM) + C*pow((rhoM/rho0),1+om)) = 0

  // To use:  Change the values as needed at the beginning of main()
  // To compile:  g++ FindRhoFromPressForJWLC.cc
  // To run:  a.out

  // This will return to the command line the density that corresponds
  // to the constants and conditions provided by the user

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>


int main()
{

  // The following are values to be changed by the user
  double press = 29.5e9;   // CJ Pressure
  double A=5.0999434e11;
  double B=4.0522056e9;
  double C=1.3233369e9;
  double R1=4.06244;
  double R2=0.94846;
  double om=0.35;
  double rho0=1717.;

  double rho_guess = 10000.0;

  double rhoM = min(rho_guess,rho0);

  double epsilon = 1.e-15;
  double delta = 1.;
  double f,df_drho,relfac=.9;
  int count = 0;

  double one_plus_omega = 1.+om;

  while(fabs(delta/rhoM)>epsilon){
    double inv_rho_rat=rho0/rhoM;
    double rho_rat=rhoM/rho0;
    double A_e_to_the_R1_rho0_over_rhoM=A*exp(-R1*inv_rho_rat);
    double B_e_to_the_R2_rho0_over_rhoM=B*exp(-R2*inv_rho_rat);
    double C_rho_rat_tothe_one_plus_omega=C*pow(rho_rat,one_plus_omega);

    f = (A_e_to_the_R1_rho0_over_rhoM +
         B_e_to_the_R2_rho0_over_rhoM + C_rho_rat_tothe_one_plus_omega) - press;

    double rho0_rhoMsqrd = rho0/(rhoM*rhoM);
    df_drho = R1*rho0_rhoMsqrd*A_e_to_the_R1_rho0_over_rhoM
            + R2*rho0_rhoMsqrd*B_e_to_the_R2_rho0_over_rhoM
            + (one_plus_omega/rhoM)*C_rho_rat_tothe_one_plus_omega;

    delta = -relfac*(f/df_drho);
    rhoM+=delta;
    rhoM=fabs(rhoM);
    if(count>=100){

      // The following is here solely to help figure out what was going on
      // at the time the above code failed to converge.  Start over with this
      // copy and print more out.
      delta = 1.;
      rhoM = 2.*rho0;
      while(fabs(delta/rhoM)>epsilon){
       double inv_rho_rat=rho0/rhoM;
       double rho_rat=rhoM/rho0;
       double A_e_to_the_R1_rho0_over_rhoM=A*exp(-R1*inv_rho_rat);
       double B_e_to_the_R2_rho0_over_rhoM=B*exp(-R2*inv_rho_rat);
       double C_rho_rat_tothe_one_plus_omega=C*pow(rho_rat,one_plus_omega);

       f = (A_e_to_the_R1_rho0_over_rhoM +
            B_e_to_the_R2_rho0_over_rhoM +
            C_rho_rat_tothe_one_plus_omega) - press;

       double rho0_rhoMsqrd = rho0/(rhoM*rhoM);
       df_drho = R1*rho0_rhoMsqrd*A_e_to_the_R1_rho0_over_rhoM
                + R2*rho0_rhoMsqrd*B_e_to_the_R2_rho0_over_rhoM
                + (one_plus_omega/rhoM)*C_rho_rat_tothe_one_plus_omega;
  
       delta = -relfac*(f/df_drho);
       rhoM+=delta;
       rhoM=fabs(rhoM);
       if(count>=150){
         //cerr << setprecision(15);
         std::cerr <<  "ERROR:ICE:JWLC::computeRhoMicro not converging. \n";
         std::cerr <<  "press= " << press << "\n";
         std::cerr <<  "delta= " << delta << " rhoM= " << rhoM << " f = " << f 
              <<" df_drho =" << df_drho << " rho_guess =" << rho_guess << "\n";
         
       }
       count++;
      }
    }
    count++;
  }
  std::cout << setprecision(15);
  std::cout << "rho = " <<  rhoM << '\n';
  std::cout << "specific volume = " <<  1./rhoM << '\n';
  std::cout << "relative specific volume = " <<  rho0/rhoM << '\n';
  
}
