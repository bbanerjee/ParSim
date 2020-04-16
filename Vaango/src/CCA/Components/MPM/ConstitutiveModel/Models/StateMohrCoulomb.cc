#include "StateMohrCoulomb.h"

using namespace Uintah;

StateMohrCoulomb::StateMohrCoulomb()
{
  stress = Vector6::Zero(); 
  strain = Vector7::Zero(); 
  plasticStrain = Vector6::Zero(); 
  microStress = Vector6::Zero();
  microStrain = Vector7::Zero();
  microPlasticStrain = Vector6::Zero();
  state = Vector3::Zero();
  microState = Vector3::Zero();
  state(0) = 1e6; // p0Star, unused
  state(1) = 1e8; // Yield suction, unused
  state(2) = 1.5; // Specific Volume, unused
}

void
StateMohrCoulomb::update(const Vector6& plasticStrainInc, 
                         const Vector7& strainInc,
                         const Vector6& stressInc, 
                         double p0StarInc)
{
  double microcheck;
  for (int i = 0; i < 7; i++) {
    microcheck = strain(i);
    microStrain(i) += strainInc(i);
    strain(i) += microStrain(i);
    microStrain(i) -= (strain(i) - microcheck);
  }

  for (int i = 0; i < 6; i++) {
    microcheck = stress(i);
    microStress(i) += stressInc(i);
    stress(i) += microStress(i);
    microStress(i) -= (stress(i) - microcheck);
  }

  for (int i = 0; i < 6; i++) {
    microcheck = plasticStrain(i);
    microPlasticStrain(i) += plasticStrainInc(i);
    plasticStrain(i) += microPlasticStrain(i);
    microPlasticStrain(i) -= (plasticStrain(i) - microcheck);
  }

  microcheck = state(0);
  microState(0) += p0StarInc;
  state(0) += microState(0);
  microState(0) -= (state(0) - microcheck);

  microcheck = state(2);
  microState(2) += state(2) * (exp(-(strainInc(0) + strainInc(1) + strainInc(2))) - 1);
  // microState(2) -= (strainInc(0)+strainInc(1)+strainInc(2))*state(2);
  state(2) += microState(2); // updated specific volume
  microState(2) -= (state(2) - microcheck);
}

double
StateMohrCoulomb::meanStress() const
{
  return firstInvariant() / 3.0;
}

double
StateMohrCoulomb::shearStress() const
{
  double J2 = secondDevInvariant();
  double shearStress = 0.0;
  if (J2 > 1.0-14) {
    shearStress = std::sqrt(3.0 * J2); 
  }
  return shearStress;
}

double
StateMohrCoulomb::firstInvariant() const
{
  double I1 =  (stress(0) + stress(1) + stress(2));
  return I1;
}

double
StateMohrCoulomb::secondInvariant() const
{
  double I2 = stress(0) * stress(1) + stress(1) * stress(2) +
              stress(2) * stress(0) - stress(3) * stress(3) -
              stress(4) * stress(4) - stress(5) * stress(5);
  return I2;
}

double
StateMohrCoulomb::thirdInvariant() const
{
  double I3 =
    stress(0) * stress(1) * stress(2) + 2 * stress(3) * stress(4) * stress(5) -
    stress(0) * stress(5) * stress(5) -     stress(1) * stress(4) * stress(4) -
    stress(2) * stress(3) * stress(3);
  return I3;
}

double
StateMohrCoulomb::firstDevInvariant() const
{
  return 0.0;
}

double
StateMohrCoulomb::secondDevInvariant() const
{
  double J2 = ((stress(0) - stress(1)) * (stress(0) - stress(1)) +
               (stress(0) - stress(2)) * (stress(0) - stress(2)) +
               (stress(1) - stress(2)) * (stress(1) - stress(2))) / 6.0 +
               (stress(3) * stress(3) + stress(4) * stress(4) + 
                stress(5) * stress(5));
  return J2;
}

double
StateMohrCoulomb::thirdDevInvariant() const
{
  double I1 = firstInvariant();
  double I2 = secondInvariant();
  double I3 = thirdInvariant();
  double J3 = I1 * I1 * I1 * 2.0 / 27.0 - I1 * I2 / 3.0 + I3;

  return J3;
}

double
StateMohrCoulomb::getTheta()
{
  double Theta;
  double Sin3Theta;
  double Third = thirdDevInvariant();
  double Second = secondDevInvariant();
  if (Second != 0)
    Sin3Theta = sqrt(3.0) * (-1.5) * Third / pow(Second, 1.5);
  else
    Sin3Theta = 0;
  if (Sin3Theta > 1)
    Sin3Theta = 1;
  if (Sin3Theta < -1)
    Sin3Theta = -1;
  Theta = (asin(Sin3Theta)) / 3.0;

  return Theta;
}

double
StateMohrCoulomb::getThetaDeg()
{
  double ThetaDeg;
  ThetaDeg = getTheta();
  ThetaDeg = ThetaDeg * 45.0 / atan(1.0);
  if (ThetaDeg < 0)
    ThetaDeg += 360;
  if (ThetaDeg > 360)
    ThetaDeg -= 360;

  return ThetaDeg;
}

double
StateMohrCoulomb::getThetaDeg_0()
{
  double ThetaDeg;
  ThetaDeg = getTheta();
  ThetaDeg = ThetaDeg * 45.0 / atan(1.0);
  ThetaDeg -= 30;
  if (ThetaDeg > 360)
    ThetaDeg -= 360;
  if (ThetaDeg < -360)
    ThetaDeg += 360;

  return ThetaDeg;
}

void
StateMohrCoulomb::read()
{
  // data will be read from file point.dta
  // algorithm to read will be analogical as in bbmmodel.read
  ifstream infile("point.dta", ios_base::in);

  // file opened
  string s;
  int slength = 0, index = 0, line = 0;
  // double temp=0;

  do {

    getline(infile, s, '\n');
    line++;
    // cout << s <<" Prep: In line no:"<<line<<endl;
    // getch ();
    if (!infile.good()) {
      cout << "Wrong Data File";
      break;
    }
  } while (s != "***Start of data***");
  // I ignore file until "Start of data";

  // Reading data - 6 stress+7 strain+ 2 state parameters =15
  double storage[15];
  for (int j = 0; j < 15; j++) {

    getline(infile, s, '\n'); // read whole line
    line++;
    // cout << s <<"In line no:"<<line<<endl;
    // getch ();
    bool notcomment = true;
    if (s == "")
      notcomment = false;
    if ((s[0] == '/') && (s[1] == '/'))
      notcomment = false; // check whether not a comment line or empty line

    if (notcomment) {
      slength = s.length(); // get length of line
      index = s.find(";");  // find where is a ; char
      if (index != 0) {
        s.erase(index, slength - index); // delete all after ;
        storage[j] = atof(s.c_str());    // converse to double
      } else
        cout << "No ; in line:" << line << " May cause errors."
             << endl; // warn about lack of ; in line
    } else
      j--;
    if (!infile.good())
      break;
  }

  // Moving data from storage to object variables

  for (int i = 0; i < 6; i++) {
    stress[i] = storage[i];
    strain[i] = storage[i + 6];
  }
  strain[6] = storage[12];
  state[0] = storage[13];
  state[1] = storage[14];
  // finished
  infile.close(); // close file
  // all done
}


void StateMohrCoulomb::write(){};

void
StateMohrCoulomb::getEigen(double Eigen[3])
{

  BBMMatrix Stress(3, 3), EigenValues(3, 3), EigenVect(3, 3);

  Stress.PutElement(1, 1, stress[0]);
  Stress.PutElement(1, 2, stress[3]);
  Stress.PutElement(1, 3, stress[4]);

  Stress.PutElement(2, 1, stress[3]);
  Stress.PutElement(2, 2, stress[1]);
  Stress.PutElement(2, 3, stress[5]);

  Stress.PutElement(3, 1, stress[4]);
  Stress.PutElement(3, 2, stress[5]);
  Stress.PutElement(3, 3, stress[2]);

  Stress.Eigen(&EigenVect, &EigenValues);
  Eigen[0] = EigenValues.getElement(1, 1);
  Eigen[1] = EigenValues.getElement(2, 2);
  Eigen[2] = EigenValues.getElement(3, 3);
}

void
StateMohrCoulomb::getEigen(double Eigen[3], BBMMatrix* EigenVectors)
{

  BBMMatrix Stress(3, 3), EigenValues(3, 3);

  Stress.PutElement(1, 1, stress[0]);
  Stress.PutElement(1, 2, stress[3]);
  Stress.PutElement(1, 3, stress[4]);

  Stress.PutElement(2, 1, stress[3]);
  Stress.PutElement(2, 2, stress[1]);
  Stress.PutElement(2, 3, stress[5]);

  Stress.PutElement(3, 1, stress[4]);
  Stress.PutElement(3, 2, stress[5]);
  Stress.PutElement(3, 3, stress[2]);

  // cout<<"getEigen: Stress matrix:"<<endl;
  // Stress.Print();

  Stress.Eigen(EigenVectors, &EigenValues);
  // cout<<"getEigen: Stress matrix after Eigen procedure:"<<endl;
  // Stress.Print();

  Eigen[0] = EigenValues.getElement(1, 1);
  Eigen[1] = EigenValues.getElement(2, 2);
  Eigen[2] = EigenValues.getElement(3, 3);
}

void
StateMohrCoulomb::getEigen(BBMMatrix* EigenValues, BBMMatrix* EigenVectors)
{

  BBMMatrix Stress(3, 3);

  Stress.PutElement(1, 1, stress[0]);
  Stress.PutElement(1, 2, stress[3]);
  Stress.PutElement(1, 3, stress[4]);

  Stress.PutElement(2, 1, stress[3]);
  Stress.PutElement(2, 2, stress[1]);
  Stress.PutElement(2, 3, stress[5]);

  Stress.PutElement(3, 1, stress[4]);
  Stress.PutElement(3, 2, stress[5]);
  Stress.PutElement(3, 3, stress[2]);

  Stress.Eigen(EigenVectors, EigenValues);
}

void
StateMohrCoulomb::setStressEigen(double Eigen[3])
{
  stress[0] = Eigen[0];
  stress[1] = Eigen[1];
  stress[2] = Eigen[2];
  stress[3] = 0.0;
  stress[4] = 0.0;
  stress[5] = 0.0;
}

bool
StateMohrCoulomb::checkIfFinite()
{
  bool F = true;
  for (int i = 0; i < 6; i++)
    if (!isfinite(stress[i])) {
      F = false;
      cout << "Stress[" << i << "]=" << stress[i] << endl;
      cout << "_finite(Stress[" << i << "])=" << isfinite(stress[i]) << endl;
    }
  for (int i = 0; i < 7; i++)
    if (!isfinite(strain[i])) {
      F = false;
      cout << "Strain[" << i << "]=" << strain[i] << endl;
      cout << "_finite(Strain[" << i << "])=" << isfinite(strain[i]) << endl;
    }
  for (int i = 0; i < 3; i++)
    if (!isfinite(state[i])) {
      F = false;
      cout << "State[" << i << "]=" << state[i] << endl;
      cout << "_finite(State[" << i << "])=" << isfinite(state[i]) << endl;
    }
  if (!F) {
    cout << "Point internal values incorrect." << endl;
    cout << "Stress:" << stress[0] << " " << stress[1] << " " << stress[2]
         << " " << stress[3] << " " << stress[4] << " " << stress[5] << " "
         << endl;
    cout << "Strain:" << strain[0] << " " << strain[1] << " " << strain[2]
         << " " << strain[3] << " " << strain[4] << " " << strain[5] << " "
         << strain[6] << endl;
    cout << "State variables:" << state[0] << " " << state[1] << " " << state[2]
         << endl;
    cout << "Stress finite:" << isfinite(stress[0]) << " "
         << isfinite(stress[1]) << " " << isfinite(stress[2]) << " "
         << isfinite(stress[3]) << " " << isfinite(stress[4]) << " "
         << isfinite(stress[5]) << " " << endl;
    cout << "strain finite:" << isfinite(strain[0]) << " "
         << isfinite(strain[1]) << " " << isfinite(strain[2]) << " "
         << isfinite(strain[3]) << " " << isfinite(strain[4]) << " "
         << isfinite(strain[5]) << " " << endl;
    cout << "State variables finite:" << isfinite(state[0]) << " "
         << isfinite(state[1]) << " " << isfinite(state[2]) << endl;
  }

  return F;
}
