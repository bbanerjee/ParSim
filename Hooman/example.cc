#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;


int main()
{

float a, b;

a=3.14567546;
b=3.14568;

if (abs(a-b)<0.000001)
   {cout<<"a and b are approximately equal to each other"<<endl;}
else
   {cout<<"a and b are not equal"<<endl;}

return 0;
}
