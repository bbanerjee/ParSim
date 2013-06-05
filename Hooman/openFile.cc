#include <iostream>
#include <fstream>
using namespace std;

template <class T>
T  **twoDynamicArray (int numberRows, int numberColumns)
{
         T **dynamicArray;
         dynamicArray= new T*[numberRows];
         for( int i=0; i<numberRows; i++)
         dynamicArray[i]= new T [numberColumns];
         return dynamicArray;
}




template <class T>
void deleteDynamicArray(T** dynArray)
{
         delete [] *dynArray;
         delete [] dynArray;
}





int main() 
{

const int sRows=35947;
const int sColumns=5;

float **surfacePoints=twoDynamicArray <float> (sRows, sColumns);

ifstream surFile;
 surFile.open("bunnySurface.txt", ios::in | ios::out);
    if (surFile.is_open())
        {
         for (int r=0; r<sRows; r++)
             {
               for (int c=0; c<sColumns; c++)
                   {
                     surFile>>surfacePoints[r][c];
          //           cout<<surfacePoints[r][c]<<"  ";
                   }
               cout<<endl;
             }
        //  cout<<endl<<surfacePoints[0][2]<<"    "<<surfacePoints[35946][2]<<endl;
          surFile.close();
         }
     else
         {
           cout<<" Couldn't open the bunnySurface.txt";
         }



/*
const int vRows=143469;
const int vColumns=4;

float **volumePoints=twoDynamicArray <float> (vRows, vColumns);

ifstream volFile;
 volFile.open("bunny_2.inp", ios::in | ios::out);
    if (volFile.is_open())
        {
         for (int r=0; r<vRows; r++)
             {
               for (int c=0; c<vColumns; c++)
                   {
                     volFile>>volumePoints[r][c];
                     cout<<volumePoints[r][c]<<"  ";
                   }
               cout<<endl;
             }
         // cout<<endl<<surfacePoints[0][2]<<"    "<<surfacePoints[35946][2]<<endl;
          volFile.close();
         }
     else
         {
           cout<<" Couldn't open the bunny_2.inp file";
         }

*/

// string line;
  // ifstream myfile;
  // myfile.open ("bunny.ply", ios::in);
  // if (myfile.is_open())
 // {
   // cout<<"The bunny.ply file is open now"<<endl;
   //    while (myfile.good())
     //   {
       //   getline (myfile, line);
       //   cout << line << endl;
       // }
       // myfile.close();
  // }
    //  else {cout<<"The bunny.ply file is not open"<<endl;}
  // myfile.close();
  // cout<<endl<<surfacePoint[0][3]<<"    "<<surfacePoint[35947][2];
return 0;
}
