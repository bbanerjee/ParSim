#include <iostream>
#include <fstream>
#include <limits>
#include <cmath>

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

//std::cout.precision(6);
//std::fixed;


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
                   }
             }
          surFile.close();
         }
     else
         {
           cout<<" Couldn't open the bunnySurface.txt";
         }



 float  **surfacePointsNumber=twoDynamicArray <float> (sRows, sColumns+1);

         for (int r=0; r<sRows; r++)
             {
                  surfacePointsNumber[r][0]=r+1;
               for (int c=1; c<sColumns+1; c++)
                   {
                     surfacePointsNumber[r][c]=surfacePoints[r][c-1];
                   }
             }

deleteDynamicArray(surfacePoints);




const int vRows=142474;
const int vColumns=4;

 float  **volumePoints=twoDynamicArray <float> (vRows, vColumns);

ifstream volFile;
 volFile.open("bunnyVolume2.txt", ios::in | ios::out);
    if (volFile.is_open())
        {
         for (int r=0; r<vRows; r++)
             {
               for (int c=0; c<vColumns; c++)
                   {
                     volFile>>volumePoints[r][c];
                   }
             }
          volFile.close();
         }
     else
         {
           cout<<" Couldn't open the bunnyVolume2.txt file";
         }


 float  **volumePointsNumber=twoDynamicArray <float> (vRows, vColumns+1);

         for (int r=0; r<vRows; r++)
             {
                  volumePointsNumber[r][0]=r+1;
               for (int c=1; c<vColumns+1; c++)
                   {
                     volumePointsNumber[r][c]=volumePoints[r][c-1];
                   }
             }

deleteDynamicArray(volumePoints);



//Finding surface points from all points


float precs=0.0001;
int number=1;

ofstream theFile ("surfacePointsInVolumePoints.txt");
 if (theFile.is_open())
  {
   for (int v=0; v<vRows; v++)
       {
            for (int s=0; s<sRows; s++)
                {
                     if (abs(volumePointsNumber[v][2]-surfacePointsNumber[s][1])<precs)
                         if (abs(volumePointsNumber[v][3]-surfacePointsNumber[s][2])<precs)
                              if (abs(volumePointsNumber[v][4]-surfacePointsNumber[s][3])<precs)
                                 {
                                    cout<<number<<"    ";
                                    number=number+1;
                                    for (int c=0; c<=vColumns; c++)
                                        {
                                          theFile<<volumePointsNumber[v][c]<<"   ";
//                                           cout<<volumePointsNumber[v][c]<<"   ";

                                        }
                                 theFile<<s+1<<endl;
                                 cout<<s+1<<endl;
                                 }
                }
       }
   theFile.close();
  } 
 else cout<<" Unable to open the surfacePointsInVolumePoints.txt";


deleteDynamicArray(volumePointsNumber);
deleteDynamicArray(surfacePointsNumber);



return 0;
}
