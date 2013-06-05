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


int * surfaceIndexPoint;
surfaceIndexPoint= new int [sRows];
float precs=0.0001;
int n=0;

ofstream theFile ("notSurfacePoints.txt");
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
                                   surfaceIndexPoint[n]=s+1;
//                                   cout<<n<<"    "<<surfaceIndexPoint[n]<<endl;
                                   n=n+1;
                                 }
                }
       }

const int numSurPoints=33059;
int m=1, q=-1;
bool b1, b2, b3;


   for (int k=1; k<=sRows; k++)
       {
//            cout<<"k="<<k<<":"<<endl;
//            p=-1;
//            while ((surfaceIndexPoint[p]==k) || (p>=numSurPoints))
//             do
           for(int p=0; p<numSurPoints; p++)
                 {
//                   if (p>=numSurPoints)
//                      break;

//                   cout<<"   p="<<p<<"    surfaceIndexPoint["<<p<<"]=" <<surfaceIndexPoint[p];
//                   p=p+1;
/*                    cout<<"   p="<<p<<"    surfaceIndexPoint["<<p<<"]=" <<surfaceIndexPoint[p];
                    b1=(surfaceIndexPoint[p]==k);
                    b2=(p>=numSurPoints);
                    b3=((surfaceIndexPoint[p]) or (p>=numSurPoints));
                    cout<<"  (surfaceIndexPoint["<<p<<"]="<<surfaceIndexPoint[p]<<"=="<<k<<")="<<b1;
                    cout<<"  ("<<p<< ">=numSurPoints="<<numSurPoints<<")="<<b2<<endl;
                    cout<<"  ( (surfaceIndexPoint["<<p<<"]="<<surfaceIndexPoint[p]<<"=="<<k<<") || (" <<p<< ">=numSurPoints="<<numSurPoints<<") )="<<b3<<endl;
*/                  if (surfaceIndexPoint[p]==k)
                       {
                         q=p;
                         break;
                       }
                  }

          if (q==-1)
              {
//                 cout<<m<<"    ";
                 theFile<<m<<"    ";
                 for (int c=0; c<sColumns-1; c++)
                    {
//                      cout<< surfacePointsNumber[k-1][c]<< "    ";
                      theFile<< surfacePointsNumber[k-1][c]<< "    ";

                    }
//                 cout<<endl;
                 theFile<<endl;
                 m=m+1;
               }
          q=-1;
       }
  }
 else cout<<" Unable to open the surfacePointsInVolumePoints.txt";


deleteDynamicArray(volumePointsNumber);
deleteDynamicArray(surfacePointsNumber);



return 0;
}
