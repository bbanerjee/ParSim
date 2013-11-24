#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <map>

#include <boost/variant.hpp>

#include <Vector3D.h>
#include <Matrix3D.h>
#include <CellIndexVector.h>
#include <CellInterpolationVector.h>

using namespace std;
using namespace BrMPM;

typedef boost::variant<std::vector<double>, std::vector<int>, std::vector<Vector3D>,
                       std::vector<Matrix3D>, std::vector<CellIndexVector>, std::vector<CellInterpolationVector> >  variant_type;
typedef std::map<std::string, variant_type> map_type;



 std::ostream&
 operator << (std::ostream& os, variant_type& var) {
   if (std::vector<double>* valDblVec = boost::get<std::vector<double> >(&var)) {
       os << "vector of doubles:" << std::endl << "[";
       for (auto iter = valDblVec->begin(); iter != valDblVec->end(); ++iter)  {
           os << *iter << "  ";
           }
       os << "]";
       } else if (std::vector<int>* valIntVec = boost::get<std::vector<int> >(&var)) {
                 os << "vector of integers:" << std::endl <<  "[";
                 for (auto iter = valIntVec->begin(); iter != valIntVec->end(); ++iter)  {
                     os << *iter << "  ";
                 }
                 os << "]";
       } else if (std::vector<Vector3D>* val3dVec = boost::get<std::vector<Vector3D> >(&var)) {
                 os << "vector of Vector3Ds:" << std::endl <<  "[";
                 for (auto iter = val3dVec->begin(); iter != val3dVec->end(); ++iter)  {
                     Vector3D vec3 = *iter;
                     os << "(" << vec3.x() << ", " << vec3.y() << ", " << vec3.z() << ")," << std::endl;
                 }
                 os << "]";
       } else if (std::vector<Matrix3D>* val3dMat = boost::get<std::vector<Matrix3D> >(&var)) {
                 os << "vector of Matrix3Ds:" << std::endl <<  "[" << std::endl;
                 for (auto iter = val3dMat->begin(); iter != val3dMat->end(); ++iter)  {
                     Matrix3D mat3 = *iter;
                     os << "[";
                     for (int ii = 0; ii < 3; ++ii) {
                         for (int jj = 0; jj < 3; ++jj) {
                             if ((ii != 2) && (jj != 2)) 
                                 os << mat3(ii, jj) << "  ";
                             else if ((ii == 2) && (jj != 2))
                                      os << mat3(ii, jj) << "  ";
                             else if ((ii == 2) && (jj == 2))
                                      os << mat3(ii, jj) << "],"; //  << std::endl;
                             else 
                                      os << mat3(ii, jj) << "  ";
                         }
                         os << std::endl;
                     }
                 }
                 os << "]";
      } else if (std::vector<CellIndexVector>* vecCellIndex = boost::get<std::vector<CellIndexVector> >(&var)) {
                os << "vector of CellIndexVectors:" << std::endl <<  "[" << std::endl;
                for (auto iter = vecCellIndex->begin(); iter != vecCellIndex->end(); ++iter)  {
                    CellIndexVector cellInd = *iter;
                    os << "[";
                    for (unsigned int i = 0; i < cellInd.size(); ++i) {
                         os << cellInd[i] << "  ";
                    }
                    os << "]," << std::endl;
                 }
                 os << "]";
      }  else if (std::vector<CellInterpolationVector>* vecCellInterpol = boost::get<std::vector<CellInterpolationVector> >(&var))                                                                                                                   {
                os << "vector of CellInterpolationVectors:" << std::endl <<  "[" << std::endl;
                for (auto iter = vecCellInterpol->begin(); iter != vecCellInterpol->end(); ++iter)  {
                    CellInterpolationVector cellInter = *iter;
                    os << "[";
                    for (unsigned int i = 0; i < cellInter.size(); ++i) {
                         os << cellInter[i] << "  ";
                    }
                    os << "]," << std::endl;
                 }
                 os << "]";            
     }
     return os;
 }


 


 template<typename T>
 void add(map_type& map, const std::string& name_map, const std::string& lable, int dwi, T& val)
      {
        std::string lable_dwi = lable + std::to_string(dwi);
        variant_type var(val);
        map[lable_dwi] = var;
       std::cout << name_map << "[" << lable_dwi << "] = " << map[lable_dwi] << std::endl;
      } 



 template<typename T>
 void get(map_type& map, const std::string& lable, int dwi, T& val)
      {
        std::string lable_dwi = lable + std::to_string(dwi);
        val = boost::get<T>(map[lable_dwi]);
        variant_type var(val);
//        std::cout << var << std::endl;
      }




/* template<typename T, typename S>
 void put(map_type& map, const std::string& lable, int dwi, const T& val, S& variant)
      {
        std::string lable_dwi = lable + std::to_string(dwi);
         variant = map[lable_dwi];        
          boost::get<T>(variant) = val;     
          map[lable_dwi] = variant;
      } */


template<typename T>
 void put(map_type& map, const std::string& lable, int dwi, const T& val)
      {
        std::string lable_dwi = lable + std::to_string(dwi);
         variant_type variant;
         variant = map[lable_dwi];
    //    variant_type var(val);
     //   T& temp  = boost::get<T>(map[lable_dwi]);
   //     temp = var;
   //     std::cout << var << std::endl;
   //     std::cout << map[lable_dwi] << std::endl;
        
          boost::get<T>(variant) = val;     
          map[lable_dwi] = variant;
      } 


 class PrintVisitor : public boost::static_visitor<std::string>
 {
 public:
    PrintVisitor(std::ostringstream& os) : d_os(os) {}    

    std::string operator () (std::vector<double>& val) const {
        d_os << "visitor vector of doubles:" << std::endl << "["; 
        for (auto iter = val.begin(); iter != val.end(); ++iter) {
             d_os << *iter << "  ";
        }
        d_os << "]";     
        return d_os.str();
    }

     std::string operator () (std::vector<int>& val) const {
        d_os << "visitor vector of integers:" << std::endl << "["; 
        for (auto iter = val.begin(); iter != val.end(); ++iter) {
             d_os << *iter << "  ";
        }
        d_os << "]";     
        return d_os.str();
    }
    
    std::string operator () (std::vector<Vector3D>& val) const {
        d_os << "visitor vector of Vector3Ds:" << std::endl <<  "[";
        for (auto iter = val.begin(); iter != val.end(); ++iter) {
            Vector3D vec3 = *iter;
            d_os << "(" << vec3.x() << ", " << vec3.y() << ", " << vec3.z() << ")," << std::endl;
        }
        d_os << "]";
        return d_os.str();
    }

    std::string operator () (std::vector<Matrix3D>& val) const {
        d_os << "visitor vector of Matrix3Ds:" << std::endl <<  "[" << std::endl;
        for (auto iter = val.begin(); iter != val.end(); ++iter) {
             Matrix3D mat3 = *iter;
             d_os << "[";
             for (int ii = 0; ii < 3; ++ii) {
                 for (int jj = 0; jj < 3; ++jj) {
                     if ((ii != 2) && (jj != 2)) 
                         d_os << mat3(ii, jj) << "  ";
                     else if ((ii == 2) && (jj != 2))
                             d_os << mat3(ii, jj) << "  ";
                     else if ((ii == 2) && (jj == 2))
                             d_os << mat3(ii, jj) << "],"; //  << std::endl;
                     else 
                             d_os << mat3(ii, jj) << "  ";
                 }
                 d_os << std::endl;
             }
        }
        d_os << "]";
        return d_os.str();
    }

   std::string operator () (std::vector<CellIndexVector>& val) const {
        d_os << "visitor vector of CellIndexVectors:" << std::endl <<  "[" << std::endl;
        for (auto iter = val.begin(); iter != val.end(); ++iter) {
            CellIndexVector cellInd = *iter;
            d_os << "[";
            for (unsigned int i = 0; i < cellInd.size(); ++i) {
                d_os << cellInd[i] << "  ";
            }
            d_os << "]," << std::endl;
         }
         d_os << "]";
         return d_os.str();            
    }


   std::string operator () (std::vector<CellInterpolationVector>& val) const {
        d_os << "visitor vector of CellInterpolationVectors:" << std::endl <<  "[" << std::endl;
        for (auto iter = val.begin(); iter != val.end(); ++iter) {
            CellInterpolationVector cellInter = *iter;
            d_os << "[";
            for (unsigned int i = 0; i < cellInter.size(); ++i) {
                d_os << cellInter[i] << "  ";
            }
            d_os << "]," << std::endl;
        }
        d_os << "]";
        return d_os.str();            
   }
    
    
 private:
 std::ostringstream& d_os;
   
};






 class ZeroVisitor : public boost::static_visitor<void>
 {
 public:
    void operator () (std::vector<double>& val) const {
        for (auto iter = val.begin(); iter != val.end(); ++iter) {
             *iter = 0.0;
        }
    }

    void operator () (std::vector<int>& val) const {
        for (auto iter = val.begin(); iter != val.end(); ++iter) {
             *iter = 0;
        }
    }

    void operator () (std::vector<Vector3D>& val) const {
        for (auto iter = val.begin(); iter != val.end(); ++iter) {
             Vector3D vec3 = *iter;
             vec3.x(0.0);
             vec3.y(0.0);
             vec3.z(0.0);
             *iter = vec3;         
        }
    }

    void operator () (std::vector<Matrix3D>& val) const {
        for (auto iter = val.begin(); iter != val.end(); ++iter) {
             Matrix3D mat3 = *iter;
             mat3.set(0.0);
             *iter = mat3;         
        }
    }

   void operator () (std::vector<CellIndexVector>& val) const {
        for (auto iter = val.begin(); iter != val.end(); ++iter) {
             CellIndexVector cellIndex = *iter;
             int size = cellIndex.size();
             std::vector<int> zeroVec(size, 0);
             cellIndex = zeroVec;
             *iter = cellIndex;         
        }
    }


   void operator () (std::vector<CellInterpolationVector>& val) const {
        for (auto iter = val.begin(); iter != val.end(); ++iter) {
             CellInterpolationVector cellInter = *iter;
             int size = cellInter.size();
             std::vector<double> zeroVector(size, 0.0);
             cellInter = zeroVector;
             *iter = cellInter;         
        }
    }
   
};

    
            


  

int main ()
 {
   std::vector<double> vec(9, 0.1);

   vec.push_back(10.5);
   vec.emplace_back(19.5);
   vec[4] = 100;

   std::vector<int> intVec(5, 2);
   intVec.push_back(4);


  const int number = 10; 
  std::vector<Vector3D> vec3D(number, Vector3D(5.2, 4.9, 6.3));
  vec3D.emplace_back(0.0, 0.0, 0.0);
  vec3D.emplace_back(0.0, 0.0, 0.0);

 variant_type v1(vec);
 variant_type v2(intVec);
 variant_type v3(vec3D);

 map_type map1;
 map1["first"] = v1;
 map1["second"] = v2;
 map1["third"] = v3;

   std::cout
   << "map1[first] = " << map1["first"] << std::endl 
   << "map1[second] = " << map1["second"] << std::endl 
   << "map1[third] = " << map1["third"] << std::endl;

 const std::string mapName = "map1";

 std::vector<Vector3D> pw(12, Vector3D(0.0));
 add(map1, mapName, "pw", 1, pw);

 std::vector<double> zeroVecDouble (7, 0.56);
 add(map1, mapName, "cIdx", 2, zeroVecDouble);

 std::vector<int> zeroVecInt (10, 7);
 add(map1, mapName, "cW", 1, zeroVecInt);

 std::vector<Matrix3D> matrix3 (5, Matrix3D(0.0));
 add(map1, mapName, "pGv", 2, matrix3);

 std::vector<int> exam_int(6, 2);
 std::vector<CellIndexVector> vec_cell (4, exam_int);
 add(map1, mapName, "cellIdx", 1, vec_cell);

 std::vector<double> exam_inter(8, 0.11);
 std::vector<CellInterpolationVector> vec_inter (6, exam_inter);
 add(map1, mapName, "cellInterpolation", 1, vec_inter);



//zero the elements and print 

for (auto iter = map1.begin(); iter != map1.end(); ++iter) {
    boost::apply_visitor(ZeroVisitor(), iter->second);
    std::ostringstream os_map1;
    PrintVisitor pv_map1(os_map1);
    std::cout << mapName << "[" << iter->first << "]: " << boost::apply_visitor(pv_map1, iter->second) << std::endl;
}
std::cout << std::endl;    
     

  Matrix3D mat;
  mat.Identity();
  std::vector<Matrix3D> vecOfMatrix(3, mat);
  
  add(map1, mapName, "pGv", 2, vecOfMatrix);
 
  std::vector<CellInterpolationVector> vecOfInter;
 get(map1, "cellInterpolation", 1, vecOfInter);
 variant_type exam(vecOfInter);
 std::cout << exam << std::endl;

 std::ostringstream os_1;
 PrintVisitor pv_1(os_1);
 std::cout << boost::apply_visitor(pv_1, exam);
 
 variant_type exam2(vecOfMatrix);
   std::cout << exam2 << std::endl;
// put(map1, "pGv", 2, exam2);
// std::cout << exam2 << std::endl;
// std::cout << map1["pGv2"] << std::endl;


 std::vector<Matrix3D> vectorMatrix(5, mat);
 boost::get<std::vector<Matrix3D> >(map1["pGv2"]) = vectorMatrix;
 std::cout << map1["pGv2"] << std::endl;
 
 std::vector<Matrix3D> vecMat(3, Matrix3D(0.0)); 

 //put(map1, "pGv", 2, vecMat, exam);
 put(map1, "pGv", 2, vecMat);
 std::cout << map1["pGv2"] << std::endl;
 
 }
       


