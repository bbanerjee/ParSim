// Function Definitions

#include <Core/Math/Matrix.h>
#include <cmath>
#include <cstdlib>
// written on Feb 13, 2013
namespace dem {

Matrix::Matrix(int i, int j)
{
  if (i < 0 || j < 0) {
    std::cerr << "Error: the dimensions of Matrix should be positive!"
              << std::endl;
    exit(-1);
  } else {
    num_row = i;
    num_col = j;

    value = new REAL[num_row * num_col];
    for (int ii = 0; ii < num_row * num_col; ii++) {
      value[ii] = 0;
    }
  }
}

Matrix::Matrix(const Matrix& A)
{

  num_row = A.num_row;
  num_col = A.num_col;

  value = new REAL[num_row * num_col];
  for (int i = 0; i < num_row * num_col; i++) {
    value[i] = A.value[i];
  }
}

Matrix::~Matrix()
{
  clear();
}

/*
std::vector<REAL> Matrix::getCol(int i) const{
        if(i<=0 || i>num_col){
                //std::cout << "Error: index exceeds!"	<< std::endl;
                exit(-1);
        }
        std::vector<REAL> result;
        result.clear();
        for(std::vector<std::vector<REAL> >::const_iterator itr=value.begin();
itr!=value.end(); itr++){
                std::vector<REAL>::const_iterator itc=(*itr).begin();
                for(int ic=0; ic!=i-1; ic++){
                        itc++;		// move to the ith element of itr row
                }
                result.push_back(*itc);
        }
        return result;
}

std::vector<REAL> Matrix::getRow(int i) const{
        if(i<=0 || i>num_row){
                //std::cout << "Error: index exceeds!"	<< std::endl;
                exit(-1);
        }
//	std::vector<std::vector<REAL> >::const_iterator itr=value.begin();
        std::vector<std::vector<REAL> >::size_type itr = 0;
        for(int ir=0; ir!=i-1; ir++){
                itr++;
        }
        return value[itr];
}


void Matrix::calcDimensions(){

    std::vector<std::vector<REAL> >::size_type row = value.size();
    std::vector<REAL>::size_type col = value[0].size();

    num_row = row;
    num_col = col;

} // calcDimensions()


void Matrix::appendRow(const std::vector<REAL> & row){
        if(row.size() != num_col && num_col != 0){	// num_col != 0 in
case that the Matrix is empty
                //std::cout << "Error: the dimesions do not match!" << std::endl;
                exit(-1);
        }
        value.push_back(row);
        num_row++;
        num_col = row.size();	// this is important to avoid the case that
at first I append a 3 elements row

        calcDimensions();
                                        // while num_col = 0, num_col should be
3
}

void Matrix::appendCol(const std::vector<REAL> & column){
        if(column.size() != num_row && num_row != 0){
                //std::cout << "Error: the dimesions do not match!" << std::endl;
                exit(-1);
        }
        std::vector<REAL> temp_row;
        if(num_row == 0){
                for(std::vector<REAL>::const_iterator itc=column.begin();
itc!=column.end(); itc++){
                        temp_row.clear();
                        temp_row.push_back(*itc);
                        this->appendRow(temp_row);
                }
        }
        else {
                int ir = 0;
                for(std::vector<std::vector<REAL> >::iterator itr=value.begin();
itr!=value.end(); itr++){	// if num_row == 0, then
means this Matrix is empty and this for loop will not be getted in which leads
to cannot appendCol to empty Matrix
                        (*itr).push_back(column[ir]);	// put the ir
element of col to the back of itr row
                        if(ir>=num_row){
                                //std::cout << "Error: dimension exceeds in
appendCol!" << std::endl;
                                exit(-1);
                        }
                        ir++;
                }
                num_col++;
        }
        num_row = column.size();
        calcDimensions();
}
*/

void
Matrix::clear()
{

  if (value != nullptr) { // not a 0x0 matrix
    delete[] value;
    num_col = 0;
    num_row = 0;

    value = nullptr;
  }
}

/*
void Matrix::getInvs(){	// at present this code can only get the inverse
of a two by two Matrix and 1x1 Matrix
        if(num_row == 2 && num_col ==2){
                REAL det;
                REAL a, b, c, d;		// [a b; c d]
                a = value[0];
                b = value[1];
                c = value[2];
                d = value[3];
                det = a*d-b*c;
//		if(det == 0){
//			//std::cout << "Error: Matrix is singular!" << std::endl;
//			exit(-1);
//		}
                value[0] = d/det;
                value[1] = -b/det;
                value[2] = -c/det;
                value[3] = a/det;

        }
        else if(num_row == 1 && num_col == 1){
                value[0] = 1.0/value[0];

        }
        else if(num_row == 3 && num_col == 3){
                double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13;
                double z14, z15, z16, z17, z18, z19, z20, z21, z22, z23, z24,
z25;

                z1 = value[0];
                z2 = value[1];
                z3 = value[2];
                z4 = value[3];
                z5 = value[4];
                z6 = value[5];
                z7 = value[6];
                z8 = value[7];
                z9 = value[8];
                z10 = -z2*z4;
                z11 = z3*z4;
                z12 = z1*z5;
                z13 = -z3*z5;
                z14 = -z1*z6;
                z15 = z2*z6;
                z16 = z13*z7;
                z17 = z15*z7;
                z18 = z2*z7;
                z19 = -z3*z7;
                z20 = -z5*z7;
                z7 = z6*z7;
                z21 = -z1*z8;
                z22 = z11*z8;
                z23 = z14*z8;
                z3 = z3*z8;
                z24 = z4*z8;
                z6 = -z6*z8;
                z1 = z1*z9;
                z8 = z10*z9;
                z25 = z12*z9;
                z2 = -z2*z9;
                z4 = -z4*z9;
                z5 = z5*z9;
                z9 = z10 + z12;
                z10 = z11 + z14;
                z11 = z13 + z15;
                z12 = z18 + z21;
                z13 = z20 + z24;
                z1 = z1 + z19;
                z8 = z16 + z17 + z22 + z23 + z25 + z8;
                z2 = z2 + z3;
                z3 = z4 + z7;
                z4 = z5 + z6;
                z5 = 1.0/z8;
                z6 = z5*z9;
                z7 = z10*z5;
                z8 = z11*z5;
                z9 = z12*z5;
                z10 = z13*z5;
                z1 = z1*z5;
                z2 = z2*z5;
                z3 = z3*z5;
                z4 = z4*z5;

                //{{z4, z2, z8},
                // {z3, z1, z7},
                // {z10, z9, z6}}


                value[0] = z4;
                value[1] = z2;
                value[2] = z8;
                value[3] = z3;
                value[4] = z1;
                value[5] = z7;
                value[6] = z10;
                value[7] = z9;
                value[8] = z6;

////std::cout << "z4: " << z4 << std::endl;
////std::cout << "Matrix.cpp, result: " << print() << std::endl;

        }
        else {
                //std::cout << "Sorry: at present this code can only get the
inverse of a 2 by 2 Matrix and 3x3 Matrix!" << std::endl;
                exit(-1);
        }

}


void Matrix::getTrans() {
        Matrix result(num_col, num_row);

        for(int i=0; i<num_row; i++){
            for(int j=0; j<num_col; j++){
                result.value[j*num_row+i] = value[i*num_col+j];
            }
        }

        num_row = result.num_row;
        num_col = result.num_col;
        for(int i=0; i<num_row*num_col; i++){
            value[i] = result.value[i];
        }

}
*/

REAL
Matrix::getNorm()
{
  if (num_row == 1) { // for column vector
    REAL norm = 0;
    for (int ic = 0; ic != num_col; ic++) {
      norm = norm + value[ic] * value[ic];
    }
    norm = sqrt(norm);
    return norm;
  } else if (num_col == 1) { // for row vector
    REAL norm = 0;
    for (int ir = 0; ir != num_row; ir++) {
      norm = norm + value[ir] * value[ir];
    }
    norm = sqrt(norm);
    return norm;
  } else {
    std::cerr << "Sorry: at present we can only get the norm of vectors!"
              << std::endl;
    exit(-1);
  }
}

void
Matrix::LU(Matrix& L, Matrix& U) const
{

  if (num_col != num_row) {
    std::cerr << "should be square Matrix in LU decomposion..." << std::endl;
    exit(-1);
  }

  L = zeros(num_row, num_col);
  U = zeros(num_row, num_col);

  //#pragma omp parallel shared(a, L, U)
  //{
  //    #pragma omp for schedule(static, 10)
  for (int k = 0; k < num_row; ++k) {
    L.value[k * num_col + k] = 1;
    for (int i = k; i < num_row; i++) {
      L.value[i * num_col + k] =
        value[i * num_col + k] / value[k * num_col + k];
      // value[i*num_col+k] = L.value[i*num_col+k];
      for (int j = k; j < num_row; j++) {
        value[i * num_col + j] =
          value[i * num_col + j] -
          L.value[i * num_col + k] * value[k * num_col + j];
      }
      //	    #pragma omp flush(a)
    }
    //	#pragma omp nowait
    for (int j = k - 1; j < num_row; j++) {
      U.value[k * num_col + j] = value[k * num_col + j];
    }
  }

  //} // end parallel

} // LU()

REAL
Matrix::getVal(int i, int j) const
{

  if (i > num_row || i < 1 || j > num_col || j < 1) {
    std::cerr << "Error: index exceeds when ()!" << std::endl;
    std::cerr << "i,j in ():\n " << i << ", " << j << std::endl;
    std::cerr << "num_row, num_col:\n " << num_row << " " << num_col
              << std::endl;
    exit(-1);
  }
  return value[(i - 1) * num_col + j - 1];

} // getVal()

Matrix&
Matrix::operator=(const Matrix& A)
{
  // need to delete Matrix *this first
  this->clear();
  this->value = new REAL[A.num_row * A.num_col];
  num_row = A.num_row;
  num_col = A.num_col;
  for (int i = 0; i < num_row * num_col; i++) {
    value[i] = A.value[i];
  }

  return *this;
}

REAL&
Matrix::operator()(int i, int j)
{
  if (i > num_row || i < 1 || j > num_col || j < 1) {
    std::cerr << "Error: index exceeds when ()!" << std::endl;
    std::cerr << "i,j in ():\n " << i << ", " << j << std::endl;
    std::cerr << "num_row, num_col:\n " << num_row << " " << num_col
              << std::endl;
    exit(-1);
  }

  return value[(i - 1) * num_col + j - 1];
}

std::string
Matrix::print() const
{
  //	//std::cout << "Matrix: " << std::endl;
  std::stringstream ss;
  for (int i = 0; i < num_row; i++) {
    for (int j = 0; j < num_col; j++) {
      ss << value[i * num_col + j] << " ";
    }
    ss << std::endl;
  }
  return ss.str();
}

// void operator += (Matrix A){
//	(*this) = (*this)+A;
//}

// non member functions
Matrix
operator+(const Matrix& A, const Matrix& B)
{
  if (A.num_row != B.num_row || A.num_col != B.num_col) {
    //std::cout << "Error: dimensions do not match!" << std::endl;
    exit(-1);
  }
  Matrix result(A.num_row, A.num_col);
  for (int i = 0; i < (A.num_row) * (A.num_col); i++) {
    result.value[i] = A.value[i] + B.value[i];
  }

  return result;
}

Matrix
operator+(REAL k, const Matrix& A)
{
  Matrix result(A.num_row, A.num_col);
  for (int i = 0; i < (A.num_row) * (A.num_col); i++) {
    result.value[i] = A.value[i] + k;
  }

  return result;
}

Matrix
operator+(const Matrix& A, REAL k)
{
  Matrix result(A.num_row, A.num_col);
  for (int i = 0; i < (A.num_row) * (A.num_col); i++) {
    result.value[i] = A.value[i] + k;
  }

  return result;
}

Matrix
operator-(const Matrix& A, const Matrix& B)
{
  if (A.num_row != B.num_row || A.num_col != B.num_col) {
    std::cerr << "Error: dimensions do not match!" << std::endl;
    exit(-1);
  }
  Matrix result(A.num_row, A.num_col);
  for (int i = 0; i < (A.num_row) * (A.num_col); i++) {
    result.value[i] = A.value[i] - B.value[i];
  }

  return result;
}

Matrix
operator-(REAL k, const Matrix& A)
{
  Matrix result(A.num_row, A.num_col);
  for (int i = 0; i < (A.num_row) * (A.num_col); i++) {
    result.value[i] = k - A.value[i];
  }

  return result;
}

Matrix
operator-(const Matrix& A, REAL k)
{
  Matrix result(A.num_row, A.num_col);
  for (int i = 0; i < (A.num_row) * (A.num_col); i++) {
    result.value[i] = A.value[i] - k;
  }

  return result;
}

Matrix operator*(const Matrix& A, const Matrix& B)
{
  if (A.num_col != B.num_row) {
    std::cerr << "Error: dimensions do not match!" << std::endl;
    exit(-1);
  }

  int M = A.num_row;
  int N = B.num_row;
  int K = B.num_col;
  Matrix result(M, K);

  for (int ii = 0; ii < M; ii++) {
    for (int jj = 0; jj < N; jj++) {
      for (int kk = 0; kk < K; kk++) {
        result.value[ii * K + kk] +=
          A.value[ii * N + jj] * B.value[jj * K + kk];
      }
    }
  }

  return result;
}

Matrix operator*(const Matrix& A, REAL k)
{
  Matrix result(A.num_row, A.num_col);
  for (int i = 0; i < (A.num_row) * (A.num_col); i++) {
    result.value[i] = A.value[i] * k;
  }

  return result;
}

Matrix operator*(REAL k, const Matrix& A)
{
  Matrix result(A.num_row, A.num_col);
  for (int i = 0; i < (A.num_row) * (A.num_col); i++) {
    result.value[i] = A.value[i] * k;
  }

  return result;
}

Matrix
operator/(const Matrix& A, REAL k)
{
  Matrix result(A.num_row, A.num_col);
  for (int i = 0; i < (A.num_row) * (A.num_col); i++) {
    result.value[i] = A.value[i] / k;
  }

  return result;
}

Matrix
operator%(Matrix& A, Matrix& r)
{ // left division, "\"

  Matrix x;
  MatrixEqnSolver(x, A, r);

  return x;

} // %

Matrix
expm(const Matrix& A)
{
  Matrix result(A.num_row, A.num_col);
  for (int i = 0; i < (A.num_row) * (A.num_col); i++) {
    result.value[i] = exp(A.value[i]);
  }

  return result;
}

// bool isnan_mat(Matrix &mat){
//
//    bool is = false;
//    for(int i=1; i<mat.num_row+1; ++i){
//	for(int j=1; j<mat.num_col+1; ++j){
//	    if( std::isnan(mat(i,j)) )
//		is = true;
//	}
//    }
//
//    return is;
//
//} // isnan_mat()

int
size(const Matrix& mat, int n)
{
  if (n != 1 && n != 2) {
    std::cerr
      << "Only can get number of rows and number of columns with size()."
      << std::endl;
    exit(-1);
  } else if (n == 1) {
    return mat.num_row;
  } else {
    return mat.num_col;
  }
} // size()

Matrix
ones(int i, int j)
{

  Matrix mat;

  if (i < 0 || j < 0) {
    std::cerr << "Error: the dimensions of Matrix should be positive!"
              << std::endl;
    exit(-1);
  } else if (i == 0 || j == 0) {
    return mat;
  } else {
    mat.num_row = i;
    mat.num_col = j;

    for (int ii = 0; ii < i * j; ii++) {
      mat.value[ii] = 1;
    }
  }

  return mat;

} // ones()

Matrix
zeros(int i, int j)
{

  Matrix mat(i, j);
  return mat;

} // zeros()

REAL
max(Matrix& mat)
{

  REAL max = mat.value[0];
  for (int i = 1; i < (mat.num_row) * (mat.num_col); i++) {
    if (mat.value[i] > max) {
      max = mat.value[i];
    }
  }

  return max;

} // max()

REAL
min(Matrix& mat)
{

  REAL min = mat.value[0];
  for (int i = 1; i < (mat.num_row) * (mat.num_col); i++) {
    if (mat.value[i] < min) {
      min = mat.value[i];
    }
  }

  return min;

} // min()

// Matrix abs(Matrix& A){
//
//    Matrix result(A.num_row, A.num_col);
//    for(int i=0; i<(A.num_row)*(A.num_col); i++){
//	result.value[i] = abs(A.value[i]);
//    }
//
//    return result;
//
//} // abs()

int
length(const Matrix& mat)
{

  int num;
  if (mat.num_row == 0 || mat.num_col == 0) {
    num = 0;
  } else {
    num = mat.num_row;
    if (num < mat.num_col)
      num = mat.num_col;
  }
  return num;

} // length()

REAL
norm(Matrix& mat)
{

  REAL val = mat.getNorm();

  return val;

} // norm()

REAL
det(Matrix& mat)
{

  if (mat.num_row == 2 && mat.num_col == 2) { // (2 x 2)

    return mat.value[0] * mat.value[3] - mat.value[1] * mat.value[2];
  } else if (mat.num_row == 1 && mat.num_col == 1) { // (1 x 1)

    return mat.value[0];
  } else if (mat.num_row == 3 && mat.num_col == 3) { // (3 x 3)

    return mat.value[0] *
             (mat.value[4] * mat.value[8] - mat.value[5] * mat.value[7]) -
           mat.value[1] *
             (mat.value[3] * mat.value[8] - mat.value[5] * mat.value[6]) +
           mat.value[2] *
             (mat.value[3] * mat.value[7] - mat.value[4] * mat.value[6]);

  } else {
    std::cerr << "Matrix dimensions problem in det()..." << std::endl;
    exit(-1);
  }

} // det()

REAL
trace(Matrix& mat)
{

  if (mat.num_row != mat.num_col) {
    std::cerr << "The matrix is not square..." << std::endl;
    exit(-1);
  }

  REAL res = 0;
  for (int i = 1; i != mat.num_row + 1; i++) {
    res += mat(i, i);
  }

  return res;

} // trace()

Matrix
linspace(REAL start, REAL stop, int num)
{

  Matrix result(1, num);

  if (num <= 0) {
    std::cerr << "Number should be positive in linspace()..." << std::endl;
    exit(-1);
  } else if (num == 1) {
    result.value[0] = stop;
  } else {

    REAL step = (stop - start) / double((num - 1));
    for (int i = 1; i < num - 1; ++i) {
      result.value[i] = start + step * (i - 1);
    }
    result.value[num - 1] = stop;
  }

  return result;

} // linspace()

void
MatrixEqnSolver(Matrix& x, Matrix& A, Matrix& r)
{

  if (A.num_row != r.num_row || r.num_col != 1) {
    std::cerr << "Dimensions do not match in MatrixEqnSolver()!" << std::endl;
    exit(-1);
  }

  Matrix L, U;
  A.LU(L, U);

  Matrix y(A.num_row, 1);
  for (int i = 1; i < A.num_row + 1; ++i) {
    REAL left = 0;
    for (int j = 1; j < i; ++j) {
      left += L(i, j) * y(j, 1);
    }
    y(i, 1) = (r(i, 1) - left) / L(i, i);
  }

  x = zeros(r.num_row, 1);

  for (int i = A.num_row; i > 0; --i) {
    REAL left = 0;
    for (int j = i + 1; j < A.num_col + 1; ++j) {
      left += U(i, j) * x(j, 1);
    }
    x(i, 1) = (y(i, 1) - left) / U(i, i);
  }

} // MatrixEqnSolver()

Matrix
inv(const Matrix& A)
{

  if (A.num_row == 2 && A.num_col == 2) {
    Matrix res(2, 2);
    REAL det;
    REAL a, b, c, d; // [a b; c d]
    a = A.value[0];
    b = A.value[1];
    c = A.value[2];
    d = A.value[3];
    det = a * d - b * c;
    if (det == 0) {
      std::cerr << "Error: Matrix is singular!" << std::endl;
      exit(-1);
    }
    res.value[0] = d / det;
    res.value[1] = -b / det;
    res.value[2] = -c / det;
    res.value[3] = a / det;

    return res;

  } else if (A.num_row == 1 && A.num_col == 1) {
    Matrix res(1, 1);
    res.value[0] = 1.0 / A.value[0];

    return res;

  } else if (A.num_row == 3 && A.num_col == 3) {
    double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13;
    double z14, z15, z16, z17, z18, z19, z20, z21, z22, z23, z24, z25;

    z1 = A.value[0];
    z2 = A.value[1];
    z3 = A.value[2];
    z4 = A.value[3];
    z5 = A.value[4];
    z6 = A.value[5];
    z7 = A.value[6];
    z8 = A.value[7];
    z9 = A.value[8];
    z10 = -z2 * z4;
    z11 = z3 * z4;
    z12 = z1 * z5;
    z13 = -z3 * z5;
    z14 = -z1 * z6;
    z15 = z2 * z6;
    z16 = z13 * z7;
    z17 = z15 * z7;
    z18 = z2 * z7;
    z19 = -z3 * z7;
    z20 = -z5 * z7;
    z7 = z6 * z7;
    z21 = -z1 * z8;
    z22 = z11 * z8;
    z23 = z14 * z8;
    z3 = z3 * z8;
    z24 = z4 * z8;
    z6 = -z6 * z8;
    z1 = z1 * z9;
    z8 = z10 * z9;
    z25 = z12 * z9;
    z2 = -z2 * z9;
    z4 = -z4 * z9;
    z5 = z5 * z9;
    z9 = z10 + z12;
    z10 = z11 + z14;
    z11 = z13 + z15;
    z12 = z18 + z21;
    z13 = z20 + z24;
    z1 = z1 + z19;
    z8 = z16 + z17 + z22 + z23 + z25 + z8;
    z2 = z2 + z3;
    z3 = z4 + z7;
    z4 = z5 + z6;
    z5 = 1.0 / z8;
    z6 = z5 * z9;
    z7 = z10 * z5;
    z8 = z11 * z5;
    z9 = z12 * z5;
    z10 = z13 * z5;
    z1 = z1 * z5;
    z2 = z2 * z5;
    z3 = z3 * z5;
    z4 = z4 * z5;

    //{{z4, z2, z8},
    // {z3, z1, z7},
    // {z10, z9, z6}}

    Matrix res(3, 3);

    res.value[0] = z4;
    res.value[1] = z2;
    res.value[2] = z8;
    res.value[3] = z3;
    res.value[4] = z1;
    res.value[5] = z7;
    res.value[6] = z10;
    res.value[7] = z9;
    res.value[8] = z6;

    return res;

  } else {
    std::cerr << "Sorry: at present this code can only get the inverse of a 2 "
                 "by 2 Matrix and 3x3 Matrix!"
              << std::endl;
    exit(-1);
  }

} // end inv

// returns the transpose of matrix
Matrix
trans(Matrix& a)
{

  int arows = a.num_row;
  int acols = a.num_col;
  Matrix res(acols, arows);

  for (int r = 1; r < arows + 1; r++) {
    for (int c = 1; c < acols + 1; c++) {
      res(c, r) = a(r, c);
    }
  }

  return res;
} // end trans

} // end of dem

/*

// used to test
using namespace periDynamics;

int main(){
        Matrix A;
        A = zeros(3,3);
        //std::cout << "A=zeros(3,3): " << A.print();

        A(1,1) = 3; A(1,2) = -4; A(1,3) = 0;
        A(2,1) = -4; A(2,2) = 2; A(2,3) = 1;
        A(3,1) = 0; A(3,2) = 1; A(3,3) = 1;

        //std::cout << "should be 3 -4 0; -4 2 1; 0 1 1" << std::endl;
        //std::cout << "dimesions of A: " << A.num_row << " " << A.num_col <<
std::endl;
        //std::cout << A.print();

        // test +
        //std::cout << "test + begin: " <<std::endl;
        Matrix Aplus;
        //std::cout << "dimesions of A: " << A.num_row << " " << A.num_col <<
std::endl;
        //std::cout << "dimesions of B: " << A.num_row << " " << A.num_col <<
std::endl;
        Aplus = A + A;
        //std::cout << "should be 6 -8 0; -8 4 2; 0 2 2" << std::endl;
        //std::cout << Aplus.print();
        // test -
        //std::cout << "test - begin: " <<std::endl;
        Matrix Aplusminus;
        Aplusminus = Aplus - A;
        //std::cout << "should be 3 -4 0; -4 2 1; 0 1 1" << std::endl;
        //std::cout << Aplusminus.print();
        // test *
        //std::cout << "test * begin: " <<std::endl;
        Matrix AA;
        AA = A*A;
        //std::cout << "A*A: " << std::endl;
        //std::cout << AA.print();
        Matrix AAA;
        AAA = AA*A;
        //std::cout << "(A*A)*A: " << std::endl;
        //std::cout << AAA.print();

        AAA = A*A*A;
        //std::cout << "A*A*A: " << std::endl;
        //std::cout << AAA.print();
        // test invs
        //std::cout << "test invs begin: " <<std::endl;

        Matrix B(2,2);
        //std::cout << "B(2,2): " << B.print();

        B(1,1) = 8; B(1,2) = -3;
        B(2,1) = -4; B(2,2) = 2;

        //std::cout << "should be 8 -3; -4 2: " << std::endl;
        //std::cout << B.print();
        Matrix Binv;


        Binv = inv(B);
        //std::cout << "Binvs: " <<std::endl;
        //std::cout << Binv.print();
//	//std::cout << "print one by one: " << Binv(1,1) << " " << Binv(1,2) <<
"; " << Binv(2,1) << " " << Binv(2,2) << std::endl;

        //std::cout << "dimesions of Binv: " << Binv.num_row << " " <<
Binv.num_col << std::endl;
        Matrix BBinv;
        BBinv = Binv*B;
        //std::cout << "BBinvs: " << std::endl;
        //std::cout << BBinv.print();
        //std::cout << "dimesions of BBinv: " << BBinv.num_row << " " <<
BBinv.num_col << std::endl;

        // test transpose
        //std::cout << "test transpose begin: " <<std::endl;
        A(1,2) = 0; A(1,3) = 12;
        //std::cout << "A: " << std::endl;
        //std::cout << A.print();
        //std::cout << "dimesions of A: " << A.num_row << " " << A.num_col <<
std::endl;
        Matrix Atran;
        Atran = trans(A);
        //std::cout << "Atran: " << std::endl;
        //std::cout << Atran.print();
        //std::cout << "dimesions of Atran: " << Atran.num_row << " " <<
Atran.num_col << std::endl;
        // test () to see if the elements can be modified
//std::cout << "point 1!" << std::endl;
        Matrix G(4,1);
//	temp_row.clear();
//	temp_row.push_back(1);
//	G.appendRow(temp_row);
//	G.appendRow(temp_row);
//	G.appendRow(temp_row);
//	G.appendRow(temp_row);

        //std::cout << "G: " << std::endl;
        //std::cout << G.print();
        //std::cout << "dimesions of G: " << G.num_row << " " << G.num_col <<
std::endl;
        G(2,1) = 1;
        //std::cout << "G: " << std::endl;
        //std::cout << G.print();
        //std::cout << "dimesions of G: " << G.num_row << " " << G.num_col <<
std::endl;
        //std::cout << "2*G: " << std::endl;
        //std::cout << (2*G).print();
        // test getNorm()
        //std::cout << "G norm: " << G.getNorm() << std::endl;
        // test exp
        //std::cout << "2G exp: " <<  std::endl;
        //std::cout << expm(G).print();
        // tesp -/+
        //std::cout << "G+1: " << std::endl;
        //std::cout << (G+1).print();
        //std::cout << "1+G: " << std::endl;
        //std::cout << (1+G).print();
        //std::cout << "G-1: " << std::endl;
        //std::cout << (G-1).print();
        //std::cout << "1-G: " << std::endl;
        G -= 1;
        G -= 2*G;
        //std::cout << (G).print();
        // test Matrix(3,4)
        Matrix F(3,4);
        //std::cout << "F: " << std::endl;
        //std::cout << F.print();
        //std::cout << "dimension: " << F.num_row << " " << F.num_col <<
std::endl;
        F(2,4) = 24;
        //std::cout << "F: " << std::endl;
        //std::cout << (2+F).print();
        //std::cout << "dimension: " << F.num_row << " " << F.num_col <<
std::endl;
        // test Matrix += 4
//	2+F;
        //std::cout << "F: " << std::endl;
//	//std::cout << (2.0+F).print();
//	//std::cout << "dimension: " << F.num_row << " " << F.num_col <<
std::endl;
        F *= G;
//	F += (2*G);
        //std::cout << "F: " << std::endl;
        //std::cout << F.print();
        //std::cout << "dimension: " << F.num_row << " " << F.num_col <<
std::endl;
        // test /
        F = F/2;
        //std::cout << "F: " << std::endl;
        //std::cout << F.print();
        //std::cout << "dimension: " << F.num_row << " " << F.num_col <<
std::endl;

        // test perator =
        Matrix C = A;
        //std::cout << "A:\n " << A.print();
        //std::cout << "C:\n " << C.print();


}

*/
