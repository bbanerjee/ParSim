#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/ColumnMatrix.h>
#include <gtest/gtest.h>

using namespace Uintah;

TEST(DenseMatrixTest, Construction) {
  DenseMatrix mat(3, 3);
  EXPECT_EQ(mat.nrows(), 3);
  EXPECT_EQ(mat.ncols(), 3);
  
  mat.zero();
  for(int i=0; i<3; ++i) {
    for(int j=0; j<3; ++j) {
      EXPECT_EQ(mat.get(i, j), 0.0);
    }
  }
}

TEST(DenseMatrixTest, PutGet) {
  DenseMatrix mat(2, 2);
  mat.put(0, 0, 1.0);
  mat.put(0, 1, 2.0);
  mat.put(1, 0, 3.0);
  mat.put(1, 1, 4.0);
  
  EXPECT_EQ(mat.get(0, 0), 1.0);
  EXPECT_EQ(mat.get(1, 1), 4.0);
}

TEST(DenseMatrixTest, Identity) {
  DenseMatrix* id = DenseMatrix::identity(3);
  for(int i=0; i<3; ++i) {
    for(int j=0; j<3; ++j) {
      if (i == j) EXPECT_EQ(id->get(i, j), 1.0);
      else EXPECT_EQ(id->get(i, j), 0.0);
    }
  }
  delete id;
}

TEST(ColumnMatrixTest, Construction) {
  ColumnMatrix col(5);
  EXPECT_EQ(col.nrows(), 5);
  EXPECT_EQ(col.ncols(), 1);
  
  col.zero();
  for(int i=0; i<5; ++i) {
    EXPECT_EQ(col[i], 0.0);
  }
}

TEST(ColumnMatrixTest, PutGet) {
  ColumnMatrix col(3);
  col.put(0, 10.0);
  col.put(1, 20.0);
  col.put(2, 30.0);
  
  EXPECT_EQ(col[0], 10.0);
  EXPECT_EQ(col.get(2), 30.0);
}

TEST(MatrixTest, MatrixVectorMult) {
  DenseMatrix mat(2, 2);
  mat.put(0, 0, 1.0); mat.put(0, 1, 2.0);
  mat.put(1, 0, 3.0); mat.put(1, 1, 4.0);
  
  ColumnMatrix x(2);
  x.put(0, 1.0);
  x.put(1, 1.0);
  
  ColumnMatrix b(2);
  int flops, memrefs;
  mat.mult(x, b, flops, memrefs);
  
  // [1 2] [1] = [3]
  // [3 4] [1] = [7]
  EXPECT_EQ(b[0], 3.0);
  EXPECT_EQ(b[1], 7.0);
}
