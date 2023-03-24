#include <gtest/gtest.h>

#include "s21_matrix_oop.h"

int main() {
  testing::InitGoogleTest();
  return RUN_ALL_TESTS();
}

namespace {

TEST(Constructor, Default) {
  auto test = S21Matrix();
  EXPECT_EQ(test.Get_Rows(), 3);
  EXPECT_EQ(test.Get_Cols(), 3);
}

TEST(Constructor, Args) {
  auto test = S21Matrix(2, 4);
  EXPECT_EQ(test.Get_Rows(), 2);
  EXPECT_EQ(test.Get_Cols(), 4);
}

TEST(Constructor, Copy) {
  auto test = S21Matrix(2, 4);
  test(0, 0) = 3;
  auto test2 = (test);
  EXPECT_EQ(test2.Get_Rows(), 2);
  EXPECT_EQ(test2.Get_Cols(), 4);
  EXPECT_EQ(test2(0, 0), 3);
}

TEST(Constructor, Move) {
  auto test = S21Matrix(2, 4);
  test(0, 0) = 3;
  auto test2 = std::move(test);
  EXPECT_EQ(test2.Get_Rows(), 2);
  EXPECT_EQ(test2.Get_Cols(), 4);
  EXPECT_EQ(test2(0, 0), 3);
}

TEST(Setter, Rows) {
  auto test = S21Matrix(2, 2);
  test(0, 0) = 3;
  test(0, 1) = 1;
  test(1, 0) = 2;
  test(1, 1) = -2;
  EXPECT_EQ(test.Get_Rows(), 2);
  test.Set_Rows(3);
  EXPECT_EQ(test.Get_Rows(), 3);
  EXPECT_EQ(test(1, 1), -2);
  test.Set_Rows(1);
  test.Set_Cols(1);
  test.Set_Rows(2);
  test.Set_Cols(2);
  EXPECT_EQ(test(1, 1), 0);
}

TEST(Setter, Cols) {
  auto test = S21Matrix(2, 1);
  test(0, 0) = 3;
  test(1, 0) = 2;
  EXPECT_EQ(test.Get_Cols(), 1);
  test.Set_Cols(10);
  EXPECT_EQ(test.Get_Cols(), 10);
  EXPECT_EQ(test(1, 9), 0);
  EXPECT_EQ(test(0, 0), 3);
}

TEST(Arithmethic, EqMatrix) {
  auto test1 = S21Matrix();
  auto test2 = S21Matrix();
  test1.Init_Matrix();
  test2.Init_Matrix();
  EXPECT_EQ(test1.EqMatrix(test2), true);
  EXPECT_EQ(test2.EqMatrix(test1), true);
  EXPECT_EQ(test1 == test2, true);
  EXPECT_EQ(test2 == test1, true);
  test1(0, 0) = 10;
  EXPECT_EQ(test1.EqMatrix(test2), false);
  EXPECT_EQ(test2.EqMatrix(test1), false);
  EXPECT_EQ(test1 == test2, false);
  EXPECT_EQ(test2 == test1, false);
}

TEST(Arithmethic, SumMatrix) {
  auto test1 = S21Matrix(2, 2);
  auto test2 = S21Matrix(2, 2);
  test1.Init_Matrix();
  test2.Init_Matrix();
  test1 += test2;
  EXPECT_EQ(test1(0, 0), 0);
  EXPECT_EQ(test1(1, 0), 2);
  EXPECT_EQ(test1(0, 1), 2);
  EXPECT_EQ(test1(1, 1), 4);
}

TEST(Arithmethic, SubMatrix) {
  auto test1 = S21Matrix(2, 2);
  auto test2 = S21Matrix(2, 2);
  test1.Init_Matrix();
  test2.Init_Matrix();
  test1 *= 3;
  test1 -= test2;
  EXPECT_EQ(test1(0, 0), 0);
  EXPECT_EQ(test1(1, 0), 2);
  EXPECT_EQ(test1(0, 1), 2);
  EXPECT_EQ(test1(1, 1), 4);
}

TEST(Arithmethic, MulMatrix) {
  auto test1 = S21Matrix(2, 2);
  auto test2 = S21Matrix(2, 2);
  test1.Init_Matrix();
  test2.Init_Matrix();
  test1 *= test2;
  EXPECT_EQ(test1(0, 0), 1);
  EXPECT_EQ(test1(1, 0), 2);
  EXPECT_EQ(test1(0, 1), 2);
  EXPECT_EQ(test1(1, 1), 5);
}

TEST(Arithmethic, MulNum) {
  auto test1 = S21Matrix(2, 2);
  auto test2 = S21Matrix(2, 2);
  test1.Init_Matrix();
  test2.Init_Matrix();
  auto test3 = 3 * test1;
  auto test4 = test2 * 3;
  EXPECT_EQ(test3(0, 0), 0);
  EXPECT_EQ(test3(1, 0), 3);
  EXPECT_EQ(test3(0, 1), 3);
  EXPECT_EQ(test3(1, 1), 6);
  EXPECT_EQ(test4(0, 0), 0);
  EXPECT_EQ(test4(1, 0), 3);
  EXPECT_EQ(test4(0, 1), 3);
  EXPECT_EQ(test4(1, 1), 6);
}

TEST(Actions, Transpose) {
  auto test1 = S21Matrix(2, 2);
  test1(0, 0) = 3;
  test1(0, 1) = 1;
  test1(1, 0) = 20;
  test1(1, 1) = -2;
  test1 = test1.Transpose();
  EXPECT_EQ(test1(0, 0), 3);
  EXPECT_EQ(test1(0, 1), 20);
  EXPECT_EQ(test1(1, 0), 1);
  EXPECT_EQ(test1(1, 1), -2);
}

TEST(Actions, Complements) {
  auto test1 = S21Matrix(2, 2);
  test1(0, 0) = 3;
  test1(0, 1) = 1;
  test1(1, 0) = 20;
  test1(1, 1) = -2;
  test1 = test1.CalcComplements();
  EXPECT_EQ(test1(0, 0), -2);
  EXPECT_EQ(test1(0, 1), -20);
  EXPECT_EQ(test1(1, 0), -1);
  EXPECT_EQ(test1(1, 1), 3);
}

TEST(Actions, Determinant) {
  auto test1 = S21Matrix(2, 2);
  test1(0, 0) = 3;
  test1(0, 1) = 1;
  test1(1, 0) = 20;
  test1(1, 1) = -2;
  EXPECT_EQ(test1.Determinant(), -26);
  test1(0, 0) = 2;
  test1(0, 1) = 1;
  test1(1, 0) = 4;
  test1(1, 1) = -2;
  EXPECT_EQ(test1.Determinant(), -8);
  test1(0, 0) = 2;
  test1(0, 1) = -5;
  test1(1, 0) = 1;
  test1(1, 1) = -2;
  EXPECT_EQ(test1.Determinant(), 1);
  test1(0, 0) = 3;
  test1(0, 1) = -1;
  test1(1, 0) = 7;
  test1(1, 1) = -2;
  EXPECT_EQ(test1.Determinant(), 1);
}

TEST(Actions, Inverse) {
  auto test1 = S21Matrix(2, 2);
  test1(0, 0) = 3;
  test1(0, 1) = 1;
  test1(1, 0) = 2;
  test1(1, 1) = -2;
  test1 = test1.InverseMatrix();
  EXPECT_EQ(test1(0, 0), 0.25);
  EXPECT_EQ(test1(0, 1), 0.25);
  EXPECT_EQ(test1(1, 0), 0.125);
  EXPECT_EQ(test1(1, 1), -0.375);

  test1 = S21Matrix(3, 3);
  test1(0, 0) = 2;
  test1(0, 1) = 5;
  test1(0, 2) = 7;
  test1(1, 0) = 6;
  test1(1, 1) = 3;
  test1(1, 2) = 4;
  test1(2, 0) = 5;
  test1(2, 1) = -2;
  test1(2, 2) = -3;

  test1 = test1.InverseMatrix();
  EXPECT_EQ(test1(0, 0), 1);
  EXPECT_EQ(test1(0, 1), -38);
  EXPECT_EQ(test1(1, 0), -1);
  EXPECT_EQ(test1(1, 1), 41);
}

}  // namespace