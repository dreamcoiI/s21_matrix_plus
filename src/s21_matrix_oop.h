

#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdio>

class S21Matrix{
//конструкторы и деструктор
S21Matrix();
S21Matrix(const int rows, const int cols);
S21Matrix(const S21Matrix& other);
S21Matrix(S21Matrix&& other);
~S21Matrix();
bool EqMatrix(const S21Matrix& other);
void SumMatrix(const S21Matrix& other);
void SubMatrix(const S21Matrix& other);
void MulNumber(const double num);
void MulMatrix(const S21Matrix& other);
//операторы перегрузки 
S21Matrix& S21Matrix ::operator=(const S21Matrix&other);
S21Matrix& S21Matrix ::operator=(S21Matrix&other) noexcept;
private:
  int rows = 0, cols = 0;
  double** matrix = nullptr;

  void MemoryAllocate();
  void MemoryFree();

};
