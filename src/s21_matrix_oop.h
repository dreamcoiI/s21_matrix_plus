

#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>

class S21Matrix {
  public:
  //конструкторы и деструктор
  S21Matrix();
  S21Matrix(const int rows, const int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other) noexcept;
  ~S21Matrix();
  int Get_Rows() const;
  int Get_Cols() const;
  void Set_Rows(const int _rows_);
  void Set_Cols(const int _cols_);
  void Print_Matrix();
  void Init_Matrix();
  bool EqMatrix(const S21Matrix& other);
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose() const;
  S21Matrix Minor_Create(int _rows_, int _cols_) const;
  double Determinant() const;
  double DoubleMatrixDeter() const ;
  S21Matrix CalcComplements() const;
  S21Matrix InverseMatrix() const;
  //операторы перегрузки
  S21Matrix& operator+(const S21Matrix& other);
  S21Matrix& operator-(const S21Matrix& other);
  S21Matrix& operator*(const S21Matrix& other);
  S21Matrix operator*(double num) const;
  bool operator==(const S21Matrix& other);
  S21Matrix& operator=(const S21Matrix& other);
  S21Matrix& operator=(S21Matrix& other) noexcept;
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix& operator*=(const S21Matrix& other);
  S21Matrix& operator*=(const double& num);
  double& operator()(const int i, const int j);
  double operator()(const int i, const int j) const;

 private:
  int rows = 0, cols = 0;
  double** matrix = nullptr;

  void MemoryAllocate();
  void MemoryFree();
};
