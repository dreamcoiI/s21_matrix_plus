#include "s21_matrix_oop.h"

//аллоцирование
void S21Matrix ::MemoryAllocate() {
  matrix = new double*[rows]();
  matrix[0] = new double[rows * cols]();
  for (auto i = 1; i != rows; i++) {
    matrix[i] = matrix[i - 1] + cols;
  }
}

//очистка
void S21Matrix ::MemoryFree() {
  if (matrix && matrix[0]) {
    delete matrix[0];
    delete matrix;
    rows = 0;
    cols = 0;
  }
}
//геттеры, сеттеры, инициализаци и отпринтовка матриц
int S21Matrix ::Get_Rows() const { return rows; }

int S21Matrix ::Get_Cols() const { return cols; }

void S21Matrix ::Set_Rows(const int _rows_) {
  if (_rows_ < 1) throw std::logic_error("Error, index is out of range");
  if (rows != _rows_) {
    S21Matrix tmp_matrix(_rows_, cols);
    for (int i = 0; i < rows && i < tmp_matrix.rows; i++) {
      for (int j = 0; j < cols && j < tmp_matrix.cols; j++) {
        tmp_matrix.matrix[i][j] = matrix[i][j];
      }
    }
    *this = std::move(tmp_matrix);
  }
}

void S21Matrix ::Set_Cols(const int _cols_) {
  if (_cols_ < 1) throw std::logic_error("Error, index is out of range");
  if (cols != _cols_) {
    S21Matrix tmp_matrix(rows, _cols_);
    for (int i = 0; i < rows && i < tmp_matrix.rows; i++) {
      for (int j = 0; j < cols && j < tmp_matrix.cols; j++) {
        tmp_matrix.matrix[i][j] = matrix[i][j];
      }
    }
    *this = std::move(tmp_matrix);
  }
}

void S21Matrix::Print_Matrix() {
  for (auto i = 0; i < rows; i++) {
    for (auto j = 0; j < cols; j++) {
      std ::cout << matrix[i][j] << " ";
    }
    std ::cout << std ::endl;
  }
  std ::cout << std ::endl;
}

void S21Matrix::Init_Matrix() {
  for (auto i = 0; i < rows; i++) {
    for (auto j = 0; j < cols; j++) {
      matrix[i][j] = i + j;
    }
  }
}
//Базовый конструктор, инициализирующий матрицу некоторой заранее заданной
//размерностью
S21Matrix ::S21Matrix() {
  rows = 3;
  cols = 3;
  MemoryAllocate();
}

//Параметризированный конструктор с количеством строк и столбцов
S21Matrix ::S21Matrix(const int rows, const int cols) : rows(rows), cols(cols) {
  if (rows < 1 || cols < 1) {
    std::logic_error("ERROR! INCORRECT MEANING");
  } else {
    MemoryAllocate();
  }
}

//конструктор копирования
S21Matrix ::S21Matrix(const S21Matrix& other)
    : rows(other.rows), cols(other.cols) {
  MemoryAllocate();
  for (auto i = 0; i < rows; i++) {
    for (auto j = 0; j < cols; j++) {
      matrix[i][j] = other.matrix[i][j];
    }
  }
}

//конструктор переноса из одной матрицы в другую
S21Matrix ::S21Matrix(S21Matrix&& other) noexcept {
  rows = other.rows;
  cols = other.cols;
  matrix = other.matrix;
  other.cols = 0;
  other.rows = 0;
  other.matrix = nullptr;
}

// деструктор
S21Matrix ::~S21Matrix() { MemoryFree(); }

//Проверяет матрицы на равенство между собой
bool S21Matrix::EqMatrix(const S21Matrix& other) {
  if (rows != other.rows || cols != other.cols) return false;
  for (auto i = 0; i < rows; i++) {
    for (auto j = 0; j < cols; j++) {
      if (fabs(matrix[i][j] - other.matrix[i][j]) > 01e-7) return false;
    }
  }
  return true;
}

//Прибавляет вторую матрицы к текущей
void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (rows != other.rows || cols != other.cols)
    throw std::logic_error("Error, incorrect values");
  for (auto i = 0; i < rows; i++) {
    for (auto j = 0; j < cols; j++) {
      matrix[i][j] += other.matrix[i][j];
    }
  }
}

//Вычитает из текущей матрицы другую
void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (rows != other.rows || cols != other.cols)
    throw std::logic_error("Error, incorrect values");
  for (auto i = 0; i < rows; i++) {
    for (auto j = 0; j < cols; j++) {
      matrix[i][j] -= other.matrix[i][j];
    }
  }
}

//Умножает текущую матрицу на число
void S21Matrix::MulNumber(const double num) {
  for (auto i = 0; i < rows; i++) {
    for (auto j = 0; j < cols; j++) {
      matrix[i][j] *= num;
    }
  }
}

//Умножает текущую матрицу на вторую
void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (rows != other.cols || cols != other.rows)
    throw std::logic_error("Error, incorrect values");
  S21Matrix Mul_matrix(rows, other.cols);
  for (auto i = 0; i < rows; i++) {
    for (auto j = 0; j < cols; j++) {
      for (auto k = 0; k < other.rows; k++) {
        Mul_matrix.matrix[i][j] += matrix[i][k] * other.matrix[k][j];
      }
    }
  }
  *this = std::move(Mul_matrix);
}

//Создает новую транспонированную матрицу из текущей и возвращает ее
S21Matrix S21Matrix ::Transpose() const {
  S21Matrix transponse(rows, cols);
  for (auto i = 0; i < rows; i++) {
    for (auto j = 0; j < cols; j++) {
      transponse.matrix[i][j] = matrix[j][i];
    }
  }
  return transponse;
}
//создаем матрицу миноров
S21Matrix S21Matrix::Minor_Create(int _rows_, int _cols_) const {
  S21Matrix minor(rows - 1, cols - 1);
  _rows_--;
  _cols_--;
  auto rows_min = 0;
  for (auto i = 0; i < rows; i++) {
    auto cols_min = 0;
    for (auto j = 0; j < cols; j++) {
      if (i != _rows_ && j != _cols_) {
        minor.matrix[rows_min][cols_min] = this->matrix[i][j];
        cols_min++;
      }
    }
    if (i != _rows_) rows_min++;
  }
  return minor;
}

//Вычисляет и возвращает определитель текущей матрицы
double S21Matrix::Determinant() const {
  if (rows != cols) throw std::logic_error("Error, rows not equal columns");
  double res = 0;
  if (rows == 1) {
    res = matrix[0][0];
  } else {
    double tmp = 0;
    int sign = 1;
    for (auto i = 0; i < cols; i++) {
      S21Matrix minor = this->Minor_Create(1, i + 1);
      if (minor.rows == 2) {
        tmp = minor.DoubleMatrixDeter();
      } else {
        tmp = minor.Determinant();
      }
      res += matrix[0][i] * tmp * sign;
      sign*=-1;
    }
  }
  return res;
}

double S21Matrix::DoubleMatrixDeter() const {
  return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
}

//Вычисляет матрицу алгебраических дополнений текущей матрицы и возвращает ее
S21Matrix S21Matrix ::CalcComplements() const {
  if (rows != cols) throw std::logic_error("Error, incorrect values");
  S21Matrix Calc_Complements(rows, cols);
  for (auto i = 0; i < rows; i++) {
    for (auto j = 0; j < cols; j++) {
      S21Matrix Minor_matrix = this->Minor_Create(i + 1, j + 1);
      Calc_Complements.matrix[i][j] =
          Minor_matrix.Determinant() * pow(-1, (i + j));
    }
  }
  return Calc_Complements;
}

S21Matrix S21Matrix ::InverseMatrix() const {
  double deter = 0;
  deter = this->Determinant();
  if (fabs(deter- 0)  < 01e-7) throw std::logic_error("Determinant is null");
  S21Matrix Inverse_Matrix = this->CalcComplements();
  Inverse_Matrix.Transpose();
  for (auto i = 0; i < rows; i++) {
    for (auto j = 0; j < rows; j++) {
      Inverse_Matrix.matrix[i][j] = Inverse_Matrix.matrix[i][j] / deter;
    }
  }
  return Inverse_Matrix;
}

//Сложение двух матриц
S21Matrix& S21Matrix ::operator+(const S21Matrix& other) {
  this->SumMatrix(other);
  return *this;
}

//Вычитание одной матрицы из другой
S21Matrix& S21Matrix ::operator-(const S21Matrix& other) {
  this->SubMatrix(other);
  return *this;
}

//Умножение матриц
S21Matrix& S21Matrix ::operator*(const S21Matrix& other) {
  this->MulMatrix(other);
  return *this;
}

//умножение матрицы на число
S21Matrix S21Matrix::operator*(double num) const {
  S21Matrix result = *this;
  result *= num;
  return result;
}

//Проверка на равенство матриц (EqMatrix)
bool S21Matrix ::operator==(const S21Matrix& other) {
  return this->EqMatrix(other);
}

//Присвоение матрице значений другой матрицы путем копирования
S21Matrix& S21Matrix ::operator=(const S21Matrix& other) {
  if (this != &other) {
    this->MemoryFree();
    cols = other.cols;
    rows = other.rows;
    MemoryAllocate();
    for (auto i = 0; i < rows; i++) {
      for (auto j = 0; j < cols; j++) {
        matrix[i][j] = other.matrix[i][j];
      }
    }
  }
  return *this;
}

//Присвоение матрице значений другой матрицы
S21Matrix& S21Matrix ::operator=(S21Matrix& other) noexcept {
  if (this != &other) {
    this->MemoryFree();
    cols = other.cols;
    rows = other.rows;
    matrix = other.matrix;
    other.matrix = nullptr;
  }
  return *this;
}

//Присвоение сложения (SumMatrix)
S21Matrix& S21Matrix ::operator+=(const S21Matrix& other) {
  this->SumMatrix(other);
  return *this;
}

//Присвоение разности (SubMatrix)
S21Matrix& S21Matrix ::operator-=(const S21Matrix& other) {
  this->SubMatrix(other);
  return *this;
}

//Присвоение умножения (MulMatrix)
S21Matrix& S21Matrix ::operator*=(const S21Matrix& other) {
  MulMatrix(other);
  return *this;
}

//Присвоение умножения (MulNumber)
S21Matrix& S21Matrix ::operator*=(const double& num) {
  MulNumber(num);
  return *this;
}

//Индексация по элементам матрицы(строки, столбцы)
double& S21Matrix::operator()(const int i, const int j) {
  if (i > rows || i < 0 || j > cols || j < 0)
    throw std::logic_error("Error, index is out of range");
  return matrix[i][j];
}

double S21Matrix::operator()(const int i, const int j) const {
  if (i > rows || i < 0 || j > cols || j < 0)
    throw std::logic_error("Error, index is out of range");
  return matrix[i][j];
}