#include "s21_matrix_oop.h"

//выделение
void S21Matrix ::MemoryAllocate() {
    matrix = new double*[rows]();
    matrix[0] = new double[rows * cols]();
    for (auto i = 1; i != rows; i++){
        matrix[i] = matrix[i-1] +cols;
    }
}

//очистка
void S21Matrix ::MemoryFree() {
    if(matrix && matrix[0]) {
        delete matrix[0];
        delete matrix;
        rows = 0;
        cols = 0;
    }
}

//Базовый конструктор, инициализирующий матрицу некоторой заранее заданной размерностью
S21Matrix ::S21Matrix() {
    rows = 3;
    cols = 3;
    MemoryAllocate();
}

//Параметризированный конструктор с количеством строк и столбцов
S21Matrix ::S21Matrix(const int rows, const int cols)
    : rows(rows), cols(cols)
 {
    if(rows< 1 || cols < 1) {
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
    other.matrix = nullptr;
    other.cols = 0;
    other.rows = 0;
}

// деструктор
S21Matrix ::~S21Matrix() {
    MemoryFree();
}

//Присвоение матрице значений другой матрицы путем копирования
S21Matrix& S21Matrix ::operator=(const S21Matrix&other) {
    if (this == &other) {
        ;
    } else {
        this->MemoryFree();
        cols = other.cols;
        rows = other.rows;
        for (auto i = 0; i < rows; i++) {
            for (auto j = 0; j < cols; j++) {
                matrix[i][j] = other.matrix[i][j];
            }
        }
    }
    return *this;
}

S21Matrix& S21Matrix ::operator=(S21Matrix&other) noexcept {
    if (this == &other) {
        ;
    } else {
        this->MemoryFree();
        cols = other.cols;
        rows = other.rows;
        matrix = other.matrix;
        other.matrix = nullptr;
    }
    return *this;
}

//Проверяет матрицы на равенство между собой 
bool S21Matrix::EqMatrix(const S21Matrix& other) {
    if (rows != other.rows || cols != other.cols) return false;
    for(auto i = 0; i < rows;i++) {
        for (auto j = 0; j < cols; j++) {
            if(fabs(matrix[i][j] - other.matrix[i][j]) > 01e-7) return false;
        }
    }
    return true;
}

//Прибавляет вторую матрицы к текущей
void S21Matrix::SumMatrix(const S21Matrix& other) {
    if (rows != other.rows || cols != other.cols) throw std::logic_error("Error, incorrect values");
    for(auto i = 0; i < rows;i++) {
        for (auto j = 0; j < cols; j++) {
            matrix[i][j] +=other.matrix[i][j];
        }
    }
}

//Вычитает из текущей матрицы другую
void S21Matrix::SubMatrix(const S21Matrix& other) {
    if (rows != other.rows || cols != other.cols) throw std::logic_error("Error, incorrect values");
    for(auto i = 0; i < rows;i++) {
        for (auto j = 0; j < cols; j++) {
            matrix[i][j] -=other.matrix[i][j];
        }
    }
}

//Умножает текущую матрицу на число
void S21Matrix::MulNumber(const double num) {
    for(auto i = 0; i < rows;i++) {
        for (auto j = 0; j < cols; j++) {
            matrix[i][j] *=num;
        }
    }
}

//Умножает текущую матрицу на вторую
void S21Matrix::MulMatrix(const S21Matrix& other) {
    if (rows != other.cols || cols != other.rows) throw std::logic_error("Error, incorrect values");
    S21Matrix finally(rows,other.cols);
    for(auto i = 0; i < rows;i++) {
        for (auto j = 0; j < cols; j++) {
            for(auto k = 0; k < other.rows;k++){
                finally.matrix[i][j] += matrix[i][k] * other.matrix[k][j];   
            }
        }
    }        
}

