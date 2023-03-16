#include "s21_matrix_oop.h"

S21Matrix ::S21Matrix() {
    row = 3;
    cols = 3;
    MemoryAllocate();
}


void S21Matrix ::MemoryAllocate() {
    matrix = new double*[row]();
    matrix[0] = new double[row * cols]();
    for (auto i = 1; i != row; i++){
        matrix[i] = matrix[i-1] +cols;
    }
}

void S21Matrix ::MemoryFree() {
    if(matrix && matrix[0]) {
        delete matrix[0];
        delete matrix;
        row = 0;
        cols = 0;
    }
}