

#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdio>

class S21Matrix{
private:
  int row = 0, cols = 0;
  double** matrix = nullptr;

  void MemoryAllocate();
  void MemoryFree();

};
