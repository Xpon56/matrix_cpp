#ifndef CPP1_S21_MATRIXPLUS_SRC_S21_MATRIX_OOP_H_
#define CPP1_S21_MATRIXPLUS_SRC_S21_MATRIX_OOP_H_

#include <cmath>
#include <cstring>
#include <iostream>

namespace s_21 {
class S21Matrix {
 public:
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other) noexcept;

  S21Matrix& operator=(const S21Matrix& other);
  S21Matrix& operator=(S21Matrix&& other) noexcept;
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix& operator*=(const S21Matrix& other);
  S21Matrix& operator*=(const double num);

  ~S21Matrix();

  int GetRows() const;
  int GetCols() const;
  void SetRows(int rows);
  void SetCols(int cols);

  S21Matrix operator+(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator*(const S21Matrix& other);
  S21Matrix operator*(double num) const;
  friend S21Matrix operator*(double, S21Matrix& matrix);
  bool operator==(const S21Matrix& other);
  double& operator()(int row, int col);
  double& operator()(int row, int col) const;

  bool EqMatrix(const S21Matrix& other);
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);

  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();

 private:
  int rows_, cols_;
  double** matrix_;

  void AllocateMemory();
  void FreeMemory();
  void CopyValues(const S21Matrix& other);
  S21Matrix Minor(int ex_row, int ex_col);
  bool IsMatrixSameDimension(S21Matrix matrix);
  bool IsMatrixSquare();
};
}  // namespace s_21

#endif