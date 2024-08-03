#include "s21_matrix_oop.h"

namespace s_21 {

S21Matrix::S21Matrix() : S21Matrix(5, 5){};

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows <= 0 || cols <= 0) {
    throw std::invalid_argument(
        "CreationError: The number of rows or cols cannot be less than 1");
  }
  AllocateMemory();
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : rows_(other.rows_), cols_(other.cols_) {
  AllocateMemory();
  CopyValues(other);
}

S21Matrix::S21Matrix(S21Matrix&& other) noexcept
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (this != &other) {
    FreeMemory();
    rows_ = other.rows_;
    cols_ = other.cols_;
    AllocateMemory();
    CopyValues(other);
  }
  return *this;
}

S21Matrix& S21Matrix::operator=(S21Matrix&& other) noexcept {
  if (this != &other) {
    FreeMemory();

    rows_ = other.rows_;
    cols_ = other.cols_;
    matrix_ = other.matrix_;

    other.rows_ = 0;
    other.cols_ = 0;
    other.matrix_ = nullptr;
  }

  return *this;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  SumMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  SubMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  MulMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const double num) {
  MulNumber(num);
  return *this;
}

S21Matrix::~S21Matrix() { FreeMemory(); }

int S21Matrix::GetRows() const { return rows_; }

int S21Matrix::GetCols() const { return cols_; }

void S21Matrix::SetRows(int rows) {
  if (rows <= 0) {
    throw std::invalid_argument(
        "SettingRowsError: The number of rows cannot be less than 1");
  }

  S21Matrix temp(rows, cols_);
  int rows_range = rows < rows_ ? rows : rows_;
  for (int i = 0; i < rows_range; i++) {
    for (int j = 0; j < cols_; j++) {
      temp.matrix_[i][j] = matrix_[i][j];
    }
  }

  *this = temp;
}

void S21Matrix::SetCols(int cols) {
  if (cols <= 0) {
    throw std::invalid_argument(
        "SettingColsError: The number of cols cannot be less than 1");
  }

  S21Matrix tmp(rows_, cols);
  int cols_range = cols < cols_ ? cols : cols_;
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_range; j++) {
      tmp.matrix_[i][j] = matrix_[i][j];
    }
  }

  *this = tmp;
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  S21Matrix res_matrix(*this);
  res_matrix.SumMatrix(other);
  return res_matrix;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  S21Matrix res_matrix(*this);
  res_matrix.SubMatrix(other);
  return res_matrix;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  S21Matrix res_matrix(*this);
  res_matrix.MulMatrix(other);
  return res_matrix;
}

S21Matrix S21Matrix::operator*(double num) const {
  S21Matrix res_matrix(*this);
  res_matrix.MulNumber(num);
  return res_matrix;
}

S21Matrix operator*(double num, S21Matrix& matrix) { return matrix * num; }

bool S21Matrix::operator==(const S21Matrix& other) { return EqMatrix(other); }

double& S21Matrix::operator()(int row, int col) {
  if (row < 0 || col < 0 || row >= rows_ || col >= cols_) {
    throw std::out_of_range("InvalidIndexError: Index is out of range");
  }

  return matrix_[row][col];
}

double& S21Matrix::operator()(int row, int col) const {
  if (row < 0 || col < 0 || row >= rows_ || col >= cols_) {
    throw std::out_of_range("InvalidIndexError: Index is out of range");
  }

  return matrix_[row][col];
}

bool S21Matrix::EqMatrix(const S21Matrix& other) {
  return std::memcmp(matrix_ + rows_, other.matrix_ + rows_, rows_ * cols_) ==
             0 &&
         IsMatrixSameDimension(other);
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (!IsMatrixSameDimension(other)) {
    throw std::range_error("SumMatrixError: Matrices of different dimensions");
  }
  for (int rowIndex = 0; rowIndex < rows_; rowIndex++) {
    for (int colIndex = 0; colIndex < cols_; colIndex++) {
      matrix_[rowIndex][colIndex] += other.matrix_[rowIndex][colIndex];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (!IsMatrixSameDimension(other)) {
    throw std::range_error("SubMatrixError: Matrices of different dimensions");
  }
  for (int rowIndex = 0; rowIndex < rows_; rowIndex++) {
    for (int colIndex = 0; colIndex < cols_; colIndex++) {
      matrix_[rowIndex][colIndex] -= other.matrix_[rowIndex][colIndex];
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int rowIndex = 0; rowIndex < rows_; rowIndex++) {
    for (int colIndex = 0; colIndex < cols_; colIndex++) {
      matrix_[rowIndex][colIndex] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ != other.rows_) {
    throw std::range_error(
        "MulMatrixError: Incorrect dimensions to multiply two matrices");
  }

  S21Matrix resultMatrix(rows_, other.cols_);

  for (int rowIndex = 0; rowIndex < rows_; rowIndex++) {
    for (int colIndex = 0; colIndex < other.cols_; colIndex++) {
      double sum = 0.0;
      for (int k = 0; k < cols_; k++) {
        sum += matrix_[rowIndex][k] * other.matrix_[k][colIndex];
      }
      resultMatrix.matrix_[rowIndex][colIndex] = sum;
    }
  }

  *this = resultMatrix;
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix transposedMatrix(cols_, rows_);
  for (int rowIndex = 0; rowIndex < rows_; rowIndex++) {
    for (int colIndex = 0; colIndex < cols_; colIndex++) {
      transposedMatrix.matrix_[colIndex][rowIndex] =
          matrix_[rowIndex][colIndex];
    }
  }

  return transposedMatrix;
}

S21Matrix S21Matrix::CalcComplements() {
  if (!IsMatrixSquare()) {
    throw std::range_error("CalcComplementsError: The matrix must be square");
  }

  S21Matrix complementsMatrix(rows_, cols_);

  if (rows_ == 1) {
    complementsMatrix.matrix_[0][0] = matrix_[0][0];
  } else {
    for (int rowIndex = 0; rowIndex < rows_; rowIndex++) {
      int sign = (rowIndex % 2 == 0) ? 1 : -1;

      for (int colIndex = 0; colIndex < cols_; colIndex++) {
        S21Matrix minor(Minor(rowIndex, colIndex));
        complementsMatrix.matrix_[rowIndex][colIndex] =
            sign * minor.Determinant();
        sign = -sign;
      }
    }
  }

  return complementsMatrix;
}

double S21Matrix::Determinant() {
  if (!IsMatrixSquare()) {
    throw std::range_error("DeterminantError: The matrix must be square");
  }

  double determinant = 0.0;

  if (rows_ == 1) {
    determinant = matrix_[0][0];
  } else if (rows_ == 2) {
    determinant = matrix_[0][0] * matrix_[1][1] - matrix_[1][0] * matrix_[0][1];
  } else {
    int sign = 1;
    for (int column = 0; column < cols_; column++) {
      S21Matrix minor(Minor(0, column));
      determinant += sign * matrix_[0][column] * minor.Determinant();
      sign = -sign;
    }
  }

  return determinant;
}

S21Matrix S21Matrix::InverseMatrix() {
  double determinant = Determinant();

  if (determinant == 0 || !IsMatrixSquare()) {
    throw std::range_error(
        "InverseError: Incompatible matrix sizes to compute inverse matrix");
  }

  S21Matrix inverseMatrix(rows_, cols_);

  if (rows_ == 1) {
    inverseMatrix.matrix_[0][0] = 1.0 / matrix_[0][0];
  } else {
    inverseMatrix = CalcComplements().Transpose();

    inverseMatrix.MulNumber(1.0 / determinant);
  }

  return inverseMatrix;
}

void S21Matrix::AllocateMemory() {
  matrix_ = new double*[rows_ + rows_ * cols_]();
  double* matrixDataStart = reinterpret_cast<double*>(matrix_ + rows_);
  int iterationLimit = (rows_ <= cols_) ? cols_ : rows_;
  for (int currentRow = 0; currentRow < iterationLimit; currentRow++) {
    matrix_[currentRow] = matrixDataStart + currentRow * cols_;
  }
  std::memset(matrixDataStart, 0, rows_ * cols_ * sizeof(double));
}

void S21Matrix::FreeMemory() {
  if (matrix_) {
    delete[] matrix_;
    matrix_ = nullptr;
  }
}

void S21Matrix::CopyValues(const S21Matrix& sourceMatrix) {
  for (int currentRow = 0; currentRow < rows_; currentRow++) {
    for (int currentCol = 0; currentCol < cols_; currentCol++) {
      matrix_[currentRow][currentCol] =
          sourceMatrix.matrix_[currentRow][currentCol];
    }
  }
}

S21Matrix S21Matrix::Minor(int excluded_row, int excluded_col) {
  S21Matrix minorMatrix(rows_ - 1, cols_ - 1);

  for (int current_row = 0, minor_row = 0; current_row < rows_; current_row++) {
    if (current_row != excluded_row) {
      for (int current_col = 0, minor_col = 0; current_col < cols_;
           current_col++) {
        if (current_col != excluded_col) {
          minorMatrix.matrix_[minor_row][minor_col] =
              matrix_[current_row][current_col];
          minor_col++;
        }
      }
      minor_row++;
    }
  }

  return minorMatrix;
}

bool S21Matrix::IsMatrixSameDimension(S21Matrix matrix) {
  return (rows_ == matrix.rows_ && cols_ == matrix.cols_);
}

bool S21Matrix::IsMatrixSquare() { return (cols_ == rows_); }

}  // namespace s_21
