#include "mkl_wrap.hpp"

#include <mkl_cblas.h>
#include <mkl_lapacke.h>
#include <cstring>

Matrix::Matrix(idx_t m, idx_t n): size(m*n), num_rows(m), num_cols(n) {
    vals = new double[size];
}

Matrix::Matrix(const Matrix& B): size(B.size), num_rows(B.num_rows), num_cols(B.num_cols) {
    vals = new double[size];
    memcpy(vals, B.vals, size*sizeof(double));
}

Matrix::Matrix(Matrix&& B): size(B.size), num_rows(B.num_rows), num_cols(B.num_cols), vals(B.vals) {
    memset(&B, 0, sizeof(B));
}

Matrix::~Matrix() {
    delete[] vals;
    memset(this, 0, sizeof(*this));
}

Matrix& Matrix::operator=(const Matrix& B) {
    if (this == &B) {
        return *this;
    }
    memcpy(vals, B.vals, size*sizeof(double));
    return *this;
}

Matrix& Matrix::operator=(Matrix&& B) {
    if (this == &B) {
        return *this;
    }
    delete[] vals;
    vals = B.vals;
    memset(&B, 0, sizeof(B));
    return *this;
}

Matrix Matrix::operator*(const Matrix& B) const {
    Matrix C(num_rows, B.num_cols);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, num_rows, B.num_cols, num_cols, 
                1.0, vals, num_cols, B.vals, B.num_cols, 0.0, C.vals, C.num_cols);
    return C;
}

Vector Matrix::operator*(const Vector& vec) const {
    Vector res(num_rows);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, num_rows, num_cols, 
                1.0, vals, num_cols, vec.vals, 1, 0.0, res.vals, 1);
    return res;
}

void Matrix::print() const {
    for (idx_t n = 0, j = 0; n < size; ++n) {
        std::cout << vals[n] << "    ";
        if (++j == num_cols) {
            std::cout << std::endl;
            j = 0;
        }
    }
}

int LineqSolve::trash[LineqSolve::trash_size];


int LineqSolve::solve(const Matrix& A, const Vector& b){
    return LAPACKE_dgesv(LAPACK_ROW_MAJOR, A.num_rows, 1, A.vals, A.num_cols, trash, b.vals, 1);
}