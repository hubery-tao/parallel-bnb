#pragma once

#include <iostream>

typedef uint32_t idx_t;

class Vector;

class Matrix {

protected:

    idx_t size;
    idx_t num_rows;
    idx_t num_cols;
    double* vals;

    friend class LineqSolve;

public:

    Matrix(idx_t m, idx_t n);
    Matrix(const Matrix& B);
    Matrix(Matrix&& B);

    ~Matrix();

    Matrix& operator=(const Matrix& B);
    Matrix& operator=(Matrix&& B);

    Matrix operator*(const Matrix& B) const;
    Vector operator*(const Vector& vec) const;

    void print() const;

    double get_pos(idx_t i, idx_t j) const {
        return vals[i*num_cols+j];
    }

    void set_pos(idx_t i, idx_t j, double x) {
        vals[i*num_cols+j] = x;
    }

    idx_t get_size() const {
        return size;
    }

    idx_t get_num_rows() const {
        return num_rows;
    }

    idx_t get_num_cols() const {
        return num_cols;
    }

};


class Vector: public Matrix {

    friend class Matrix;
    friend class LineqSolve;

private:
    using Matrix::get_num_rows;
    using Matrix:: get_num_cols;

public:

    Vector(idx_t n): Matrix(n, 1) {}
    Vector(const Vector& vec): Matrix(vec) {}
    Vector(Vector&& vec): Matrix(vec) {}

    double get_pos(idx_t i) const {
        return vals[i];
    }

    void set_pos(idx_t i, double x) const {
        vals[i] = x;
    }

};

class LineqSolve {

    static constexpr idx_t trash_size = 1024;
    static int trash[trash_size];
    
public:

    LineqSolve() = delete;

    // solve Ax = b, and replace b with x, return positive if singular
    static int solve(const Matrix& A, const Vector& b);

};