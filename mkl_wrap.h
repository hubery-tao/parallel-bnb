#pragma once

#include <stdint.h>
#include <mkl_cblas.h>
#include <mkl_lapacke.h>
#include <assert.h>

typedef uint32_t idx_t;

#define IDX_MAX 1024
extern int* int_trash;

void mkl_wrap_init();

struct Vector{
    idx_t n;
    double* vals;
};

struct Matrix {
    idx_t n_rows;
    idx_t n_cols;
    double* vals;
};

void vec_ini(struct Vector* A, idx_t n);
void mat_ini(struct Matrix* A, idx_t m, idx_t n);

static inline double vec_get_pos(const struct Vector* v, idx_t i) {
    return v->vals[i];
}

static inline void vec_set_pos(struct Vector* v, idx_t i, double x) {
    v->vals[i] = x;
}

static inline double mat_get_pos(const struct Matrix* A, idx_t i, idx_t j) {
    return A->vals[A->n_rows*i+j];
}

static inline void mat_set_pos(struct Matrix* A, idx_t i, idx_t j, double x) {
    A->vals[A->n_rows*i+j] = x;
}

void vec_print(const struct Vector* v);
void mat_print(const struct Matrix* A);

// C = A*B
static inline void mat_mul_mat(const struct Matrix* A, const struct Matrix* B, struct Matrix* C) {

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->n_rows, B->n_cols, A->n_cols, 1.0, A->vals, A->n_rows, B->vals, B->n_rows, 0.0, C->vals, C->n_rows);
}

// y = Ax
static inline void mat_mul_vec(const struct Matrix* A, const struct Vector* x, struct Vector* y) {
    cblas_dgemv(CblasRowMajor, CblasNoTrans, A->n_rows, A->n_cols, 1.0, A->vals, A->n_rows, x->vals, 1, 0.0, y->vals, 1);
}

// solve Ax = b, A:n*n, b:n*1, b is modified to x
static inline int lineq_solve(const struct Matrix* A, struct Vector* b) {
    return LAPACKE_dgesv(LAPACK_ROW_MAJOR, A->n_rows, 1, A->vals, A->n_rows, int_trash, b->vals, 1);
}