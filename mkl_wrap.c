#include "mkl_wrap.h"
#include <stdio.h>
#include <stdlib.h>

int* int_trash;

void mkl_wrap_init() {
    int_trash = (int*)malloc(IDX_MAX*sizeof(int));
}

void vec_ini(struct Vector* A, idx_t n) {
    A->n = n;
    A->vals = (double*)malloc(n*sizeof(double));
}

void mat_ini(struct Matrix* A, idx_t m, idx_t n) {
    A->n_rows = m;
    A->n_cols = n;
    A->vals = (double*)malloc(m*n*sizeof(double));
}

void vec_print(const struct Vector* v) {
    for (idx_t i = 0; i < v->n; ++i) {
        printf("%lf    ", vec_get_pos(v, i));
    }
    printf("\n");
}

void mat_print(const struct Matrix* A) {
    for (idx_t i = 0; i < A->n_rows; ++i) {
        for (idx_t j = 0; j < A->n_cols; ++j) {
            printf("%lf    ", mat_get_pos(A, i, j));
        }
        printf("\n");
    }
}