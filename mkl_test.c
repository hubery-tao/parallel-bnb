#include "mkl_wrap.h"
#include <stdio.h>

int main() {

    mkl_wrap_init();

    struct Matrix A;
    mat_ini(&A, 1, 1);
    A.vals[0] = 12;

    mat_print(&A);

    struct Vector x;
    vec_ini(&x, 1);
    x.vals[0] = 3;

    vec_print(&x);

    // struct Vector y;
    // vec_ini(&y, 2);

    int info = lineq_solve(&A, &x);

    printf("%d\n", info);

    vec_print(&x);

    return 0;
}
