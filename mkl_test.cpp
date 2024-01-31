#include "mkl_wrap.hpp"
#include <stdio.h>

int main() {

    Matrix A(2,2);
    for (idx_t i = 0; i < A.get_num_rows(); ++i) {
        for (idx_t j = 0; j < A.get_num_cols(); ++j) {
            A.set_pos(i, j, i*A.get_num_rows()+j+1);
        }
    }
    A.print();

    (A*A).print();

    Vector vec(2);
    for (idx_t i = 0; i < vec.get_size(); ++i) {
        vec.set_pos(i, i+1);
    }

    Vector b = A*vec;

    LineqSolve::solve(A, b);

    b.print();

    return 0;
}
