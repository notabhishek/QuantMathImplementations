#include <iostream>
#include "gaussj.h"

int main()
{
    // 1) Demonstration: Solve A X = B for small 3x3 example
    {
        Matrix A = {
            {2, 1, 0},
            {1, 2, 1},
            {0, 1, 2}
        };
        // Suppose we have 2 right-hand sides in B:
        // e.g. columns: b1 = (5,7,12), b2 = (1,2,3)
        Matrix B = {
            {5, 1},
            {7, 2},
            {12, 3}
        };

        std::cout << "=== Solve System A X = B ===\n";
        printMatrix(A, "Original A");
        printMatrix(B, "Original B (RHS)");

        bool success = gaussj(A, B);
        if (!success) {
            std::cout << "Matrix A was singular!\n\n";
        } else {
            printMatrix(A, "Transformed A (should be A^-1)");
            printMatrix(B, "Solution X (for each column of original B)");
        }
    }

    // 2) Demonstration: Invert a 3x3 matrix directly
    {
        Matrix Ainv = {
            {2, 1, 0},
            {1, 2, 1},
            {0, 1, 2}
        };
        std::cout << "=== Invert A directly ===\n";
        printMatrix(Ainv, "Original Matrix (to invert)");

        bool success = gaussj_invert(Ainv);
        if (!success) {
            std::cout << "Matrix was singular, could not invert!\n\n";
        } else {
            printMatrix(Ainv, "Inverted A");
       }

        std::cout << "=== Invert this Matrix again ===\n";
        bool success2 = gaussj_invert(Ainv);
        if (!success2) {
            std::cout << "Matrix was singular, could not invert!\n\n";
        } else {
            printMatrix(Ainv, "Inverted A");
        }

    }

    return 0;
}
