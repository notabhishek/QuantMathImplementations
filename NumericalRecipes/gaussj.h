#pragma once 

#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include <string>

using Matrix = std::vector<std::vector<double>>;


/*
 * Utility to print any n x m matrix 
 */
 void printMatrix(const Matrix &M, const std::string &label = "Matrix")
 {
    std::cout << label << ":\n";
    for (const auto &row : M) {
        for (auto val : row) {
            std::cout << val << "\t";
        }
        std::cout << "\n";
    }
    std::cout << std::endl;
}

/*
 * Perform Gauss-Jordan elimination on a matrix 'a' (n x n)
 * with corresponding right-hand side 'b' (n x m).
 * Based on Numerical Recipes, 3rd Edition.
 *
 * On return:
 *  - 'a' is (ideally) transformed into its inverse if it was nonsingular
 *  - 'b' is transformed into the solution matrix for the original system
 *    A_original * X = B_original.
 * 
 * Returns 'true' if the system is nonsingular, 'false' otherwise.
 */
bool gaussj(Matrix& a, Matrix& b) {
    const int n = static_cast<int>(a.size());
    
    // Verify that the dimensions of a and b are compatible
    assert(n == static_cast<int>(a[0].size()));
    assert(n == static_cast<int>(b.size()));
    
    const int m = static_cast<int>(b[0].size());


    std::vector<bool> pivotUsed(n, false); // Tracks which columns have been used as pivot columns
    std::vector<int> pivotRowIdx(n); // Index i tracks the original row number of i'th pivot element
    std::vector<int> pivotColIdx(n); // Index i tracks the original col number of i'th pivot element

    for(int pivotStep = 0; pivotStep < n; ++pivotStep) {
        // Find the pivot element by searching for the largest element 
        // in an unused pivot column
        
        double maxValue = 0.0;
        int pivotRow = -1;
        int pivotCol = -1;

        for(int c = 0; c < n; ++c) {
            if(!pivotUsed[c]) {
                for(int r = 0; r < n; ++r) {
                    if(!pivotUsed[r] && std::abs(a[r][c]) > maxValue) {
                        maxValue = std::abs(a[r][c]);
                        pivotRow = r;
                        pivotCol = c;
                    }
                }
            }
        }

        if(pivotRow == -1 || pivotCol == -1) {
            // The matrix is singular
            return false;
        }
        
        // Track the pivotStep'th pivot element 
        pivotUsed[pivotCol] = true;
        pivotRowIdx[pivotStep] = pivotRow;
        pivotColIdx[pivotStep] = pivotCol;


        // Move the pivot to a diagonal
        if(pivotRow != pivotCol) {
            // Swap rows pivotRow and pivotCol to take pivot to (pivotCol, pivotCol)
            for(int c = 0; c < n; ++c) {
                std::swap(a[pivotRow][c], a[pivotCol][c]);
            }
            for(int c = 0; c < m; ++c) {
                std::swap(b[pivotRow][c], b[pivotCol][c]);
            }
        }

        const double pivotValue = a[pivotCol][pivotCol];
        assert(pivotValue != 0.0);
        const double pivotInv = 1.0 / pivotValue;

        a[pivotCol][pivotCol] = 1.0; // We know this is 1.0
        // Normalize the pivot row to make a[pivotCol][pivotCol] = 1
        for(int c = 0; c < n; ++c) {
            a[pivotCol][c] *= pivotInv;
        }
        for(int c = 0; c < m; ++c) {
            b[pivotCol][c] *= pivotInv;
        }

        // Elimination
        // Make the other rows zero in the pivot column
        for(int r = 0; r < n; ++r) {
            if(r != pivotCol) {
                // We apply transformation R_r = R_r + factor * R_pivotCol to make a[r][pivotCol] = 0 
                // a[r][pivotCol] + factor * 1 = 0; 
                // factor = -a[r][pivotCol]
                double factor = a[r][pivotCol];
                a[r][pivotCol] = 0.0; // We know this is 0.0
                for(int c = 0; c < n; ++c) {
                    a[r][c] -= factor * a[pivotCol][c];
                }
                for(int c = 0; c < m; ++c) {
                    b[r][c] -= factor * b[pivotCol][c];
                }
            }
        }
        // printMatrix(a, "Matrix A after step " + std::to_string(pivotStep));
        // printMatrix(b, "Matrix B after step " + std::to_string(pivotStep));
    }

    // Standard guass jordan choses ith pivot in row i, but we chose pivot in pivotColIdx[i] in row pivotRowIdx[i]
    // We need to unscrable the rows to get the actual inverse in a 
    for(int step = n - 1; step >= 0; --step) { // in reverse order
        if(pivotRowIdx[step] != pivotColIdx[step]) {
            // Swap columns pivotRowIdx[step] and pivotColIdx[step]
            for(int r = 0; r < n; ++r) {
                std::swap(a[r][pivotRowIdx[step]], a[r][pivotColIdx[step]]);
            }
        }
    }

    return true;
}

// If we are only interested in the inverse of a
bool gaussj(Matrix& a) {
    Matrix b(a.size()); // dummy Matrix with zero columns
    return gaussj(a, b);
}

bool gaussj_invert(Matrix& a) {
    return gaussj(a);
}