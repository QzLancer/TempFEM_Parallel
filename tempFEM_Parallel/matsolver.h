#ifndef MATSOLVER_H
#define MATSOLVER_H
#include "datatype.h"
#if !defined(ARMA_32BIT_WORD)
#define ARMA_32BIT_WORD
#endif
#include <armadillo>

using namespace arma;

class MatSolver
{
public:
    MatSolver();
    virtual double *solveMatrix(umat locs, mat vals, vec F, int size) = 0;
    virtual double *solveMatrix_LU(vec F) = 0;
    virtual ~MatSolver() = default;
};

#endif // MATSOLVER_H
