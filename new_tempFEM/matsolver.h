#ifndef MATSOLVER_H
#define MATSOLVER_H
#include "datatype.h"
#if !defined(ARMA_32BIT_WORD)
#define ARMA_32BIT_WORD
#endif
#include <armadillo>
#include "slu_ddefs.h"

class MatSolver
{
public:
    MatSolver();
};

#endif // MATSOLVER_H
