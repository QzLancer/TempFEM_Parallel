#ifndef SERIALMATSOLVER_H
#define SERIALMATSOLVER_H

#include "matsolver.h"
#include "slu_ddefs.h"

class SerialMatSolver : public MatSolver
{
public:
    SerialMatSolver();
    virtual double * solveMatrix(umat locs, mat vals, vec F, int size) override;
    virtual double * solveMatrix_LU(vec F) override;
    virtual ~SerialMatSolver() override;

private:
    SuperMatrix sluA;
    NCformat *Astore;
    double   *a;
    int      *asub, *xa;
    int      *perm_c; /* column permutation vector */
    int      *perm_r; /* row permutations from partial pivoting */
    SuperMatrix L;      /* factor L */
    SuperMatrix U;      /* factor U */
    SuperMatrix B;
    int      nrhs, ldx, info, m, n, nnz;
    double   *rhs;
    mem_usage_t   mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;
    trans_t trans;
};

#endif // SERIALMATSOLVER_H
