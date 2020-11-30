#ifndef PARALLELMATSOLVER_H
#define PARALLELMATSOLVER_H

#include "matsolver.h"
#include "slu_mt_ddefs.h"

class ParallelMatSolver : public MatSolver
{
public:
    ParallelMatSolver();
    virtual double * solveMatrix_LU(vec F) override;
    virtual double * solveMatrix(umat locs, mat vals, vec F, int size) override;
    virtual double * solveMatrix1(umat locs, mat vals, vec F, int size) override;
    virtual ~ParallelMatSolver() override;

private:
    double *a, *rhs;
    int *asub, *xa;
    int m, n, nnz, nrhs, nprocs, info;
    SuperMatrix A, B, L, U;
    int_t      *perm_r; /* row permutations from partial pivoting */
    int_t      *perm_c; /* column permutation vector */
    int_t      panel_size, permc_spec, relax;
    superlu_memusage_t superlu_memusage;
    Gstat_t  Gstat1;

    yes_no_t refact, usepr;
    trans_t trans;
    double diag_pivot_thresh, drop_tol;
    void *work;
    int lwork;
    double *utime;
    flops_t *ops, flopcnt;
    SuperMatrix AC;
    superlumt_options_t options;
    fact_t fact;

    //solve1新增变量
    SuperMatrix sluX;
    equed_t     equed;
    int_t       ldx;
    double      *rhsb, *rhsx, *xact;
    double      *R, *C;
    double      *ferr, *berr;
    double      u, rpg, rcond;
};

#endif // PARALLELMATSOLVER_H
