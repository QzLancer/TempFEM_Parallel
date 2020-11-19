#include "parallelmatsolver.h"

ParallelMatSolver::ParallelMatSolver()
{

}

double *ParallelMatSolver::solveMatrix_LU(vec F)
{
    for (int i = 0; i < m; i++){
        rhs[i] = F(i);
    }
    dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);

    dgstrs(trans, &L, &U, perm_r, perm_c, &B, &Gstat1, &info);
    double *sol, *res;
    res = doubleMalloc(m);
    if (info == 0 || info == n + 1) {
        sol = (double*)((DNformat*)B.Store)->nzval;
    } else if (info > 0 && lwork == -1) {
        printf("dgssv() error returns INFO= " IFMT "\n", info);
        if ( info <= n ) { /* factorization completes */
            superlu_dQuerySpace(nprocs, &L, &U, panel_size, &superlu_memusage);
            printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions " IFMT "\n",
           superlu_memusage.for_lu/1e6,
           superlu_memusage.total_needed/1e6,
           superlu_memusage.expansions);
        }
    }

    for(int i = 0; i < n; ++i){
        res[i] = sol[i];
    }
    return res;
}

double *ParallelMatSolver::solveMatrix(umat locs, mat vals, vec F, int size)
{
    nprocs = 8;

    //将数据转化成列压缩存储形式
    m = size; n = size;
    std::map<int, std::map<int, double>> mapmatrix; //mapmatrix[列编号][行编号][值]
    int row, col;
    double val;
    for(int i = 0; i < vals.size(); ++i){
        row = locs(0,i);
        col = locs(1,i);
        val = vals(i);
        if(mapmatrix.count(col) == 0){
            std::map<int, double> temp;
            temp[row] = val;
            mapmatrix[col] = temp;
        }else{
            if(mapmatrix[col].count(row) == 0){
                mapmatrix[col][row] = val;
            }else{
                mapmatrix[col][row] += val;
            }
        }
    }

    nnz = 0;
    xa = intMalloc(size+1);
    xa[0] = 0;
    for(std::map<int, std::map<int, double>>::iterator m = mapmatrix.begin(); m != mapmatrix.end(); ++m){
        nnz += m->second.size();
        xa[(m->first)+1] = nnz;
    }
    asub = intMalloc(nnz);
    a = doubleMalloc(nnz);
    int i = 0;
    for(std::map<int, std::map<int, double>>::iterator m = mapmatrix.begin(); m != mapmatrix.end(); ++m){
        for(std::map<int, double>::iterator n = m->second.begin(); n != m->second.end(); ++n){
            asub[i] = n->first;
            a[i] = n->second;
//            cout << "a:" << i << " = " << a[i] << endl;
            ++i;
        }
    }

    dCreate_CompCol_Matrix(&A, n, n, nnz, a, asub, xa,SLU_NC, SLU_D, SLU_GE);
    //将内存拷贝过来
    //memmove(rhs, unknown_b, 5*sizeof(double));
    nrhs = 1;
    if (!(rhs = doubleMalloc(n * nrhs))) SUPERLU_ABORT("Malloc fails for rhs[].");
    for (int i = 0; i < n; i++){
        rhs[i] = F(i);
    }
    dCreate_Dense_Matrix(&B, n, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);

    if ( !(perm_r = intMalloc(n)) ) SUPERLU_ABORT("Malloc fails for perm_r[].");
    if ( !(perm_c = intMalloc(n)) ) SUPERLU_ABORT("Malloc fails for perm_c[].");

    permc_spec = 1;
    get_perm_c(permc_spec, &A, perm_c);

    refact = NO;
    trans = NOTRANS;
    panel_size = sp_ienv(1);
    relax = sp_ienv(3);
    diag_pivot_thresh = 1.0;
    usepr = NO;
    drop_tol = 0.0;
    work = NULL;
    lwork = 0;
//    lwork = 10000000000000;
//    work = SUPERLU_MALLOC(lwork);
//    printf("Use work space of size LWORK = " IFMT " bytes\n", lwork);
//    if ( !work ) {
//        SUPERLU_ABORT("DLINSOLX: cannot allocate work[]");
//    }

    fact = DOFACT;

    n = A.ncol;
    StatAlloc(n, nprocs, panel_size, relax, &Gstat1);
    StatInit(n, nprocs, &Gstat1);
    utime = Gstat1.utime;
    ops = Gstat1.ops;

    pdgstrf_init(nprocs, fact, trans, refact, panel_size, relax, diag_pivot_thresh, usepr, drop_tol, perm_c, perm_r, work, lwork, &A, &AC, &options, &Gstat1);

    pdgstrf(&options, &AC, perm_r, &L, &U, &Gstat1, &info);


    pdgssv(nprocs, &A, perm_c, perm_r, &L, &U, &B, &info);
    double *sol, *res;
    res = doubleMalloc(m);
    if ( info == 0 ) {
        sol = (double*)((DNformat*)B.Store)->nzval;
    } else {
        printf("dgssv() error returns INFO= " IFMT "\n", info);
        if ( info <= n ) { /* factorization completes */
            superlu_dQuerySpace(nprocs, &L, &U, panel_size, &superlu_memusage);
            printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions " IFMT "\n",
           superlu_memusage.for_lu/1e6,
           superlu_memusage.total_needed/1e6,
           superlu_memusage.expansions);
        }
    }
    for(int i = 0; i < size; ++i){
        res[i] = sol[i];
    }
    return res;
}

ParallelMatSolver::~ParallelMatSolver()
{
    SUPERLU_FREE(rhs);
    //    SUPERLU_FREE(xact);
    SUPERLU_FREE(perm_r);
    SUPERLU_FREE(perm_c);
    //    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    Destroy_SuperMatrix_Store(&A);
    free(a);
    free(asub);
    free(xa);
}
