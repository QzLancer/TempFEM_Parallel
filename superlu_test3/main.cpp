#include "slu_mt_ddefs.h"
//#include "SuperLU_MT.h"

#include <iostream>

void c_bridge_pdgssv_(int_t *nprocs, int_t *n, int_t *nnz, int_t *nrhs,
         double *values, int_t *rowind, int_t *colptr,
         double *b, int_t *ldb, int_t *info);

using namespace std;

int main()
{
     double   *a, *rhs;
     double   s, u, p, e, r, l;
     int      *asub, *xa;
     int m, n, nnz, nrhs, i, nprocs, info;

     cout << "abc" << endl;
     /* Initialize matrix A. */
     m = n = 5;
     nnz = 12;
     s = 19.0; u = 21.0; p = 16.0; e = 5.0; r = 18.0; l = 12.0;
     a[0] = s; a[1] = l; a[2] = l; a[3] = u; a[4] = l; a[5] = l;
     a[6] = u; a[7] = p; a[8] = u; a[9] = e; a[10]= u; a[11]= r;
     asub[0] = 0; asub[1] = 1; asub[2] = 4; asub[3] = 1;
     asub[4] = 2; asub[5] = 4; asub[6] = 0; asub[7] = 2;
     asub[8] = 0; asub[9] = 3; asub[10]= 3; asub[11]= 4;
     xa[0] = 0; xa[1] = 3; xa[2] = 6; xa[3] = 8; xa[4] = 10; xa[5] = 12;

     nrhs = 1;
     for (i = 0; i < m; ++i) rhs[i] = 1.0;

     nprocs = 8;

     c_bridge_pdgssv_(&nprocs, &n, &nnz, &nrhs
                      , a, asub, xa
                     ,rhs, &m, &info);

    return 0;
}

void c_bridge_pdgssv_(int_t *nprocs, int_t *n, int_t *nnz, int_t *nrhs,
         double *values, int_t *rowind, int_t *colptr,
         double *b, int_t *ldb, int_t *info)
{
    SuperMatrix A, B, L, U;
    SCformat *Lstore;
    NCformat *Ustore;
    int_t      *perm_r; /* row permutations from partial pivoting */
    int_t      *perm_c; /* column permutation vector */
    int_t      panel_size, permc_spec;
    superlu_memusage_t superlu_memusage;


    dCreate_CompCol_Matrix(&A, *n, *n, *nnz, values, rowind, colptr,
               SLU_NC, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&B, *n, *nrhs, b, *ldb, SLU_DN, SLU_D, SLU_GE);

    if ( !(perm_r = intMalloc(*n)) ) SUPERLU_ABORT("Malloc fails for perm_r[].");
    if ( !(perm_c = intMalloc(*n)) ) SUPERLU_ABORT("Malloc fails for perm_c[].");

    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   permc_spec = 0: natural ordering
     *   permc_spec = 1: minimum degree ordering on structure of A'*A
     *   permc_spec = 2: minimum degree ordering on structure of A'+A
     *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
     */
    permc_spec = 1;
    get_perm_c(permc_spec, &A, perm_c);

    panel_size = sp_ienv(1);

    pdgssv(*nprocs, &A, perm_c, perm_r, &L, &U, &B, info);

    if ( *info == 0 ) {

    Lstore = (SCformat *) L.Store;
    Ustore = (NCformat *) U.Store;
        printf("#NZ in factor L = " IFMT "\n", Lstore->nnz);
        printf("#NZ in factor U = " IFMT "\n", Ustore->nnz);
        printf("#NZ in L+U = " IFMT "\n", Lstore->nnz + Ustore->nnz - L.ncol);

    superlu_dQuerySpace(*nprocs, &L, &U, panel_size, &superlu_memusage);
    printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions " IFMT "\n",
           superlu_memusage.for_lu/1e6, superlu_memusage.total_needed/1e6,
           superlu_memusage.expansions);

    } else {
    printf("dgssv() error returns INFO= " IFMT "\n", *info);
    if ( info <= n ) { /* factorization completes */
        superlu_dQuerySpace(*nprocs, &L, &U, panel_size, &superlu_memusage);
        printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions " IFMT "\n",
           superlu_memusage.for_lu/1e6,
           superlu_memusage.total_needed/1e6,
           superlu_memusage.expansions);
    }
    }

    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    Destroy_SuperMatrix_Store(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
}
