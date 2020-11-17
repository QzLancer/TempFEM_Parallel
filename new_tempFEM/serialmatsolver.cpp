#include "serialmatsolver.h"

SerialMatSolver::SerialMatSolver()
{

}

double *SerialMatSolver::solveMatrix(umat locs, mat vals, vec F, int size)
{
    //将数据转化成列压缩存储形式
    options.ColPerm = COLAMD;
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

    set_default_options(&options);
    dCreate_CompCol_Matrix(&sluA, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
    Astore = (NCformat *)sluA.Store;
    //        printf("Dimension %dx%d; # nonzeros %d\n", sluA.nrow, sluA.ncol, Astore->nnz);

    nrhs = 1;
    if (!(rhs = doubleMalloc(m * nrhs))) ABORT("Malloc fails for rhs[].");
    //将内存拷贝过来
    //memmove(rhs, unknown_b, 5*sizeof(double));
    for (int i = 0; i < m; i++){
        rhs[i] = F(i);
    }
    dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);

    if (!(perm_c = intMalloc(n))) ABORT("Malloc fails for perm_c[].");
    if (!(perm_r = intMalloc(m))) ABORT("Malloc fails for perm_r[].");

    /* Initialize the statistics variables. */
    StatInit(&stat);

    dgssv(&options, &sluA, perm_c, perm_r, &L, &U, &B, &stat, &info);
    double *sol, *res;
    res = doubleMalloc(m);
    if (info == 0) {
        //            std::ofstream myTemp3D("../tempFEM/test/Temp3D.txt");
        /* This is how you could access the solution matrix. */
        sol = (double*)((DNformat*)B.Store)->nzval;
        //            myTemp3D.close();
    }else {
        cout << "info = " << info << endl;
    }

    for(int i = 0; i < size; ++i){
        res[i] = sol[i];
    }

    //        qDebug() << "Matrix solver finish.";

    return res;
}

double *SerialMatSolver::solveMatrix_LU(vec F)
{
    trans = NOTRANS;

    //将内存拷贝过来
    //memmove(rhs, unknown_b, 5*sizeof(double));
    for (int i = 0; i < m; i++){
        rhs[i] = F(i);
    }
    dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);

    double *sol, *res;
    res = doubleMalloc(m);
    dgstrs(trans, &L, &U, perm_r, perm_c, &B, &stat, &info);

    if (info == 0) {
        //            std::ofstream myTemp3D("../tempFEM/test/Temp3D.txt");
        /* This is how you could access the solution matrix. */
        sol = (double*)((DNformat*)B.Store)->nzval;
        //            myTemp3D.close();
    }else {
        cout << "info = " << info << endl;
    }

    for(int i = 0; i < m; ++i){
        res[i] = sol[i];
    }

    //        qDebug() << "Matrix solver finish.";

    return res;
}

SerialMatSolver::~SerialMatSolver()
{
    SUPERLU_FREE(rhs);
    //    SUPERLU_FREE(xact);
    SUPERLU_FREE(perm_r);
    SUPERLU_FREE(perm_c);
    //    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    Destroy_SuperMatrix_Store(&sluA);
    free(a);
    free(asub);
    free(xa);
}
