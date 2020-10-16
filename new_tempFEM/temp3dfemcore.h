#ifndef TEMP3DFEMCORE_H
#define TEMP3DFEMCORE_H
#include "datatype.h"
#if !defined(ARMA_32BIT_WORD)
#define ARMA_32BIT_WORD
#endif
#include <armadillo>
#include "slu_ddefs.h"
#include <iomanip>
#include <omp.h>

class Temp3dfemcore
{
public:
    Temp3dfemcore(const char *fn);
    bool load3DFEMCOMSOL();
    void generateMetisMesh(int part);
    void preCalculation();
    void setCondition();
    void NRSolve();
    double TtoCond(double T);

private:
    //直接从COMSOL分网数据中读取到的信息
    const char *m_COMSOLMesh;
    int m_num_pts{0};
    int m_num_VtxEle{0};
    int m_num_EdgEle{0};
    int m_num_TriEle{0};
    int m_num_TetEle{0};
    C3DNode *mp_3DNode;
    CVtxElement *mp_VtxEle;
    CEdgElement *mp_EdgEle;
    CTriElement *mp_TriEle;
    CTetElement *mp_TetEle;
    const double Precision{1e-8};

    //区域分解相关信息

};

#endif // TEMP3DFEMCORE_H
