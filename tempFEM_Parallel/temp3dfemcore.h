#ifndef TEMP3DFEMCORE_H
#define TEMP3DFEMCORE_H
#include "datatype.h"
#if !defined(ARMA_32BIT_WORD)
#define ARMA_32BIT_WORD
#endif
#include <armadillo>
#include "matsolver.h"

#include <iomanip>
#include <omp.h>

using namespace arma;

class Temp3dfemcore
{
public:
    Temp3dfemcore();
    Temp3dfemcore(const char *fn);
    ~Temp3dfemcore();
    bool load3DFEMCOMSOL();
    void preCalculation();
    void setCondition(Modeltype TYPE);
    void NRSolve();
    double AirTtoCond(double T);
    double CopperTtoCond(double T);
    double IronTtoCond(double T);
    double Iron304TtoCond(double T);
    double KaptonTtoCond(double T);
    void outputResult();

    //DDTLM相关
    bool GenerateMetisMesh(int partition);

protected:
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
    char m_METISMesh[256];
    int m_num_part{0};
    int *m_tpartTable;    //保存三角形单元在第几个分区
    int *m_npartTable;    //保存节点在第几个分区

    int nprocs;

    MatSolver *solver;
};

#endif // TEMP3DFEMCORE_H
