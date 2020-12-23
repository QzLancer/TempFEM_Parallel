#include "temp3dfemcore.h"
#include "parallelmatsolver.h"

#include <vector>
#include <map>
#include <omp.h>

#pragma execution_character_set("utf-8")

using namespace arma;

Temp3dfemcore::Temp3dfemcore()
{

}

Temp3dfemcore::Temp3dfemcore(const char *fn):
    m_COMSOLMesh(fn),
    mp_3DNode(nullptr),
    mp_VtxEle(nullptr),
    mp_EdgEle(nullptr),
    mp_TriEle(nullptr),
    mp_TetEle(nullptr),
    m_tpartTable(new int(0)),
    m_npartTable(new int(0)),
    nprocs(8),
    solver(new ParallelMatSolver)
{

}

Temp3dfemcore::~Temp3dfemcore()
{
//    SUPERLU_FREE(rhs);
//    //    SUPERLU_FREE(xact);
//    SUPERLU_FREE(perm_r);
//    SUPERLU_FREE(perm_c);
//    //    Destroy_CompCol_Matrix(&A);
//    Destroy_SuperMatrix_Store(&B);
//    Destroy_SuperNode_Matrix(&L);
//    Destroy_CompCol_Matrix(&U);
//    Destroy_SuperMatrix_Store(&sluA);
//    free(a);
//    free(asub);
//    free(xa);
    delete solver;
    free(m_tpartTable);
    free(m_npartTable);
    delete m_COMSOLMesh;
    delete [] mp_3DNode;
    delete [] mp_VtxEle;
    delete [] mp_TriEle;
    delete [] mp_TetEle;

}

bool Temp3dfemcore::load3DFEMCOMSOL()
{
    FILE *fp = nullptr;
    int err;
    char ch[256];
    err = fopen_s(&fp, m_COMSOLMesh, "r");
    if (!fp) {
        cout << "Error: openning file!" << endl;
        return 1;
    }
    //--------------Read the head-----------------------------
    for (int i = 0; i < 18; i++) {
        fgets(ch, 256, fp);
    }
    //-----------------mesh point-----------------------------
    if (fscanf_s(fp, "%d # number of mesh points\n", &m_num_pts) != 1) {
        cout << "Error: reading num_bdr_ns!" << endl;
        return 1;
    }
    else cout << m_num_pts << "number of mesh points." << endl;
    mp_3DNode = new C3DNode[m_num_pts];
    int pts_ind;//the beginning of the points index

    if (fscanf_s(fp, "%d # lowest mesh point index\n", &pts_ind) != 1) {
        cout << "Error: reading pts_ind!" << endl;
        return 1;
    }
    fgets(ch, 256, fp);

    for (int i = pts_ind; i < m_num_pts; i++) {
        //读取x,y坐标
        if (fscanf_s(fp, "%lf %lf %lf\n", &(mp_3DNode[i].x), &(mp_3DNode[i].y), &(mp_3DNode[i].z)) != 3) {
            cout << "Error: reading mesh point!" << endl;
            return 1;
        }
        //                else{
        //                    qDebug() << mp_3DNode[i].x << " " << mp_3DNode[i].y << " " << mp_3DNode[i].z << "\n";
        //                }
    }
    //---------------vertexnode-------------------------------
    for (int i = 0; i < 7; i++)
        fgets(ch, 256, fp);
    int num_vtx_ns;

    if (fscanf_s(fp, "%d # number of nodes per element\n", &num_vtx_ns) != 1) {
        cout << "Error: reading num_vtx_ns!" << endl;
        return 1;
    }

    if (fscanf_s(fp, "%d # number of elements\n", &m_num_VtxEle) != 1) {
        cout << "Error: reading m_num_VtxEle!" << endl;
        return 1;
    }
    else cout << m_num_VtxEle <<  "number of Vertex elements." << endl;
    fgets(ch, 256, fp);
    //    qDebug() << m_num_VtxEle;
    mp_VtxEle = new CVtxElement[m_num_VtxEle];
    for (int i = 0; i < m_num_VtxEle; i++) {
        if (fscanf_s(fp, "%d \n", &((mp_VtxEle + i)->n)) != 1) {
            cout << "Error: reading vertex element points!" << endl;
            return 1;
        }
    }
    //---------------vertexdomain-------------------------------
    for (int i = 0; i < 2; i++)
        fgets(ch, 256, fp);
    for (int i = 0; i < m_num_VtxEle; i++){
        if (fscanf_s(fp, "%d \n", &(mp_VtxEle[i].domain)) != 1){
            cout << "Error: reading vertex domain!" << endl;
            return 1;
        }
        else{
            mp_VtxEle[i].domain++;
        }
    }
    //----------------edgnode-----------------------------------
    for (int i = 0; i < 6; i++){
        fgets(ch, 256, fp);
    }
    if(fscanf_s(fp, "%d # number of elements\n", &m_num_EdgEle) != 1){
        cout <<  "Error: reading m_num_EdgEle" << endl;
        return 1;
    }
    else cout << m_num_EdgEle << "number of Edge elements." << endl;
    mp_EdgEle = new CEdgElement[m_num_EdgEle];
    fgets(ch, 256, fp);
    for (int i = 0; i < m_num_EdgEle; i++){
        if(fscanf_s(fp, "%d %d\n", &mp_EdgEle[i].n[0], &mp_EdgEle[i].n[1]) != 2){
            cout << "Error: reading edg element points!" << endl ;
            return 1;
        }
    }
    //----------------edgdomain---------------------------------
    for(int i = 0; i < 2; i++){
        fgets(ch, 256, fp);
    }
    for (int i = 0; i < m_num_EdgEle; i++){
        if(fscanf_s(fp, "%d \n", &mp_EdgEle[i].domain) != 1){
            cout << "Error: reading edgdomain!" << endl;
            return 1;
        }
        else{
            mp_EdgEle[i].domain++;
            //            qDebug() << "mp_EdgEle: " << mp_EdgEle[i].domain;
        }
    }
    //----------------trinode-----------------------------------
    for (int i = 0; i < 6; i++){
        fgets(ch, 256, fp);
    }
    if(fscanf_s(fp, "%d # number of elements\n", &m_num_TriEle) != 1){
        cout << "Error: reading m_num_TriEle!" << endl;
        return 1;
    }
    else cout << m_num_TriEle << "number of Triangle elements." << endl;
    fgets(ch, 256, fp);
    mp_TriEle = new CTriElement[m_num_TriEle];
    for (int i = 0; i < m_num_TriEle; i++){
        if (fscanf_s(fp, "%d %d %d \n", &mp_TriEle[i].n[0], &mp_TriEle[i].n[1], &mp_TriEle[i].n[2]) != 3) {
            cout << "Error: reading elements points!" << endl;
            return 1;
        }
    }
    //----------------tridomain---------------------------------
    for (int i = 0; i < 2; i++){
        fgets(ch, 256, fp);
    }
    for (int i = 0; i < m_num_TriEle; i++){
        if(fscanf_s(fp, "%d \n", &mp_TriEle[i].domain) != 1){
            cout << "Error: reading tridomain!" << endl;
            return 1;
        }
        else{
            mp_TriEle[i].domain++;
        }
    }
    //----------------tetnode-----------------------------------
    for (int i = 0; i < 6; i++){
        fgets(ch, 256, fp);
    }
    if(fscanf_s(fp, "%d # number of elements\n", &m_num_TetEle) != 1){
        cout << "Error: reading m_num_TriEle!" << endl;
        return 1;
    }
    else cout << m_num_TetEle << "number of Tetrahedron elements." << endl;
    fgets(ch, 256, fp);
    mp_TetEle = new CTetElement[m_num_TetEle];
    for (int i = 0; i < m_num_TetEle; i++){
        if (fscanf_s(fp, "%d %d %d %d\n", &mp_TetEle[i].n[0], &mp_TetEle[i].n[1], &mp_TetEle[i].n[2], &mp_TetEle[i].n[3]) != 4) {
            cout << "Error: reading Tet elements points!" << endl;
            return 1;
        }
    }
    //----------------tetdomain-----------------------------------

    for (int i = 0; i < 2; i++){
        fgets(ch, 256, fp);
    }
    for (int i = 0; i < m_num_TetEle; i++){
        if(fscanf_s(fp, "%d \n", &mp_TetEle[i].domain) != 1){
            cout << "Error: reading tetdomain!" << endl;
            return 1;
        }
        else{
            //            mp_TetEle[i].domain++;
        }
    }

    fclose(fp);

    return 0;

}


void Temp3dfemcore::preCalculation()
{
    //I:计算所有四面体单元中的p,q,r,s,volume

    omp_set_num_threads(nprocs);
#pragma omp parallel for
    for(int i = 0; i < m_num_TetEle; ++i){
        mat p0 = zeros<mat>(3,3);
        mat v0 = zeros<mat>(4,4);

        for(int j = 0; j < 4; ++j){
            mp_TetEle[i].x[j] = mp_3DNode[mp_TetEle[i].n[j]].x;
            mp_TetEle[i].y[j] = mp_3DNode[mp_TetEle[i].n[j]].y;
            mp_TetEle[i].z[j] = mp_3DNode[mp_TetEle[i].n[j]].z;
        }
        //求解p
        p0 = {{mp_TetEle[i].x[1], mp_TetEle[i].y[1], mp_TetEle[i].z[1]},
              {mp_TetEle[i].x[2], mp_TetEle[i].y[2], mp_TetEle[i].z[2]},
              {mp_TetEle[i].x[3], mp_TetEle[i].y[3], mp_TetEle[i].z[3]}};
        mp_TetEle[i].p[0] = det(p0);
        p0 = {{mp_TetEle[i].x[0], mp_TetEle[i].y[0], mp_TetEle[i].z[0]},
              {mp_TetEle[i].x[2], mp_TetEle[i].y[2], mp_TetEle[i].z[2]},
              {mp_TetEle[i].x[3], mp_TetEle[i].y[3], mp_TetEle[i].z[3]}};
        mp_TetEle[i].p[1] = -det(p0);
        p0 = {{mp_TetEle[i].x[0], mp_TetEle[i].y[0], mp_TetEle[i].z[0]},
              {mp_TetEle[i].x[1], mp_TetEle[i].y[1], mp_TetEle[i].z[1]},
              {mp_TetEle[i].x[3], mp_TetEle[i].y[3], mp_TetEle[i].z[3]}};
        mp_TetEle[i].p[2] = det(p0);
        p0 = {{mp_TetEle[i].x[0], mp_TetEle[i].y[0], mp_TetEle[i].z[0]},
              {mp_TetEle[i].x[1], mp_TetEle[i].y[1], mp_TetEle[i].z[1]},
              {mp_TetEle[i].x[2], mp_TetEle[i].y[2], mp_TetEle[i].z[2]}};
        mp_TetEle[i].p[3] = -det(p0);

        //求解q
        p0 = {{1, mp_TetEle[i].y[1], mp_TetEle[i].z[1]},
              {1, mp_TetEle[i].y[2], mp_TetEle[i].z[2]},
              {1, mp_TetEle[i].y[3], mp_TetEle[i].z[3]}};
        mp_TetEle[i].q[0] = -det(p0);
        p0 = {{1, mp_TetEle[i].y[0], mp_TetEle[i].z[0]},
              {1, mp_TetEle[i].y[2], mp_TetEle[i].z[2]},
              {1, mp_TetEle[i].y[3], mp_TetEle[i].z[3]}};
        mp_TetEle[i].q[1] = det(p0);
        p0 = {{1, mp_TetEle[i].y[0], mp_TetEle[i].z[0]},
              {1, mp_TetEle[i].y[1], mp_TetEle[i].z[1]},
              {1, mp_TetEle[i].y[3], mp_TetEle[i].z[3]}};
        mp_TetEle[i].q[2] = -det(p0);
        p0 = {{1, mp_TetEle[i].y[0], mp_TetEle[i].z[0]},
              {1, mp_TetEle[i].y[1], mp_TetEle[i].z[1]},
              {1, mp_TetEle[i].y[2], mp_TetEle[i].z[2]}};
        mp_TetEle[i].q[3] = det(p0);

        //求解r
        p0 = {{mp_TetEle[i].x[1], 1, mp_TetEle[i].z[1]},
              {mp_TetEle[i].x[2], 1, mp_TetEle[i].z[2]},
              {mp_TetEle[i].x[3], 1, mp_TetEle[i].z[3]}};
        mp_TetEle[i].r[0] = -det(p0);
        p0 = {{mp_TetEle[i].x[0], 1, mp_TetEle[i].z[0]},
              {mp_TetEle[i].x[2], 1, mp_TetEle[i].z[2]},
              {mp_TetEle[i].x[3], 1, mp_TetEle[i].z[3]}};
        mp_TetEle[i].r[1] = det(p0);
        p0 = {{mp_TetEle[i].x[0], 1, mp_TetEle[i].z[0]},
              {mp_TetEle[i].x[1], 1, mp_TetEle[i].z[1]},
              {mp_TetEle[i].x[3], 1, mp_TetEle[i].z[3]}};
        mp_TetEle[i].r[2] = -det(p0);
        p0 = {{mp_TetEle[i].x[0], 1, mp_TetEle[i].z[0]},
              {mp_TetEle[i].x[1], 1, mp_TetEle[i].z[1]},
              {mp_TetEle[i].x[2], 1, mp_TetEle[i].z[2]}};
        mp_TetEle[i].r[3] = det(p0);

        //求解s
        p0 = {{mp_TetEle[i].x[1], mp_TetEle[i].y[1], 1},
              {mp_TetEle[i].x[2], mp_TetEle[i].y[2], 1},
              {mp_TetEle[i].x[3], mp_TetEle[i].y[3], 1}};
        mp_TetEle[i].s[0] = -det(p0);
        p0 = {{mp_TetEle[i].x[0], mp_TetEle[i].y[0], 1},
              {mp_TetEle[i].x[2], mp_TetEle[i].y[2], 1},
              {mp_TetEle[i].x[3], mp_TetEle[i].y[3], 1}};
        mp_TetEle[i].s[1] = det(p0);
        p0 = {{mp_TetEle[i].x[0], mp_TetEle[i].y[0], 1},
              {mp_TetEle[i].x[1], mp_TetEle[i].y[1], 1},
              {mp_TetEle[i].x[3], mp_TetEle[i].y[3], 1}};
        mp_TetEle[i].s[2] = -det(p0);
        p0 = {{mp_TetEle[i].x[0], mp_TetEle[i].y[0], 1},
              {mp_TetEle[i].x[1], mp_TetEle[i].y[1], 1},
              {mp_TetEle[i].x[2], mp_TetEle[i].y[2], 1}};
        mp_TetEle[i].s[3] = det(p0);

        //求解volume
        v0 = {{1, mp_TetEle[i].x[0], mp_TetEle[i].y[0], mp_TetEle[i].z[0]},
              {1, mp_TetEle[i].x[1], mp_TetEle[i].y[1], mp_TetEle[i].z[1]},
              {1, mp_TetEle[i].x[2], mp_TetEle[i].y[2], mp_TetEle[i].z[2]},
              {1, mp_TetEle[i].x[3], mp_TetEle[i].y[3], mp_TetEle[i].z[3]}};
        mp_TetEle[i].Volume = det(v0)/6;
    }

    //II:计算所有三角形单元的Area
#pragma omp parallel for
    for(int i = 0; i < m_num_TriEle; ++i){
        for(int j = 0; j < 3; ++j){
            mp_TriEle[i].x[j] = mp_3DNode[mp_TriEle[i].n[j]].x;
            mp_TriEle[i].y[j] = mp_3DNode[mp_TriEle[i].n[j]].y;
            mp_TriEle[i].z[j] = mp_3DNode[mp_TriEle[i].n[j]].z;
            //                qDebug() << "mp_TriEle " << i << " x " << j << " = " << mp_TriEle[i].x[j];
            //                qDebug() << "mp_TriEle " << i << " y " << j << " = " << mp_TriEle[i].y[j];
            //                qDebug() << "mp_TriEle " << i << " z " << j << " = " << mp_TriEle[i].z[j];
        }
        vec ab = zeros<vec>(3);
        vec ac = zeros<vec>(3);
        vec s0 = zeros<vec>(3);
        ab = {mp_TriEle[i].x[1]-mp_TriEle[i].x[0], mp_TriEle[i].y[1]-mp_TriEle[i].y[0], mp_TriEle[i].z[1]-mp_TriEle[i].z[0]};
        ac = {mp_TriEle[i].x[2]-mp_TriEle[i].x[0], mp_TriEle[i].y[2]-mp_TriEle[i].y[0], mp_TriEle[i].z[2]-mp_TriEle[i].z[0]};
        s0 = cross(ab, ac);
        mp_TriEle[i].Area = sqrt(s0(0)*s0(0) + s0(1)*s0(1) + s0(2)*s0(2))/2;
        //            qDebug() << "mp_TriEle" << i << "area = " << mp_TriEle[i].Area;
        //           ab.print("ab = ");
        //           ac.print("ac = ");
    }
}

void Temp3dfemcore::setCondition(Modeltype model)
{
    switch (model) {
    case MODEL1:
        //热源设置
        //        std::ofstream mytetsource("../tempFEM/test/tetsource.txt");
        for(int i = 0; i < m_num_TetEle; ++i){
            if((mp_TetEle[i].domain == 7) | (mp_TetEle[i].domain == 9)){
                mp_TetEle[i].source = 500000;
            }else{
                mp_TetEle[i].source = 0;
            }
            //         mytetsource << "mp_TetEle " << i << " domain = " << mp_TetEle[i].source << endl;
        }
        //热导率设置
        for(int i = 0; i < m_num_TetEle; ++i){
            if((mp_TetEle[i].domain == 7) | (mp_TetEle[i].domain == 9)){
                mp_TetEle[i].cond = 386.6;
                mp_TetEle[i].Material = Copper;
                mp_TetEle[i].LinearFlag = 1;
            }
            else if((mp_TetEle[i].domain == 2) | (mp_TetEle[i].domain == 4) | (mp_TetEle[i].domain == 10) | (mp_TetEle[i].domain == 12)| (mp_TetEle[i].domain == 13)){
                mp_TetEle[i].cond = 80.32;
                mp_TetEle[i].Material = Iron;
                mp_TetEle[i].LinearFlag = 1;
            }
            else if((mp_TetEle[i].domain == 1) | (mp_TetEle[i].domain == 6) | (mp_TetEle[i].domain == 8)){
                mp_TetEle[i].cond = 0.26;
                mp_TetEle[i].Material = Nylon;
                mp_TetEle[i].LinearFlag = 1;
            }
            else if((mp_TetEle[i].domain == 3) | (mp_TetEle[i].domain == 5) | (mp_TetEle[i].domain == 11)){
    //            mp_TetEle[i].cond = 0.03;   //初始猜测，用于Y0的计算
                mp_TetEle[i].cond = 0.01;
                mp_TetEle[i].Material = Air;
                mp_TetEle[i].LinearFlag = 0;
            }
        }
        //第三类边界条件设置
        for(int i = 0; i < m_num_TriEle; ++i){
            if((mp_TriEle[i].domain == 1) | (mp_TriEle[i].domain == 2) | (mp_TriEle[i].domain == 3) | (mp_TriEle[i].domain == 4) |
                    (mp_TriEle[i].domain == 5) | (mp_TriEle[i].domain == 6) | (mp_TriEle[i].domain == 91) | (mp_TriEle[i].domain == 92) |
                    (mp_TriEle[i].domain == 93) | (mp_TriEle[i].domain == 137) | (mp_TriEle[i].domain == 143) | (mp_TriEle[i].domain == 183)){
                mp_TriEle[i].bdr = 3;
                mp_TriEle[i].h = 20;
                mp_TriEle[i].Text = 293.15;
            }
        }
        break;

    case MODEL2:
        //热源设置
        //        std::ofstream mytetsource("../tempFEM/test/tetsource.txt");
        for(int i = 0; i < m_num_TetEle; ++i){
            if(mp_TetEle[i].domain == 3){
                mp_TetEle[i].source = 200000;
            }else if(mp_TetEle[i].domain == 12 || mp_TetEle[i].domain == 15){
                mp_TetEle[i].source = 100000;
            }else{
                mp_TetEle[i].source = 0;
            }
            //         mytetsource << "mp_TetEle " << i << " domain = " << mp_TetEle[i].source << endl;
        }
        //热导率设置
        for(int i = 0; i < m_num_TetEle; ++i){
            if((mp_TetEle[i].domain == 3) | (mp_TetEle[i].domain == 9) | (mp_TetEle[i].domain == 10) | (mp_TetEle[i].domain == 11) |
                    (mp_TetEle[i].domain == 12) | (mp_TetEle[i].domain == 13) | (mp_TetEle[i].domain == 15) | (mp_TetEle[i].domain == 16) |
                    (mp_TetEle[i].domain == 17)){
                mp_TetEle[i].cond = 386.6;
                mp_TetEle[i].Material = Copper;
                mp_TetEle[i].LinearFlag = 1;
            }

//            if((mp_TetEle[i].domain == 9) | (mp_TetEle[i].domain == 10) | (mp_TetEle[i].domain == 11) |
//                    (mp_TetEle[i].domain == 13) | (mp_TetEle[i].domain == 16) |
//                    (mp_TetEle[i].domain == 17)){
//                mp_TetEle[i].cond = 386.6;
//                mp_TetEle[i].Material = Iron;
//                mp_TetEle[i].LinearFlag = 0;
//            }
//            else if(mp_TetEle[i].domain == 3 | mp_TetEle[i].domain == 12 | mp_TetEle[i].domain == 15){
//                mp_TetEle[i].cond = 386.6;
//                mp_TetEle[i].Material = Iron;
//                mp_TetEle[i].LinearFlag = 1;
//            }

            else if((mp_TetEle[i].domain == 1) | (mp_TetEle[i].domain == 6) | (mp_TetEle[i].domain == 8)){
                mp_TetEle[i].cond = 76.2;
                mp_TetEle[i].Material = Iron;
                mp_TetEle[i].LinearFlag = 0;
            }
            else if((mp_TetEle[i].domain == 18) | (mp_TetEle[i].domain == 21) | (mp_TetEle[i].domain == 22)){
                mp_TetEle[i].cond = 15.33;
                mp_TetEle[i].Material = Iron304;
                mp_TetEle[i].LinearFlag = 0;
            }
            else if((mp_TetEle[i].domain == 4) | (mp_TetEle[i].domain == 20)){
                mp_TetEle[i].cond = 16.7;
                mp_TetEle[i].Material = Ceramics;
                mp_TetEle[i].LinearFlag = 1;
            }
            else if(mp_TetEle[i].domain == 2){
                mp_TetEle[i].cond = 1.76;
                mp_TetEle[i].Material = Kapton;
                mp_TetEle[i].LinearFlag = 0;
            }
            else if((mp_TetEle[i].domain == 5) | (mp_TetEle[i].domain == 7) | (mp_TetEle[i].domain == 14) | (mp_TetEle[i].domain == 19) |
                        (mp_TetEle[i].domain == 20)){
    //            mp_TetEle[i].cond = 0.03;   //初始猜测，用于Y0的计算
                mp_TetEle[i].cond = 0.01;
                mp_TetEle[i].Material = Air;
                mp_TetEle[i].LinearFlag = 0;
            }
        }
        //第三类边界条件设置
        for(int i = 0; i < m_num_TriEle; ++i){
            if((mp_TriEle[i].domain == 1) | (mp_TriEle[i].domain == 2) | (mp_TriEle[i].domain == 3) | (mp_TriEle[i].domain == 4) |
                    (mp_TriEle[i].domain == 5) | (mp_TriEle[i].domain == 6) | (mp_TriEle[i].domain == 7) | (mp_TriEle[i].domain == 8) |
                    (mp_TriEle[i].domain == 9) | (mp_TriEle[i].domain == 10) | (mp_TriEle[i].domain == 29) | (mp_TriEle[i].domain == 30) |
                    (mp_TriEle[i].domain == 63) | (mp_TriEle[i].domain == 64) | (mp_TriEle[i].domain == 67) | (mp_TriEle[i].domain == 68) |
                    (mp_TriEle[i].domain == 77) | (mp_TriEle[i].domain == 78) | (mp_TriEle[i].domain == 81) | (mp_TriEle[i].domain == 82) |
                    (mp_TriEle[i].domain == 138) | (mp_TriEle[i].domain == 139) | (mp_TriEle[i].domain == 162) | (mp_TriEle[i].domain == 163) |
                    (mp_TriEle[i].domain == 164) | (mp_TriEle[i].domain == 165) | (mp_TriEle[i].domain == 166) | (mp_TriEle[i].domain == 176) |
                    (mp_TriEle[i].domain == 181) | (mp_TriEle[i].domain == 183) | (mp_TriEle[i].domain == 198) | (mp_TriEle[i].domain == 216) |
                    (mp_TriEle[i].domain == 233) | (mp_TriEle[i].domain == 245) | (mp_TriEle[i].domain == 254) | (mp_TriEle[i].domain == 286) |
                    (mp_TriEle[i].domain == 288) | (mp_TriEle[i].domain == 314) | (mp_TriEle[i].domain == 320) | (mp_TriEle[i].domain == 324) |
                    (mp_TriEle[i].domain == 325) | (mp_TriEle[i].domain == 329) | (mp_TriEle[i].domain == 330) | (mp_TriEle[i].domain == 331)){
                mp_TriEle[i].bdr = 3;
                mp_TriEle[i].h = 25;
                mp_TriEle[i].Text = 293.15;
            }
        }
        break;
    }

}

void Temp3dfemcore::NRSolve()
{
    const double PI = 3.14159265358979323846;
    int pos = 0;

    //1.确定第三类边界条件三角形单元数量
    int numbdr = 0;
    for(int i = 0; i < m_num_TriEle; i++){
        if(mp_TriEle[i].bdr == 3) numbdr++;
    }
    cout << "numbdr = " << numbdr << endl;
    //2.确定线性四面体单元和非线性四面体单元的数量
    int num_linear = 0;
    int num_nonlinear = 0;
    for(int i = 0; i < m_num_TetEle; ++i){
        if(mp_TetEle[i].LinearFlag == 1){
            ++num_linear;
        }else{
            ++num_nonlinear;
        }
    }
    cout << "number of linear element = " << num_linear << endl;
    cout << "number of nonlinear element = " << num_nonlinear << endl;
    //3.求解每个单元的系数矩阵
    CTetResistMatrix* TetResist = (CTetResistMatrix*)malloc(m_num_TetEle * sizeof(CTetResistMatrix));
    for(int k = 0; k < m_num_TetEle; ++k){
        for(int i = 0; i < 4; ++i){
            for(int j = 0; j < 4; ++j){
                TetResist[k].C[i][j] = (mp_TetEle[k].q[i]*mp_TetEle[k].q[j]+mp_TetEle[k].r[i]*mp_TetEle[k].r[j]+mp_TetEle[k].s[i]*mp_TetEle[k].s[j])/(36*mp_TetEle[k].Volume);
            }
        }
    }
    //4.线性单元装配过程
    umat locs(2, 16*m_num_TetEle+9*numbdr);
    mat vals(1, 16*m_num_TetEle+9*numbdr);
    vec F = zeros<vec>(m_num_pts);
    double *Va = (double*)calloc(m_num_pts, sizeof(double));
    double *Va_old = (double*)calloc(m_num_pts, sizeof(double));
    for(int i = 0; i < m_num_pts; ++i){
        Va[i] = 273.15;
        Va_old[i] = 0;
    }
    double St;
    double Ft;
    double Se;
    double Fe;
    //四面体单元装配过程
    for(int k = 0; k < m_num_TetEle; ++k){
        for(int i = 0; i < 4; ++i){
            for(int j = 0; j < 4; ++j){
                if(mp_TetEle[k].LinearFlag == 1){
                    St = mp_TetEle[k].cond*TetResist[k].C[i][j];
                    locs(0, pos) = mp_TetEle[k].n[i];
                    locs(1, pos) = mp_TetEle[k].n[j];
                    vals(0, pos) = St;
                    ++pos;
                }
            }
            Ft = mp_TetEle[k].source*mp_TetEle[k].Volume/4;
            F(mp_TetEle[k].n[i]) = F(mp_TetEle[k].n[i]) + Ft;
        }
    }

    //三角形单元装配
    for(int k = 0; k < m_num_TriEle; ++k){
        for(int i = 0; i < 3; ++i){
            if(mp_TriEle[k].bdr == 3){
                for(int j = 0; j < 3; ++j){
                    if(i == j){
                        Se = mp_TriEle[k].h*mp_TriEle[k].Area/6;
                    }else{
                        Se = mp_TriEle[k].h*mp_TriEle[k].Area/12;
                    }
                    locs(0, pos) = mp_TriEle[k].n[i];
                    locs(1, pos) = mp_TriEle[k].n[j];
                    vals(0, pos) = Se;
                    ++pos;
                }
                Fe = mp_TriEle[k].h*mp_TriEle[k].Text*mp_TriEle[k].Area/3;
                F(mp_TriEle[k].n[i]) = F(mp_TriEle[k].n[i]) + Fe;
            }

        }
    }

    //5.迭代过程
    double t1, t2;
    double time;
    int MAX_ITER = 200;
    int pos1 = pos;
    vec F1 = F;
    for(int iter = 0; iter < MAX_ITER; ++iter){
        t1 = SuperLU_timer_();
        pos = pos1;
        F = F1;
        for(int k = 0; k < m_num_TetEle; ++k){
            if(mp_TetEle[k].LinearFlag == 0){
                int n[4];
                for(int i = 0; i < 4; ++i){
                    n[i] = mp_TetEle[k].n[i];
                }
                double T = (Va[n[0]]+Va[n[1]]+Va[n[2]]+Va[n[3]])/4;
                //根据材料的类型计算热导率
                double Cond;
                double CondPartialT;
                switch(mp_TetEle[k].Material){
                case(Copper):
                    Cond = CopperTtoCond(T);
                    CondPartialT = (CopperTtoCond(T+0.01)-CopperTtoCond(T))/0.01;
                    break;
                case(Iron):
                    Cond = IronTtoCond(T);
                    CondPartialT = (IronTtoCond(T+0.01)-IronTtoCond(T))/0.01;
                    break;
                case(Air):
                    Cond = AirTtoCond(T);
                    CondPartialT = (AirTtoCond(T+0.01)-AirTtoCond(T))/0.01;
                    break;
                case(Iron304):
                    Cond = Iron304TtoCond(T);
                    CondPartialT = (Iron304TtoCond(T+0.01)-Iron304TtoCond(T))/0.01;
                    break;
                case(Kapton):
                    Cond = KaptonTtoCond(T);
                    CondPartialT = (KaptonTtoCond(T+0.01)-KaptonTtoCond(T))/0.01;
                    break;
                }
//                                cout << "Cond = " << Cond << endl;
                //                qDebug() << "CondPartialT = " << CondPartialT;
                for(int i = 0; i < 4; ++i){
                    double FJ = 0;
                    double J1 = CondPartialT*(TetResist[k].C[i][0]*Va[n[0]]+TetResist[k].C[i][1]*Va[n[1]]+TetResist[k].C[i][2]*Va[n[2]]+TetResist[k].C[i][3]*Va[n[3]]);
                    for(int j = 0; j < 4; ++j){
                        St = Cond*TetResist[k].C[i][j] + J1;
                        locs(0, pos) = mp_TetEle[k].n[i];
                        locs(1, pos) = mp_TetEle[k].n[j];
                        vals(0, pos) = St;
                        ++pos;
                        FJ = FJ+J1*(Va[n[j]]);
                    }
                    F(n[i]) = F(n[i]) + FJ;
                }
            }
        }

        t2 = SuperLU_timer_();

        cout << "Jacobi Assemble time: " << t2 - t1 << endl;

//        Va = solveMatrix(locs, vals, F, m_num_pts);
        Va = solver->solveMatrix(locs, vals, F, m_num_pts);
        //6.判断收敛性
        double inner_error = 1;
        double a0 = 0, b = 0;
        for(int i = 0; i < m_num_pts; ++i){
            a0 += (Va_old[i] - Va[i])*(Va_old[i] - Va[i]);
            b += Va[i] * Va[i];
        }
        cout << "a0 = " << a0 << endl;
        cout << "b = " << b << endl;
        inner_error = sqrt(a0)/sqrt(b);
        cout << "inner_error = " << inner_error << endl;
        if(inner_error > Precision){
            t2 = SuperLU_timer_();
            time = t2-t1;
            cout << "iter step " << iter << " time = " << time << endl;
            for(int i = 0; i < m_num_pts; ++i){
                Va_old[i] = Va[i];
            }
        }else{
            t2 = SuperLU_timer_();
            time = t2-t1;
            cout << "iter step " << iter << " time = " << time << endl;
            for(int i = 0; i < m_num_pts; ++i){
                mp_3DNode[i].V = Va[i];
            }
            break;
        }
    }
    ///////////输出结果
    char fpath[256];
    sprintf(fpath,"../result/Temp3DNR_%d.txt",m_num_TetEle);
    std::ofstream mytemp(fpath);
    //    double temp[15076];
    mytemp << "x,y,z,Temp(K)" << endl;
    for(int i = 0; i < m_num_pts; i++){
        mytemp << mp_3DNode[i].x << "," << mp_3DNode[i].y << "," << mp_3DNode[i].z << "," << Va[i] << endl;
    }

    cout << "Ok." << endl;
}

double Temp3dfemcore::AirTtoCond(double T)
{
    return -0.00227583562+1.15480022e-4*T-7.90252856e-8*T*T+4.11702505e-11*T*T*T-7.43864331e-15*T*T*T*T;

}

double Temp3dfemcore::CopperTtoCond(double T)
{
    if(T >= 1 && T < 40){
        return 12.55868+36.66487*T+1.387207*T*T-0.07168113*T*T*T+6.99799E-4*T*T*T*T;
    }else if(T >= 40 && T < 70){
        return 2174.919-45.25448*T+0.3738471*T*T-9.504397E-4*T*T*T;
    }else if(T >= 70 && T < 100){
        return 2545.87-67.53869*T+0.8176488*T*T-0.004470238*T*T*T+9.22619E-6*T*T*T*T;
    }else if(T >= 100 && T < 300){
        return 555.4-2.116905*T+0.008971429*T*T-1.266667E-5*T*T*T;
    }

    return 423.7411-0.3133575*T+0.001013916*T*T-1.570451E-6*T*T*T+1.06222E-9*T*T*T*T-2.64198E-13*T*T*T*T*T;

}

double Temp3dfemcore::IronTtoCond(double T)
{
    if(T >= 1 && T < 20){
        return 161.9832*T+4.851772*T*T-0.7906936*T*T*T+0.01675048*T*T*T*T;
    }else if(T >= 20 && T < 68.5){
        return 2218.339+23.88314*T-5.081026*T*T+0.1373453*T*T*T-0.001490064*T*T*T*T+5.882957E-6*T*T*T*T*T;
    }else if(T >= 68.5 && T < 155.0){
        return 1682.9-46.91385*T+0.5423243*T*T-0.002845324*T*T*T+5.647112E-6*T*T*T*T;
    }else if(T >= 155.0 && T < 1183.0){
        return 155.6355-0.4470943*T+8.9922E-4*T*T-9.859391E-7*T*T*T+4.862049E-10*T*T*T*T-7.562752E-14*T*T*T*T*T;
    }

    return -5.097463+0.03983086*T-9.936321E-6*T*T;

}

double Temp3dfemcore::Iron304TtoCond(double T)
{
    if(T >= 1 && T < 45){
        return -0.03740871+0.06460546*T+0.003720604*T*T-8.390067E-5*T*T*T+6.006594E-7*T*T*T*T;
    }else if(T >= 45 && T < 293){
        return -1.031521+0.1813807*T-0.001088656*T*T+3.411681E-6*T*T*T-3.988389E-9*T*T*T*T;
    }

    return 6.742253+0.02864915*T;

}

double Temp3dfemcore::KaptonTtoCond(double T)
{
    if(T >= 1 && T < 45){
        return -0.001372384+0.005601653*T+2.082966E-6*T*T-5.05445E-9*T*T*T;
    }

    return -0.007707532+0.005769136*T+5.622796E-7*T*T-4.329984E-10*T*T*T;

}
