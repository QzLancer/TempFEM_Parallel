#include "tlmcore.h"

#include <QDebug>

TLMCore::TLMCore(const char *fn):
    Temp3dfemcore(fn)
{

}

TLMCore::~TLMCore()
{

}

void TLMCore::TLMSolve1()
{
    const double PI = 3.14159265358979323846;
    int pos = 0;

    //1.确定第三类边界条件三角形单元数量
    int numbdr = 0;
    for(int i = 0; i < m_num_TriEle; i++){
        if(mp_TriEle[i].bdr == 3) numbdr++;
    }
    qDebug() << "numbdr = " << numbdr;

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
    qDebug() << "number of linear element = " << num_linear;
    qDebug() << "number of nonlinear element = " << num_nonlinear;

    //3.求解每个单元的系数矩阵
    CTetResistMatrix* TetResist = (CTetResistMatrix*)malloc(m_num_TetEle * sizeof(CTetResistMatrix));
    CTetConnanceMatrix* TetY0 = (CTetConnanceMatrix*)malloc(m_num_TetEle * sizeof(CTetConnanceMatrix));
    for(int k = 0; k < m_num_TetEle; ++k){
        for(int i = 0; i < 4; ++i){
            for(int j = 0; j < 4; ++j){
                TetResist[k].C[i][j] = (mp_TetEle[k].q[i]*mp_TetEle[k].q[j]+mp_TetEle[k].r[i]*mp_TetEle[k].r[j]+mp_TetEle[k].s[i]*mp_TetEle[k].s[j])/(36*mp_TetEle[k].Volume);
                TetY0[k].Y0[i][j] = -mp_TetEle[i].cond * TetResist[k].C[i][j];
//                TetY0[k].Y0[i][j] = 0.0001;
            }
        }
    }

    //4.初始化单元入射/反射电压矩阵
    CTetVoltageMatrix* TetVi = (CTetVoltageMatrix*)malloc(m_num_TetEle * sizeof(CTetVoltageMatrix));
    for(int k = 0; k < m_num_TetEle; ++k){
        for(int i = 0; i < 4; ++i){
            for(int j = 0; j < 4; ++j){
                TetVi[k].Vi[i][j] = 0;
            }
        }
    }

    //5.线性单元装配
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

    //6. 非线性单元左侧系数矩阵装配
    double Yt;
    double Y0 = 0.3;
    for(int k = 0; k < m_num_TetEle; ++k){
        if(mp_TetEle[k].LinearFlag == 0){
            for(int i = 0; i < 4; ++i){
                double Yii = 0;
                for(int j = 0; j < 4; ++j){
//                    Yt = -TetY0[k].Y0[i][j];
//                    locs(0, pos) = mp_TetEle[k].n[i];
//                    locs(1, pos) = mp_TetEle[k].n[j];
//                    vals(0, pos) = Yt;
//                    ++pos;

                    if(i != j){
                        if((TetResist[k].C[i][j] < 0)){
                            Yt = -TetY0[k].Y0[i][j];
                            Yii -= Yt;
                            locs(0, pos) = mp_TetEle[k].n[i];
                            locs(1, pos) = mp_TetEle[k].n[j];
                            vals(0, pos) = Yt;
                            ++pos;
                        }
                        else{
                            locs(0, pos) = mp_TetEle[k].n[i];
                            locs(1, pos) = mp_TetEle[k].n[j];
                            vals(0, pos) = 0;
                            ++pos;
                        }
                    }
                }
                locs(0, pos) = mp_TetEle[k].n[i];
                locs(1, pos) = mp_TetEle[k].n[i];
                vals(0, pos) = Yii;
                ++pos;
            }
        }
    }

    //7.第一次求解，由于Vi=0，右侧附加电流也为0
    Va_old = solveMatrix(locs, vals, F, m_num_pts);

    char Vaoldpath[256];
    sprintf(Vaoldpath,"../result/Temp3DTLM_%d_old.txt",m_num_TetEle);
    std::ofstream myoldtemp(Vaoldpath);
    for(int i = 0; i < m_num_pts; i++){
        myoldtemp << mp_3DNode[i].x << " " << mp_3DNode[i].y << " " << mp_3DNode[i].z << " " << Va_old[i] << endl;
    }
    qDebug() << "First solve finish!";

    //8. 后续迭代求解 反射->入射
    const int ITER_MAX = 200;
    for(int iter = 0; iter < ITER_MAX; ++iter){
        vec I = zeros<vec>(m_num_pts);
        for(int k = 0; k < m_num_TetEle; ++k){
            if(mp_TetEle[k].LinearFlag == 0){
                double avgT = (Va_old[mp_TetEle[k].n[0]] + Va_old[mp_TetEle[k].n[1]] + Va_old[mp_TetEle[k].n[2]] + Va_old[mp_TetEle[k].n[3]])/4;
//                double cond = TtoCond(avgT);
                double cond = 0.03;
                for(int i = 0; i < 4; ++i){
                    for(int j = 0; j < 4; ++j){
                        if(i != j){
                            double Si = cond * TetResist[k].C[i][j];
                            double Vd  = Va_old[mp_TetEle[k].n[i]] - Va_old[mp_TetEle[k].n[j]];
                            if(Si < 0){
                                double Vr = Vd - TetVi[k].Vi[i][j];
                                double Vc = 2*Vr*TetY0[k].Y0[i][j]/(TetY0[k].Y0[i][j] - Si);
                                TetVi[k].Vi[i][j] = Vc - Vr;
                                I(mp_TetEle[k].n[i]) += 2*TetY0[k].Y0[i][j]*TetVi[k].Vi[i][j];

//                                double Vr = Va_old[mp_TetEle[k].n[i]] - Va_old[mp_TetEle[k].n[j]] - TetVi[k].Vi[i][j];
//                                double Vc = 2*Vr*Y0/(Y0 - Si);
//                                TetVi[k].Vi[i][j] = Vc - Vr;
//                                I(mp_TetEle[k].n[i]) += 2*Y0*TetVi[k].Vi[i][j];
                            }
                            else{
                                I(mp_TetEle[k].n[i]) = I(mp_TetEle[k].n[i]) + Si*Vd;
                            }
                        }
                    }
                }
            }
        }

        vec F2 = F + I;

        //分解后的LU直接计算
        Va = triangleSolve(F2);
//        Va = solveMatrix(locs, vals, F2, m_num_pts);

        ////////////测试，输出I
        char fpath[256];
        sprintf(fpath,"../result/I_%d.txt",m_num_TetEle);
        std::ofstream myI(fpath);
        for(int i = 0; i < m_num_pts; i++){
            myI << I[i] << endl;
        }
        sprintf(fpath,"../result/Temp3DTLM_%d.txt",m_num_TetEle);;
        //    double temp[15076];
        std::ofstream mytemp(fpath);
        for(int i = 0; i < m_num_pts; i++){
            mytemp << mp_3DNode[i].x << " " << mp_3DNode[i].y << " " << mp_3DNode[i].z << " " << Va[i] << endl;
        }


        //判断收敛性
        double inner_error = 1;
        double a0 = 0, b = 0;
        for(int i = 0; i < m_num_pts; ++i){
            a0 += (Va_old[i] - Va[i])*(Va_old[i] - Va[i]);
            b += Va[i] * Va[i];
        }
        inner_error = sqrt(a0)/sqrt(b);
        qDebug() << "inner_error = " << inner_error;
        if(inner_error > Precision){
            qDebug() << "iter step " << iter;
            for(int i = 0; i < m_num_pts; ++i){
                Va_old[i] = Va[i];
            }
        }else{
            qDebug() << "iter step " << iter;
            for(int i = 0; i < m_num_pts; ++i){
                mp_3DNode[i].V = Va[i];
            }
            break;
        }
    }


//    //输出结果
//    char fpath[256];
//    sprintf(fpath,"../result/Temp3DTLM_%d.txt",m_num_TetEle);
//    std::ofstream mytemp(fpath);
//    //    double temp[15076];
//    for(int i = 0; i < m_num_pts; i++){
//        mytemp << mp_3DNode[i].x << " " << mp_3DNode[i].y << " " << mp_3DNode[i].z << " " << Va[i] << endl;
//    }

    qDebug()<<"Ok.";
    //释放内存
    free(TetVi);
    free(TetY0);
    free(TetResist);

}
