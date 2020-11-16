//#include "mainwindow.h"
#include "temp3dfemcore.h"
#include "tlmcore.h"

#include <QApplication>
#include <QDebug>

enum Demo{
    NR3D,
    DDTLM3D,
    TLM1
};

void NR3d(){
    Temp3dfemcore *temp = new Temp3dfemcore("../new_tempFEM/model/mesh_contactor3D_212662.mphtxt");
    temp->load3DFEMCOMSOL();
    temp->preCalculation();
    temp->setCondition();
    temp->NRSolve();
}

void DDTLM3d(int part){
    Temp3dfemcore *temp = new Temp3dfemcore("../new_tempFEM/model/mesh_contactor3D_212662.mphtxt");
    temp->load3DFEMCOMSOL();
    temp->preCalculation();
    temp->setCondition();
    temp->GenerateMetisMesh(part);
    temp->DDTLMSolve();
}

void TLM3d1(){
    TLMCore *temp = new TLMCore("../new_tempFEM/model/mesh_contactor3D_212662.mphtxt");
    temp->load3DFEMCOMSOL();
    temp->preCalculation();
    temp->setCondition();
    temp->TLMSolve1();
}

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    Demo showwhat = TLM1;
    switch (showwhat) {
    case NR3D:
        NR3d();
        break;
    case DDTLM3D:
        DDTLM3d(16);
        break;
    case TLM1:
        TLM3d1();
        break;
    }
    qDebug() << "Success!!";
    return a.exec();
}

