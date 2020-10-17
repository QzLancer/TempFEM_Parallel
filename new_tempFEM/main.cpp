#include "mainwindow.h"
#include "temp3dfemcore.h"

#include <QApplication>
#include <QDebug>

enum Demo{
    NR3D,
    DDTLM3D
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

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    Demo showwhat = DDTLM3D;
    switch (showwhat) {
    case NR3D:
        NR3d();
        break;
    case DDTLM3D:
        DDTLM3d(8);
        break;
    }
    qDebug() << "Success!!";
    return a.exec();
}

