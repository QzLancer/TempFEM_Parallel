#include "mainwindow.h"
#include "temp3dfemcore.h"

#include <QApplication>
#include <QDebug>

enum Demo{
    DDTLM3D
};

void NR3d(){
    Temp3dfemcore *temp = new Temp3dfemcore("../new_tempFEM/model/mesh_contactor3D_132441.mphtxt");
    temp->load3DFEMCOMSOL();
    temp->preCalculation();
    temp->setCondition();
    temp->NRSolve();
}

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    Demo showwhat = DDTLM3D;
    switch (showwhat) {
    case DDTLM3D:
        NR3d();
        break;
    }
    qDebug() << "Success!!";
    return a.exec();
}

