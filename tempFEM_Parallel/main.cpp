#include <iostream>

#include "tlmcore.h"

using namespace std;

void TLM3d1(){
    TLMCore *temp = new TLMCore("../new_tempFEM/model/mesh_contactor3D_455944.mphtxt");
    temp->load3DFEMCOMSOL();
    temp->preCalculation();
    temp->setCondition();
    temp->TLMSolve1();
}

void NR3d(){
    Temp3dfemcore *temp = new Temp3dfemcore("../new_tempFEM/model/mesh_contactor3D_212662.mphtxt");
    temp->load3DFEMCOMSOL();
    temp->preCalculation();
    temp->setCondition();
    temp->NRSolve();
}

int main()
{
    TLM3d1();
    return 0;
}
