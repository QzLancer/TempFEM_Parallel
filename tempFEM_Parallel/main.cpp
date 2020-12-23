#include <iostream>

#include "tlmcore.h"
#include "slu_mt_ddefs.h"

using namespace std;

Modeltype model = MODEL2;

void TLM3d1(){
    TLMCore *temp = new TLMCore("../tempFEM_Parallel/model/mesh_relay3D_713573.mphtxt");
    temp->load3DFEMCOMSOL();
    temp->preCalculation();
    temp->setCondition(model);
    temp->TLMSolve1();
}

void NR3d(){
    Temp3dfemcore *temp = new Temp3dfemcore("../tempFEM_Parallel/model/mesh_relay3D_713573.mphtxt");
    temp->load3DFEMCOMSOL();
    temp->preCalculation();
    temp->setCondition(model);
    temp->NRSolve();
}

int main()
{
    double t1 = SuperLU_timer_();
    TLM3d1();
    double t2 = SuperLU_timer_() - t1;
    cout << endl << "Total time: " << t2 << endl;
    return 0;
}
