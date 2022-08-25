#include <iostream>

#include "tlmcore.h"
#include "slu_mt_ddefs.h"

using namespace std;

Modeltype model = MODEL2;

void TLM3d(){
    TLMCore *temp = new TLMCore("../tempFEM_Parallel/model/mesh_relay3D_118129.mphtxt");
    temp->load3DFEMCOMSOL();
    temp->preCalculation();
    temp->setCondition(model);
    temp->TLMSolve1();
}

void NR3d(){
    Temp3dfemcore *temp = new Temp3dfemcore("../tempFEM_Parallel/model/mesh_relay3D_118129.mphtxt");
    temp->load3DFEMCOMSOL();
    temp->preCalculation();
    temp->setCondition(model);
    temp->NRSolve();
}

void AdaptiveTLM3d(){
    TLMCore *temp = new TLMCore("../tempFEM_Parallel/model/mesh_relay3D_118129.mphtxt");
    temp->load3DFEMCOMSOL();
    temp->preCalculation();
    temp->setCondition(model);
    temp->AdaptiveTLMSolve();
}

int main()
{
    double t1 = SuperLU_timer_();
    NR3d();
    double t2 = SuperLU_timer_() - t1;
    cout << endl << "Total time: " << t2 << endl;
    return 0;
}
