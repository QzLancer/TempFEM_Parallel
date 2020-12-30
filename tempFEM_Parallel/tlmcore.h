#ifndef TLMCORE_H
#define TLMCORE_H

#include "temp3dfemcore.h"

class TLMCore : public Temp3dfemcore
{
public:
    TLMCore();
    TLMCore(const char *fn);
    ~TLMCore();
    void TLMSolve1();   //TLM解决方案1，选取热导率初值带入
    void AdaptiveTLMSolve();

};

#endif // TLMCORE_H
