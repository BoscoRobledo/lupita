#include "UndirectedModel.h"
const int UndirectedModel::getOrder() const
{
    return k;
}
const int UndirectedModel::getDimensionality() const
{
    return d;
}
double** UndirectedModel::getcovM()
{
    return covM;
}
double** UndirectedModel::getcorrM()
{
    return corrM;
}

UndirectedModel::UndirectedModel(double** _corrM, double** _covM, int _d, int _k)
{
    corrM=_corrM;
    covM=_covM;
    d=_d;
    k=_k;
}

ostream & operator<<(ostream & out, UndirectedModel & g)
{
    out<<g.structure;
    return out;
}

const int UndirectedModel::getBestVertex() const
{
    return bestV;
}
