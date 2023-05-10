#include "order1.h"
void DW(double** W,double** covM,double** corrM,double* perm, int D, int k)
{
    for(int c=0;c<D;c++)
    {
        W[c][c]=1.0/covM[c][c];
        perm[c]=c;
    }
}
