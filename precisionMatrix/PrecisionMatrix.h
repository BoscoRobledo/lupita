#ifndef PRECMATRIX_HPP_INCLUDED
#define PRECMATRIX_HPP_INCLUDED

void recursiveW(double** W,double** covM, const int n, int* perm);
void getW(double** W, double** covM, double** corrM, double* perm, int D, int k);
void WChL(double** W,double** covM,double** corrM, double* perm, int D);

#endif // PRECMATRIX_HPP_INCLUDED
