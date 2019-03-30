#ifndef PRECMATRIX_HPP_INCLUDED
#define PRECMATRIX_HPP_INCLUDED

void vectorSqMatrixProduct(double* x,double** A,int n, double* b);
double vectorDotProduct(double* x, double* y, int n);
void recursiveW(double** W,double** covM, const int n, int* perm);
void getW(double** W, double** covM, double** corrM, double* perm, int D, int k);
void WChL(double** W,double** covM,double** corrM, double* perm, int D);
#endif // PRECMATRIX_HPP_INCLUDED
