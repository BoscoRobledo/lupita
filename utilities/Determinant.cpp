#include <cmath>
#include <cstdlib>
#include <cstdio>

/*
 * Computes determinant from given matrix using Cholesky decomposition
 *  A: matrix
 *  n: Size of matrix.
 *  TODO: Check if log(det) works as well, in order to get smaller data
 *  returns: determinant.
 */
double choleskyDeterminant(double **A, int n)
{
    double** L=(double**)malloc(n*sizeof(double*));
    for(int c=0;c<n;c++)
        L[c]=(double*)malloc(n*sizeof(double));
    // Cholesky Decomposition from matrix A in L
    for (int i=0;i<n;i++)
        for (int j=0;j<=i;j++)
            {
                double s=0;
                for (int k=0;k<j;k++)
                    s += L[i][k] * L[j][k];
                L[i][j]=(i==j)?sqrt(A[i][i]-s):(1.0/L[j][j]*(A[i][j]-s));
            }
    // Compute determinant
    double det=1;
    for(int c=0;c<n;c++)
        det*=L[c][c]*L[c][c];
    for(int c=0;c<n;c++)
        free(L[c]);
    free(L);
    return det;
}

/*
 * Computes determinant from given matrix using the best strategy for smaller matrices
 *  D: matrix
 *  n: Size of matrix.
 *  TODO: Check if log(det) works as well, in order to get smaller data
 *  returns: determinant.
 */
double determinant(double** D, int order)
{
    int i,j;
    double sum = 0.0;
    if(order==1)
        return D[0][0];
    if(order==2)
    {
        sum = D[0][0]*D[1][1]-D[0][1]*D[1][0];
        return sum;
    }
    if(order>2 && order<6)
    {
        double** aux=new double*[order];
        for(int c=0;c<order;c++)
            aux[c]=new double[order];
        for(int p=0;p<order;p++)
        {
            int h=0, k=0;
            for(i=1;i<order;i++)
            {
                for(j=0;j<order;j++)
                {
                    if(j==p)
                        continue;
                    aux[h][k] = D[i][j];
                    k++;
                    if(k == order-1)
                    {
                        h++;
                        k = 0;
                    }
                }
            }
            sum = sum + D[0][p]*pow(-1,p)*determinant(aux,order-1);
        }
        for(int c=0;c<order;c++)
            delete[] aux[c];
        delete[] aux;
        return sum;
    }
    else
    {
        return choleskyDeterminant(D,order);
    }
}
