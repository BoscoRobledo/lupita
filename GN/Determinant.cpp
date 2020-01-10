#include <cmath>
#include <cstdlib>
#include <cstdio>

double choleskyDeterminant(double **A, int n)
{
    double** L=(double**)malloc(n*sizeof(double*));
    for(int c=0;c<n;c++)
        L[c]=(double*)malloc(n*sizeof(double));
    for (int i=0;i<n;i++)
        for (int j=0;j<=i;j++)
            {
                double s=0;
                for (int k=0;k<j;k++)
                    s += L[i][k] * L[j][k];
                L[i][j]=(i==j)?sqrt(A[i][i]-s):(1.0/L[j][j]*(A[i][j]-s));
            }
    double det=1;
    for(int c=0;c<n;c++)
        det*=L[c][c]*L[c][c];
    for(int c=0;c<n;c++)
        free(L[c]);
    free(L);
    if(det>1e10)
        printf("Esto se esta haciendo muy grande %lf.",det);
    return det;
}


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
