#include <iostream>
#include <iomanip>
//#include "sampling/Sampler.h"
#include "GN/PrecMatrix.h"
using namespace std;
typedef void (*Wfun)(double**,double**,double**,double*,int, int);
int main()
{
    int D; //dimension of the problem
    int order;   //order of graph structure used in the tree probability approx (1. All variables are independent, 2. Chow & Liu Tree, 3. 3-t-CherryJunctionTree )

    cin>>D>>order;
    double *mu=new double[D];
    double **covM=new double*[D];               /* Population Covariance Matrix */
    double **corrM=new double*[D];               /* Population Correlation Matrix */
    double *perm=new double[D];
    double **W=new double*[D];
    for(int c=0;c<D;c++)
    {
        W[c]=new double[D];
        covM[c]=new double[D];
        corrM[c]=new double[D];
        for(int d=0;d<D;d++)
            W[c][d]=0.0;
    }
    for(int c=0;c<D;c++)
        cin>>mu[c];

    for(int c=0;c<D;c++)
        for(int d=0;d<D;d++)
            cin>>covM[c][d];

    for(int c=0;c<D;c++)
        for(int d=0;d<D;d++)
            cin>>corrM[c][d];
    for(int c=0;c<D;c++)
    {
        for(int d=0;d<D;d++)
            cout<<setprecision(12)<<setw(15)<<covM[c][d]<<" ";
        cout<<endl;
    }
    cout<<"Correlación"<<endl;
    for(int c=0;c<D;c++)
    {
        for(int d=0;d<D;d++)
            cout<<setprecision(12)<<setw(15)<<corrM[c][d]<<" ";
        cout<<endl;
    }


    getW(W,covM,corrM,perm,D,2);

    for(int c=0;c<D;c++)
        cout<<perm[c]<<" ";
    cout<<endl;
    for(int c=0;c<D;c++)
    {
        for(int d=0;d<D;d++)
            cout<<setprecision(12)<<setw(15)<<W[c][d]<<" ";
        cout<<endl;
    }
    return 0;
}




/*
    int nsamples=10000;
    double **samples=new double*[nsamples];
    for(int c=0;c<nsamples;c++)
        samples[c]=new double[D];
    getSamplesTChJT(mu,covM,corrM,samples,D,nsamples,order);
    for(int c=0;c<nsamples;c++)
    {
        for(int d=0;d<D;d++)
            cout<<setprecision(12)<<setw(15)<<samples[c][d]<<" ";
        cout<<endl;
    }*/

    /*Graph<int> CHL(false);
    for(int c=0;c<=6;c++)
        CHL.AddVertex(c);
    CHL.AddEdge(0,1,23.3,false);
    CHL.AddEdge(1,2,12.3,false);
    CHL.AddEdge(2,3,22.3,false);
    CHL.AddEdge(3,4,45.3,false);
    CHL.AddEdge(3,5,34.3,false);
    CHL.AddEdge(0,1,23.3,false);
    CHL.AddEdge(0,2,12.3,false);
    CHL.AddEdge(0,3,22.3,false);
    CHL.AddEdge(3,4,45.3,false);
    CHL.AddEdge(2,5,34.3,false);
    CHL.AddEdge(2,6,32.3,false);
for(int c=0;c<D;c++)
    cout<<setprecision(3)<<setw(7)<<mu[c]<<" ";
    cout<<endl;

    for(int c=0;c<D;c++)
    {
        for(int d=0;d<D;d++)
            cout<<setprecision(3)<<setw(7)<<covM[c][d]<<" ";
        cout<<endl;
    }

    for(int c=0;c<D;c++)
    {

        for(int d=0;d<D;d++)
            cout<<setprecision(12)<<setw(15)<<corrM[c][d]<<" ";
        cout<<endl;
    }

    cout<<CHL;*/




