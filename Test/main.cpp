#include <iostream>
#include <iomanip>
#include "../UndirectedModels/ChowLiuTree/ChowLiuTree.h"
#include "../UndirectedModels/tCherryJunctionTree/tCherryJunctionTree.h"
using namespace std;

int main()
{
    int D; //dimension of the problem
    int order;   //order of graph structure used in the tree probability approx (1. All variables are independent, 2. Chow & Liu Tree, 3. 3-t-CherryJunctionTree )

    cin>>D>>order;
    double *mu=new double[D];
    double **covM=new double*[D];               /* Population Covariance Matrix */
    double **corrM=new double*[D];               /* Population Correlation Matrix */
    for(int c=0;c<D;c++)
    {
        covM[c]=new double[D];
        corrM[c]=new double[D];
    }

    for(int c=0;c<D;c++)
        cin>>mu[c];

    for(int c=0;c<D;c++)
        for(int d=0;d<D;d++)
            cin>>covM[c][d];

    for(int c=0;c<D;c++)
        for(int d=0;d<D;d++)
            cin>>corrM[c][d];

    ChowLiuTree chowLiuTree(corrM,covM,D);
    chowLiuTree.build();

    cout<<"Construccion de Chow Liu tree"<<endl<<chowLiuTree;

    tCherryJunctionTree tchjtbase(&chowLiuTree,corrM,covM,D);
    tchjtbase.build(2);
    cout<<"Construccion de tCherry Junction Tree k=2 desde chow&liu tree"<<endl<<tchjtbase;

    tCherryJunctionTree tchjt(corrM,covM,D);
    tchjt.build(2);
    cout<<"Construccion de tCherry Junction Tree k=2"<<endl<<tchjt;
    tchjt.increaseOrder();
    cout<<"Actualizacion de orden k=3"<<endl<<tchjt;


    return 0;
}
