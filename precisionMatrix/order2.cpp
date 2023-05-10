#include "Graph.h"
#include "Edge.h"
#include "InEdge.h"
#include "TCHLib.h"
#include <iomanip>
void WChL(double** W,double** covM,double** corrM, double* perm, int D, int k)
{
    /*for(int c=0;c<D;c++)
    {
        for(int d=0;d<D;d++)
            cout<<setprecision(3)<<setw(7)<<corrM[c][d]<<" ";
        cout<<endl;
    }*/

    Graph<int,int> G,CHLT;
    for(int c=0;c<D;c++)
        G.AddVertex(c);
    for(int c=0;c<D;c++)
        for(int d=0;d<D;d++)
            if(c!=d)
                G.AddEdge(c,d,d,corrM[c][d]*corrM[c][d],false);  //build the complete graph from a full correlation matrix
    //cout<<G;
    Edge mx=CHLT.fillMSTFromGraph(G); //Get MST from the graph previously constructed. Undirected graph
    int root=0;
    if(covM[mx.GetStartID()][mx.GetStartID()]<=covM[mx.GetDestID()][mx.GetDestID()])  //Pick a root and turn the graph to a DAG
        root=mx.GetStartID();
    else
        root=mx.GetDestID();
    Graph<int,int> DAGC(true);
    DAGC.fillDAG(root,CHLT);
    //cout<<endl<<"DAG Chow&Liu"<<endl;
    //cout<<DAGC<<endl;

    vector<int> permut=DAGC.DepthFirstSearch(root); //The DFS in a tree gives a permutation of nodes, which turns out to be the tree ancestral ordering.
    //for(int c=0;c<D;c++)
     //   cout<<permut[c]<<" ";


    vector<double> beta(D,0.0);
    vector<double> cv(D,0.0);
    for(int pc=0;pc<D;pc++)                             //Computing conditional variances, and regression coefs, for order 2 only
    {
        int c=permut[pc];
        vector<InEdge<int>> ep=DAGC.GetVertex(c).GetIngoingEdges();
        if(ep.size()==0)
        {
            cv[pc]=covM[c][c];
        }
        else
        {
            int p=ep[0].GetOriginID();
            beta[pc]=covM[p][c]/covM[p][p];
            cv[pc]=covM[c][c]-((covM[p][c]*covM[p][c])/covM[p][p]);
        }
    }
    vector<int> pos(D);
    for(int c=0;c<D;c++)
        pos[permut[c]]=c;
    for(int pc=0;pc<D;pc++)
    {
        int c=permut[pc];
        vector<InEdge<int>> ep=DAGC.GetVertex(c).GetIngoingEdges();
        W[pc][pc]=1.0/cv[pc];
        if(ep.size()>0)
        {
            int pp=pos[ep[0].GetOriginID()];
            W[pc][pp]=W[pp][pc]=-beta[pc]/cv[pc];
            W[pp][pp]=W[pp][pp]+(beta[pc]*beta[pc])/cv[pc];
        }
    }
    for(int c=0;c<D;c++)
        perm[c]=(double)permut[c];

    /*cout<<"W CHLT"<<endl;
    for(int c=0;c<D;c++)
        cout<<permut[c]<<"\t";
    cout<<endl<<endl;
    for(int c=0;c<D;c++)
    {
        for(int d=0;d<D;d++)
            cout<<setprecision(8)<<setw(12)<<W[c][d]<<"\t";
        cout<<endl;
    }*/


    //TCHLib tl(CHLT,0,cM,cV);
    //tl.increaseOrder();
    //cout<<tl;
}


    /*cout<<endl<<"Beta:    ";
    for(int c=0;c<D;c++)
        cout<<setprecision(3)<<setw(7)<<beta[c]<<" ";
    cout<<endl<<"CondCov: ";
    for(int c=0;c<D;c++)
        cout<<setprecision(3)<<setw(7)<<cv[c]<<" ";*/

