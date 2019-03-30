#include "TCHULib.h"
#include "PrecMatrix.h"
#include "Graph.h"
#include "Edge.h"
#include "InEdge.h"
#include <iomanip>

void vectorSqMatrixProduct(double* x,double** A,int n, double* b)
{
    for(int c=0;c<n;c++)
    {
        b[c]=0;
        for(int d=0;d<n;d++)
            b[c]+=x[d]*A[d][c];
    }
}

double vectorDotProduct(double* x, double* y, int n)
{
    double res=0.0;
    for(int c=0;c<n;c++)
        res+=x[c]*y[c];
    return res;
}



void recursiveW(double** W,double** covM, const int n, int* perm)
{
    if(n==1)
        W[0][0]=1.0/covM[perm[n-1]][perm[n-1]];
    else
    {
        recursiveW(W,covM,n-1,perm);
        double* beta=new double[n-1];
        double* sigmaxY=new double[n-1];
        for(int c=0;c<(n-1);c++)
            sigmaxY[c]=covM[perm[c]][perm[n-1]];
        vectorSqMatrixProduct(sigmaxY,W,n-1,beta);
        double condV=covM[perm[n-1]][perm[n-1]]-vectorDotProduct(beta,sigmaxY,n-1);
        for(int c=0;c<(n-1);c++)
        {
            for(int d=0;d<(n-1);d++)
                W[c][d]+=(beta[c]*beta[d])/condV;
            W[c][n-1]=W[n-1][c]=-beta[c]/condV;
        }
        W[n-1][n-1]=1.0/condV;
        delete []sigmaxY;
        delete []beta;
    }
}

void DFSaccumW(TCHULib& chT, double** W, int &windx, bool* visited, int vertex_id)
{
    visited[vertex_id]=true;
    for (OutEdge<Separator>& e : chT.getTCHJT()->GetVertex(vertex_id).GetOutgoingEdges())
    {
        int neighbor_id = e.GetDestID();
        if (!visited[neighbor_id])
        {
            int currV=e.GetData().v[0];
            vector<InEdge<int>> parentset=chT.getDAG()->GetVertex(currV).GetIngoingEdges();
            int pssize=(int)parentset.size();
            vector<int> AO=chT.getPerm();
            vector<int> AOPS;
            vector<int> mp;
            int i=0;
            for(int d : AO)       ///Non-efficient. Sort ingoing edges according to Ancestral ordering. Solution. Ensure ingoing edge is sorted when creating it.
            {
                for(InEdge<int>& e : parentset)
                    if(d==e.GetOriginID())
                    {
                        AOPS.push_back(d);
                        mp.push_back(i);
                        break;
                    }
                i++;
                if((int)AOPS.size()==(chT.getOrder()-1))
                    break;
            }
            double **psW =new double*[pssize];
            double* covPS=new double[pssize];
            double* beta=new double[pssize];
            for(int c=0;c<pssize;c++)
            {
                psW[c]=new double[pssize];
                covPS[c]=chT.getcovM()[AOPS[c]][currV];
            }
            recursiveW(psW,chT.getcovM(),pssize,AOPS.data()); ///Get parent set precision matrix
            vectorSqMatrixProduct(covPS,psW,pssize,beta);
            double condV=chT.getcovM()[currV][currV]-vectorDotProduct(beta,covPS,pssize);
            for(int c=0;c<pssize;c++)
            {
                for(int d=0;d<pssize;d++)
                    W[mp[c]][mp[d]]+=(beta[c]*beta[d])/condV;
                W[windx][mp[c]]=W[mp[c]][windx]=-beta[c]/condV;
                W[windx][windx]=1.0/condV;
            }
            /*printf("Precision Matrix %d\n", windx);
            for(int c=0;c<chT.getN();c++)
            {
                for(int d=0;d<chT.getN();d++)
                    printf("%.8lf ",W[c][d]);
                printf("\n");
            }
            cout<<endl<<endl;
            for(int c=0;c<chT.getN();c++)
            {
                for(int d=0;d<chT.getN();d++)
                    cout<<setprecision(8)<<setw(12)<<W[c][d]<<"\t";
                cout<<endl;
            }*/
            windx++;
            delete[] covPS;
            delete[] beta;
            for(int c=0;c<pssize;c++)
                delete[] psW[c];
            delete[] psW;

            DFSaccumW(chT,W,windx,visited,neighbor_id);
        }
    }
}

void getWfromTChJT(TCHULib& chT, double ** W)
{
    //const int n= chT.getN();
    const int k= chT.getOrder();
    double** covM= chT.getcovM();
    vector<int> perm=chT.getPerm();
    ///find W of the complete graph defined in the root.
    recursiveW(W,covM,k,perm.data());
    /*printf("Precision Matrix root\n");
    for(int c=0;c<chT.getN();c++)
    {
        for(int d=0;d<chT.getN();d++)
            printf("%.8lf ",W[c][d]);
        printf("\n");
    }*/
    ///Do a DFS, and for each dominant vertex of each node of the cherry tree acumulate W. Parents of each dominating vertex are all connected.
    bool* visited=new bool[chT.getTCHJT()->VertexCount()];
    for(int c=0;c<chT.getTCHJT()->VertexCount();c++)
        visited[c]=false;
    int wIndx=k;
    DFSaccumW(chT,W,wIndx,visited,chT.getRoot());
}




void WChL(double** W,double** covM,double** corrM, double* perm, int D)
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
        /*printf("Precision Matrix %d\n",pc);
        for(int c=0;c<D;c++)
        {
            for(int d=0;d<D;d++)
                printf("%.8lf ",W[c][d]);
            printf("\n");
        }*/
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


void getW(double** W, double** covM,double** corrM, double* perm, int D, int k)
{
    for(int c=0;c<D;c++)
        for(int d=0;d<D;d++)
            W[c][d]=0.0;
    vector<int> p;
    if(k==0)
        cout<<"Can't build a t-Cherry Junction Tree of order 0";
    else if(k==1)                           ///All variables are independent
    {
        for(int c=0;c<D;c++)
            p.push_back(c);
        for(int c=0;c<D;c++)
            W[c][c]=1.0/covM[c][c];
    }
    else if(k==D)                                   ///Complete covariance matrix inverse
    {
        for(int c=0;c<D;c++)
            p.push_back(c);
        recursiveW(W,covM,D,p.data());
    }
    else
    {
        Graph<int,int> G,CHLT;      ///Build first a chow&liu tree
        for(int c=0;c<D;c++)
            G.AddVertex(c);
        for(int c=0;c<D;c++)
            for(int d=0;d<D;d++)
                if(c!=d)
                    G.AddEdge(c,d,d,corrM[c][d]*corrM[c][d],false);  ///build a complete graph from a full correlation matrix
        Edge mx=CHLT.fillMSTFromGraph(G);                       ///Extract MST
        int root=0;
        if(covM[mx.GetStartID()][mx.GetStartID()]<=covM[mx.GetDestID()][mx.GetDestID()])  //Pick a root and turn the graph to a DAG
            root=mx.GetStartID();
        else
            root=mx.GetDestID();

        TCHULib ch3(CHLT,root,corrM,covM,false);     ///Convert Chow&Liu tree into a 2nd order t-cherry juncton tree
        for(int c=3;c<=k;c++)
            ch3.increaseOrder();                                ///increase order till kth order
        ch3.buildDAG();
        getWfromTChJT(ch3,W);
        p=ch3.getPerm();
        //for(int c=0;c<ch3.getN();c++)
        //    cout<<p[c]<<"\t";
        //cout<<endl<<"DAG tchjt 2"<<endl;

        //cout<<*ch3.getDAG();
    }

    for(int c=0;c<D;c++)
        perm[c]=(double)p[c];

}



/*
cout<<"W TCHJT 2"<<endl;
    f
    cout<<endl<<endl;
    for(int c=0;c<ch3.getN();c++)
    {
        for(int d=0;d<ch3.getN();d++)
            cout<<setprecision(8)<<setw(12)<<W[c][d]<<"\t";
        cout<<endl;
    }
*/
