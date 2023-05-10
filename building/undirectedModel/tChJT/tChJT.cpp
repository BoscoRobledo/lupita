#include "../../../graph/Graph.h"
#include "tChJT.h"
#include "DeterminantCalculator.h"
#include <queue>
#include <iostream>
#include <fstream>
#include <dirent.h>

///-----Properties
double tChJT::getWeight()
{
    return weight;
}

double tChJT::getModelDifferentialEntropy()
{
    return diffEntropy;
}

///-----Properties

///ctor - dtor

void tChJT::initVars()
{
    weight=0.0;
    diffEntropy=0.0;
    k=0;
}


tChJT::tChJT(ChLT* chowLiuTree, double** corrM, double** covM, int d) : UndirectedModel(corrM,covM,d,0), baseModel(chowLiuTree)
{
    initVars();
}

tChJT::tChJT(double** corrM, double** covM, int d) : UndirectedModel(corrM,covM,d,0)
{
    initVars();
    baseModel=new ChLT(corrM,covM,d);
    baseModel->build();
}

tChJT::~tChJT()
{
}

///ctor - dtor


///-----Construction of t=2 Cherry Junction Tree from a Chow&Liu Tree

void tChJT::Getdonating_V_ariable(int CD, vector<int> &s, int CA)
{
    for (int vD : tchjt.GetVertexData(CD).nodes)
    {
        bool exists=false;
        for (int vS : tchjt.GetVertexData(CA).nodes)
            if(vS==vD)
            {
                exists=true;
                break;
            }
        if(!exists)
            s.push_back(vD);
    }
}

void tChJT::DFSChLT(int vertex_id, vector<bool> & visited, bool root, int parent_cluster)
{
    visited[vertex_id] = true;
    for (Edge<double> e : baseModel->structure.GetVertex(vertex_id).GetOutgoingEdges())
    {
        int neighbor_id = e.GetDestinationID();
        if (!visited[neighbor_id])
        {
            Cluster cl;
            cl.nodes.push_back(vertex_id);
            cl.nodes.push_back(neighbor_id);
            cl.deleted=false;
            cl.updated=false;
            cl.weight=-0.5*log(1.0-corrM[vertex_id][neighbor_id]*corrM[vertex_id][neighbor_id]);
            cl.diffEntropy=0.5*log(2.0*M_PI*M_E*2.0*M_PI*M_E*(covM[vertex_id][vertex_id]*covM[neighbor_id][neighbor_id]-covM[vertex_id][neighbor_id]*covM[neighbor_id][vertex_id]));
            weight+=cl.weight;
            diffEntropy+=cl.diffEntropy;
            int cluster_id=tchjt.AddVertex(cl);
            if(!root)
            {
                //create edges in both ways
                Separator s,s2;
                s.nodes.push_back(vertex_id);
                s2.nodes.push_back(vertex_id);
                //select proper Xv in forward separator
                if(vertex_id==tchjt.GetVertexData(parent_cluster).nodes[0])
                    s.Xv.push_back(tchjt.GetVertexData(parent_cluster).nodes[1]);
                else
                    s.Xv.push_back(tchjt.GetVertexData(parent_cluster).nodes[0]);
                s2.Xv.push_back(neighbor_id);
                s.Xu=s2.Xv[0];
                s2.Xu=s.Xv[0];

                //initialization of weight-determinant calculation variables
                double* mrow=new double[1];
                s.corrMCholDec.push_back(mrow);
                s2.corrMCholDec.push_back(mrow);
                mrow=new double[1];
                s.covMCholDec.push_back(mrow);
                s2.covMCholDec.push_back(mrow);
                s.corrMCholDec[0][0]=sqrt(corrM[vertex_id][vertex_id]);
                s.covMCholDec[0][0]=sqrt(covM[vertex_id][vertex_id]);
                s.corrMDet=corrM[vertex_id][vertex_id];
                s.covMDet=covM[vertex_id][vertex_id];
                s.diffEntropy=s2.diffEntropy=0.0;
                s.complete=s2.complete=true;
                s.weight=s2.weight=0.0;


                tchjt.AddEdge(parent_cluster,cluster_id,s);
                tchjt.AddEdge(cluster_id,parent_cluster,s2);
            }
            else
            {
                parent_cluster=cluster_id;
                bestV=cluster_id;
            }
            root=false;
            DFSChLT(neighbor_id,visited,false,cluster_id);
        }
    }
}

void tChJT::build(int _k)
{
    vector<bool> visited(d,false);
    DFSChLT(baseModel->getBestVertex(),visited,true,-1);
    k=baseModel->getOrder();
    for(int kappa=3;kappa<=k;kappa++)
        increaseOrder();

}


///-----Construction from 2 tCherry Junction Tree (Chow&Liu Tree)


///-----Order Update

//Equation 10 implementation
void tChJT::setW(PotentialUpdate &gamma,Edge<Separator>& e)
{
    vector<int> nodes;
    //ToDo: This can be done in O(1)
    gamma.corrMCholDec.assign(e.GetData().corrMCholDec.begin(),e.GetData().corrMCholDec.end());
    gamma.covMCholDec.assign(e.GetData().covMCholDec.begin(),e.GetData().covMCholDec.end());
    nodes.assign(e.GetData().nodes.begin(),e.GetData().nodes.end());

    //this section makes this function O(k^2) complexity
        //Lambda_S_ij U X_u
    gamma.corrMCholDec.push_back(new double[k]);
    gamma.covMCholDec.push_back(new double[k]);
    int nodes_size=k-1;
    nodes.push_back(gamma.Xu);
    for(int j=0;j<=nodes_size;j++)
    {
        double sumCorr=0.0, sumCov=0.0;
        for (int i=0;i<j;i++)
        {
            sumCorr += gamma.corrMCholDec[nodes_size][i] * gamma.corrMCholDec[j][i];
            sumCov += gamma.covMCholDec[nodes_size][i] * gamma.covMCholDec[j][i];
        }
        gamma.corrMCholDec[nodes_size][j]=(nodes_size==j)?sqrt(corrM[nodes[nodes_size]][nodes[nodes_size]]-sumCorr):((1.0/gamma.corrMCholDec[j][j])*(corrM[nodes[nodes_size]][nodes[j]]-sumCorr));
        gamma.covMCholDec[nodes_size][j]=(nodes_size==j)?sqrt(covM[nodes[nodes_size]][nodes[nodes_size]]-sumCov):((1.0/gamma.covMCholDec[j][j])*(covM[nodes[nodes_size]][nodes[j]]-sumCov));
    }
    double num=gamma.corrMCholDec[nodes_size][nodes_size]*gamma.corrMCholDec[nodes_size][nodes_size];
    //Lambda_S_ij U X_v U X_u
    nodes[nodes_size]=gamma.Xv;
    for(int j=0;j<=nodes_size;j++)
    {
        double sumCorr=0.0, sumCov=0.0;
        for (int i=0;i<j;i++)
        {
            sumCorr += gamma.corrMCholDec[nodes_size][i] * gamma.corrMCholDec[j][i];
            sumCov += gamma.covMCholDec[nodes_size][i] * gamma.covMCholDec[j][i];
        }
        gamma.corrMCholDec[nodes_size][j]=(nodes_size==j)?sqrt(corrM[nodes[nodes_size]][nodes[nodes_size]]-sumCorr):((1.0/gamma.corrMCholDec[j][j])*(corrM[nodes[nodes_size]][nodes[j]]-sumCorr));
        gamma.covMCholDec[nodes_size][j]=(nodes_size==j)?sqrt(covM[nodes[nodes_size]][nodes[nodes_size]]-sumCov):((1.0/gamma.covMCholDec[j][j])*(covM[nodes[nodes_size]][nodes[j]]-sumCov));
    }
    gamma.corrMCholDec.push_back(new double[k+1]);
    gamma.covMCholDec.push_back(new double[k+1]);
    nodes_size=k;
    nodes.push_back(gamma.Xu);
    for(int j=0;j<=nodes_size;j++)
    {
        double sumCorr=0.0, sumCov=0.0;
        for (int i=0;i<j;i++)
        {
            sumCorr += gamma.corrMCholDec[nodes_size][i] * gamma.corrMCholDec[j][i];
            sumCov += gamma.covMCholDec[nodes_size][i] * gamma.covMCholDec[j][i];
        }
        gamma.corrMCholDec[nodes_size][j]=(nodes_size==j)?sqrt(corrM[nodes[nodes_size]][nodes[nodes_size]]-sumCorr):(1.0/gamma.corrMCholDec[j][j]*(corrM[nodes[nodes_size]][nodes[j]]-sumCorr));
        gamma.covMCholDec[nodes_size][j]=(nodes_size==j)?sqrt(covM[nodes[nodes_size]][nodes[nodes_size]]-sumCov):(1.0/gamma.covMCholDec[j][j]*(covM[nodes[nodes_size]][nodes[j]]-sumCov));
    }
    double den=gamma.corrMCholDec[nodes_size][nodes_size]*gamma.corrMCholDec[nodes_size][nodes_size];
    gamma.wInc=0.5*log(num/den);

}

    ///--Priority queue creation
void tChJT::FillPriorityQueue()
{
    vector<bool> visited(tchjt.VertexCount(),false);
    DFSPQ(getBestVertex(),visited);
}


void tChJT::DFSPQ(int vertex_id, vector<bool> & visited)
{
    visited[vertex_id] = true;
    for (Edge<Separator> e : tchjt.GetVertex(vertex_id).GetOutgoingEdges())
    {
        int neighbor_id = e.GetDestinationID();
        if (!visited[neighbor_id])
        {
            PotentialUpdate gamma,gamma2;
            gamma.CA=gamma2.CD=vertex_id;
            gamma.CD=gamma2.CA=neighbor_id;
            gamma.Xv=gamma2.Xu=e.GetData().Xv[0];
            gamma2.Xv=gamma.Xu=e.GetData().Xu;
            gamma.updateCase=gamma2.updateCase=4;
            setW(gamma,e);
            //ToDo. this can be replicated to some extent from previous calculus
            setW(gamma2,e);
            Gamma.push(gamma);
            Gamma.push(gamma2);
            DFSPQ(neighbor_id,visited);
        }
    }
}
///--Priority queue creation

bool tChJT::isValid(PotentialUpdate& T)
{
    //If donor cluster does not exist or active cluster has already k+1 vars, continue
    if(tchjt.GetVertex(T.CA).GetData().updated || tchjt.GetVertex(T.CD).GetData().deleted)
        return false;
    //If both clusters size are k+1 then no potential update step exists
    if((int)(tchjt.GetVertex(T.CA).GetData().nodes.size())>k && (int)(tchjt.GetVertex(T.CD).GetData().nodes.size())>k)
        return false;
    return true;
}


void tChJT::increaseOrder()
{
    FillPriorityQueue();
    //Priority queue printing
    priority_queue<PotentialUpdate> ax;
    cout<<"Valores de cola de prioridad"<<endl<<"CA-CD-XV-XU-WINC"<<endl;
    while(!Gamma.empty())
    {
        PotentialUpdate t=Gamma.top();
        ax.push(t);
        Gamma.pop();
        cout<<t.CA<<" "<<t.CD<<" "<<t.Xv<<" "<<t.Xu<<" "<<t.wInc<<endl;
    }
    while(!ax.empty())
    {
        PotentialUpdate t=ax.top();
        Gamma.push(t);
        ax.pop();
    }
    //-------------

    int updatableClusters=tchjt.VertexCount();
    double maxClusterWeight = 0.0;
    while(updatableClusters > 0)
    {
        PotentialUpdate T = Gamma.top();
        Gamma.pop();
        cout<<"Potential update\n"<<T.CA<<" "<<T.CD<<" "<<T.Xv<<" "<<T.Xu<<" "<<T.wInc<<": ";
        ///If maximum weight change is valid, add var v to active cluster
        if(isValid(T))
        {
            cout<<"Taken!!!\n";
            tchjt.GetVertex(T.CA).GetData().nodes.push_back(T.Xv);
            tchjt.GetVertex(T.CA).GetData().weight=T.wSXvXu;
            tchjt.GetVertex(T.CA).GetData().updated=true;
            updatableClusters--;
            if(T.updateCase==4)
            {
                for (Edge<Separator>& e : tchjt.GetVertex(T.CD).GetOutgoingEdges())       ///Link donor cluster neighbors to CA'
                    if(!(tchjt.GetVertex(e.GetDestinationID()).GetData().deleted) && e.GetDestinationID()!=T.CA)
                    {
                        tchjt.AddEdge(e.GetDestinationID(),T.CA,e.GetData(),false);  ///In case 4 new separators are always of size k-1
                        tchjt.AddEdge(T.CA,e.GetDestinationID(),e.GetData(),false);   ///So keeping the data in old ones works fine
                        if(!(tchjt.GetVertex(e.GetDestinationID()).GetData().updated))
                        {
                            PotentialUpdate gamma;
                            gamma.Xv=T.Xu;                ///Now variable Xu can can be donated from active cluster to the neighbors of CD
                            gamma.CA=e.GetDestinationID();
                            gamma.CD=T.CA;
                            gamma.Xu=T.Xu;
                            gamma.updateCase=2;
                            setW(gamma, e);  ///Donor cluster (active cluster) is now size k+1, and its active neighbor cluster is size K, so case 2 always happen
                            Gamma.push(gamma);
                        }
                    }
                ///case 4 happens, so "drop" donor cluster
                tchjt.GetVertex(T.CD).GetData().updated=true;
                tchjt.GetVertex(T.CD).GetData().deleted=true;
                updatableClusters--;
            }
            else if(T.updateCase == 2)
            {
                for (Edge<Separator>& e : tchjt.GetVertex(T.CA).GetOutgoingEdges())       ///Find the separators between CA and CD and add the variable v
                    if(e.GetDestinationID()==T.CD)
                    {
                        e.GetData().nodes.push_back(T.Xv);
                        e.GetData().weight=T.wS;
                        e.GetData().corrMCholDec=T.corrMCholDec;
                        e.GetData().covMCholDec=T.covMCholDec;
                    }

                for (Edge<Separator>& e : tchjt.GetVertex(T.CD).GetOutgoingEdges())       ///Find the separator between CA and CD and add the variable v
                    if(e.GetDestinationID()==T.CA)
                    {
                        e.GetData().nodes.push_back(T.Xv);
                        e.GetData().weight=T.wS;
                        e.GetData().corrMCholDec=T.corrMCholDec;
                        e.GetData().covMCholDec=T.covMCholDec;
                    }
            }
            for (Edge<Separator>& e : tchjt.GetVertex(T.CA).GetOutgoingEdges())       ///Add potential steps generated from adding the variable v to the neighbors of active cluster and deleting CD if happened
                if(!(tchjt.GetVertex(e.GetDestinationID()).GetData().deleted) && !(tchjt.GetVertex(e.GetDestinationID()).GetData().updated))  ///If new or old neighbor of CA is still active, there is a potential step to add
                {
                    PotentialUpdate gamma;
                    gamma.Xv=T.Xv;                ///Now variable v can can be donated from active cluster to its active neighbors
                    gamma.CA=e.GetDestinationID();
                    gamma.CD=T.CA;
                    gamma.Xu=T.Xu;
                    gamma.updateCase=2;
                    setW(gamma, e);  ///Donor cluster (active cluster) is now size k+1, and its active neighbor cluster is size K, so case 2 always happen
                    Gamma.push(gamma);
                }
            if(maxClusterWeight<tchjt.GetVertex(T.CA).GetData().weight)
            {
                rootCluster=T.CA;
                maxClusterWeight=tchjt.GetVertex(T.CA).GetData().weight;
            }
        }
        else
        {
            cout<<"Dropped :(\n";
        }
    }
    k++;
    //ProcessBUDS();
    //clean();
}


/*void tChJT::ProcessBUDS()
{
    vector<bool> visited(tcjt->VertexCount(),false);
    vector<pair<int,pair<int,vector<int>>>> buds;
    DFSGetBUDS(root,visited,buds);
    for(pair<int,pair<int,vector<int>>>& pp : buds)
    {
        int vertex_id=pp.first, neighbor_id=pp.second.first;
        vector<int> v1,v2;
        Getdonating_V_ariable(vertex_id,v1,neighbor_id);
        Getdonating_V_ariable(neighbor_id,v2,vertex_id);
        double maxW=-1e40, ws1=0.0,ws2=0.0,wc=0.0;

        int p1,p2;
        for(int c : v1)
            for(int d : v2)
            {
                double wi=0.0,wsi1=0.0,wsi2=0.0,wci=0.0;
                getWBUD(wi,ws1,ws2,wc,pp.second.second,c,d);
                if(wi>maxW)
                {
                    maxW=wi;
                    ws1=wsi1;
                    ws2=wsi2;
                    wc=wci;
                    p1=c;
                    p2=d;
                }
            }
        Cherry vn;
        vn.active=vn.exists=true;
        for(int i : pp.second.second)
            vn.nodes.push_back(i);
        vn.nodes.push_back(p1);
        vn.nodes.push_back(p2);
        vn.weight=wc;
        int v_id=tcjt->AddVertex(vn);
        Separator nS1,nS2,nS3,nS4;
        for(int i : pp.second.second)
        {
            nS1.nodes.push_back(i);
            nS2.nodes.push_back(i);
            nS3.nodes.push_back(i);
            nS4.nodes.push_back(i);             ///Build the four separators
        }
        nS1.nodes.push_back(p1);
        nS2.nodes.push_back(p1);
        nS3.nodes.push_back(p2);
        nS4.nodes.push_back(p2);

        int p01=v1[0], p02=v2[0];
        if(p01==p1)
            p01=v1[1];
        if(p02==p2)
            p02=v2[1];

        nS1.dVA=p01;
        nS1.v.push_back(p2);
        nS1.complete=true;
        nS1.weight=ws1;
        tcjt->AddEdge(vertex_id,v_id,nS1,0.0,false);

        nS2.dVA=p2;
        nS2.v.push_back(p01);
        nS2.complete=true;
        nS2.weight=ws1;
        tcjt->AddEdge(v_id,vertex_id,nS2,0.0,false);

        nS3.dVA=p1;
        nS3.v.push_back(p02);
        nS3.complete=true;
        nS3.weight=ws2;
        tcjt->AddEdge(v_id,neighbor_id,nS3,0.0,false);

        nS4.dVA=p02;
        nS4.v.push_back(p1);
        nS4.complete=true;
        nS4.weight=ws2;
        tcjt->AddEdge(neighbor_id,v_id,nS4,0.0,false);
    }
}



void tChJT::DFSGetBUDS(int vertex_id, vector<bool>& visited, vector<pair<int,pair<int,vector<int>>>>& buds)
{
    visited[vertex_id] = true;
    for (OutEdge<Separator>& e : tcjt->GetVertex(vertex_id).GetOutgoingEdges())
    {
        int neighbor_id = e.GetDestID();
        if (!visited[neighbor_id] && tcjt->GetVertexData(neighbor_id).exists)
        {
            if((int)e.GetData().nodes.size()<(k-1))
                buds.push_back(make_pair(vertex_id,make_pair(neighbor_id,e.GetData().nodes)));
            DFSGetBUDS(neighbor_id,visited,buds);
        }

    }
}*/



/*void tChJT::getWBUD(double &w, double &wS1, double &wS2, double &wC12, vector<int>& sep, int v1, int v2)
{
    double dC,dSuv2, dSuv1, dS;
    int sCA=sep.size()+2, ssep=sep.size();
    double** aux=new double*[sCA];
    for(int c=0;c<sCA;c++)
        aux[c]=new double[sCA];
    //separator
    for(int c=0;c<ssep;c++)
            for(int d=0;d<ssep;d++)
                aux[c][d]=corrM[sep[c]][sep[d]];
    dS = determinant(aux,ssep);
    //separator union v1
    for(int c=0;c<ssep;c++)
        aux[ssep][c]=aux[c][ssep]=corrM[sep[c]][v1];
    aux[ssep][ssep]=corrM[v1][v1];
    dSuv1=determinant(aux,ssep+1);
    wS1=-0.5*log(dSuv1);
    //separator union v1
    for(int c=0;c<ssep;c++)
        aux[ssep][c]=aux[c][ssep]=corrM[sep[c]][v2];
    aux[ssep][ssep]=corrM[v2][v2];
    dSuv2=determinant(aux,ssep+1);
    wS2=-0.5*log(dSuv2);

    //New Cluster
    for(int c=0;c<ssep;c++)
        aux[ssep+1][c]=aux[c][ssep+1]=corrM[sep[c]][v1];
    aux[ssep+1][ssep]=aux[ssep][ssep+1]=corrM[v2][v1];
    aux[ssep+1][ssep+1]=corrM[v2][v2];
    dC=determinant(aux,ssep+2);
    wC12=-0.5*log(dC);
    if((dSuv1*dSuv2)/(dC*dS)<-1e40)
    {
        dC=0.0;
        for(int c=0;c<n;c++)
        {
            for(int d=0;d<n;d++)
                cout<<corrM[c][d]<<" ";
            cout<<endl;
        }
        for(int c=0;c<n;c++)
        {
            for(int d=0;d<n;d++)
                cout<<covM[c][d]<<" ";
            cout<<endl;
        }
        //separator
        for(int c=0;c<ssep;c++)
                for(int d=0;d<ssep;d++)
                    aux[c][d]=corrM[sep[c]][sep[d]];
        cout<<"separator: ";
        for(int c=0;c<ssep;c++)
            cout<<sep[c]<<" ";
        cout<<endl;
        for(int c=0;c<ssep;c++)
        {
            for(int d=0;d<ssep;d++)
                cout<<aux[c][d]<<" ";
            cout<<endl;
        }



        //separator union v1
        for(int c=0;c<ssep;c++)
            aux[ssep][c]=aux[c][ssep]=corrM[sep[c]][v1];
        aux[ssep][ssep]=corrM[v1][v1];
        dSuv1=determinant(aux,ssep+1);
        cout<<"separator union v1:"<<v1<<endl;
        for(int c=0;c<ssep+1;c++)
        {
            for(int d=0;d<ssep+1;d++)
                cout<<aux[c][d]<<" ";
            cout<<endl;
        }
        //separator union v1
        for(int c=0;c<ssep;c++)
            aux[ssep][c]=aux[c][ssep]=corrM[sep[c]][v2];
        aux[ssep][ssep]=corrM[v2][v2];
        dSuv2=determinant(aux,ssep+1);
        cout<<"separator union v2:"<<v2<<endl;
        for(int c=0;c<ssep+1;c++)
        {
            for(int d=0;d<ssep+1;d++)
                cout<<aux[c][d]<<" ";
            cout<<endl;
        }

        //New Cluster
        for(int c=0;c<ssep;c++)
            aux[ssep+1][c]=aux[c][ssep+1]=corrM[sep[c]][v1];
        aux[ssep+1][ssep]=aux[ssep][ssep+1]=corrM[v2][v1];
        aux[ssep+1][ssep+1]=corrM[v2][v2];
        dC=determinant(aux,ssep+2);

        cout<<"separator union v1 v2:"<<v1<<""<<v2<<endl;
        for(int c=0;c<ssep+2;c++)
        {
            for(int d=0;d<ssep+2;d++)
                cout<<aux[c][d]<<" ";
            cout<<endl;
        }


    }

    for(int c=0;c<sCA;c++)
    delete[] aux[c];
        delete[] aux;
    w=(dSuv1*dSuv2)/(dC*dS);
}*/

/*
void tChJT::clean()
{
    aux= new Graph<Cherry,Separator>(true);
    weight=0.0;
    vector<int> m(tcjt->VertexCount(),-1);   ///vector for mapping clusters from the constructed cherry tree to the new and clean
    double maxWCluster = 0.0;
    for(int c=0;c<tcjt->VertexCount();c++)
    {
        if(tcjt->GetVertexData(c).exists)
        {
            Cherry cn;
            cn.active=cn.exists=true;
            for(int n : tcjt->GetVertexData(c).nodes)
                cn.nodes.push_back(n);
            cn.weight=tcjt->GetVertexData(c).weight;
            weight+=cn.weight;
            m[tcjt->GetVertex(c).GetID()]=aux->AddVertex(cn);
            if (cn.weight>maxWCluster)
            {
                maxWCluster = cn.weight;
                root = m[tcjt->GetVertex(c).GetID()];
            }
        }
    }
    //cout<<"Suma de Peso de clusters:"<<weight<<endl;
    for(int c=0;c<tcjt->VertexCount();c++)
    {
        if(tcjt->GetVertexData(c).exists)
            for(OutEdge<Separator>& s : tcjt->GetVertex(c).GetOutgoingEdges())
                if(tcjt->GetVertexData(s.GetDestID()).exists && (int)s.GetData().nodes.size()==(k-1)) ///at this point, every valid separator must have k-1 vars
                {
                    Separator as;
                    as.complete=true;
                    for(int i:s.GetData().nodes)
                        as.nodes.push_back(i);
                    if(s.GetData().complete)
                    {
                        as.v.push_back(s.GetData().v[0]);
                        as.dVA=s.GetData().dVA;
                    }
                    else
                    {
                        vector<int> ev;
                        Getdonating_V_ariable(s.GetDestID(),ev,c);
                        as.v.push_back(ev[0]);
                        Getdonating_V_ariable(c,ev,s.GetDestID());
                        as.dVA=ev[1];
                    }
                    as.weight=s.GetData().weight;
                    weight-=(as.weight/2.0);
                    aux->AddEdge(m[c],m[s.GetDestID()],as,0.0,false);
                }
    }
    //root=m[root];
    if (root<0)
        cout<<"Mala raiz: " <<root<<endl;
    delete tcjt;
    tcjt=aux;
    aux=nullptr;
    while(!pq.empty())
        pq.pop();
}*/



///-----Order Update










ostream & operator<<(ostream & out, tChJT & g)
{

    g.Print(out);
    return out;
}


void tChJT::Print(ostream & out)
{
    out << "W= "<<getWeight()<<endl<<"Clusters: \n";
    for(int c=0;c<tchjt.VertexCount();c++)
    {
        if(!(tchjt.GetVertex(c).GetData().deleted))
        {
            out<<"# "<<tchjt.GetVertex(c).GetID()<<":( ";
            for (int d : tchjt.GetVertex(c).GetData().nodes)
                out << d << " ";
            out<<"), W= "<<tchjt.GetVertex(c).GetData().weight<<"\n";
        }

    }

    out << "\n\n";
    out << "Out Separators: \n";
    for(int c=0;c<tchjt.VertexCount();c++)
    {
        if(!(tchjt.GetVertex(c).GetData().deleted))
        {
            Vertex<Cluster,Separator> v=tchjt.GetVertex(c);
            out << v.GetID() << "-> {";
            for (Edge<Separator> e : v.GetOutgoingEdges())
            {
                out << " " << tchjt.GetVertex(e.GetDestinationID()).GetID()<<" - C:"<<e.GetData().complete<<" - DV:"<<e.GetData().Xu<<" - V:";
                for(int nod : e.GetData().Xv)
                    out << "," << nod;
                out<<" - NODES:(";
                for (int nod : e.GetData().nodes)
                    out << " " << nod;
                out<<")| ";
            }
            out << "}\n";
        }
    }
    out << "\n\n";
    out << "In Separators \n";
    for(int c=0;c<tchjt.VertexCount();c++)
    {
        if(!(tchjt.GetVertex(c).GetData().deleted))
        {
            Vertex<Cluster,Separator> v=tchjt.GetVertex(c);
            out << v.GetID() << "<-{";
            for (Edge<Separator> e : v.GetIngoingEdges())
            {
                out << " " << tchjt.GetVertex(e.GetOriginID()).GetID()<<"- C:"<<e.GetData().complete<<" - DV:"<<e.GetData().Xu<<" - V:";
                for(int nod : e.GetData().Xv)
                    out << "," << nod;
                out<<" - NODES: (";
                for (int nod : e.GetData().nodes)
                    out << " " << nod;
                out<<")| ";
            }
            out << "}\n";
        }
    }
    out << "\n";
}
