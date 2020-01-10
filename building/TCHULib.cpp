#include "Graph.h"
#include "Determinant.h"
#include "TCHULib.h"

#include <queue>
#include <iostream>
#include <fstream>
#include <dirent.h>

///-----Properties
const int TCHULib::getOrder() const
{
    return k;
}
double TCHULib::getWeight()
{
    return weight;
}
const int TCHULib::getRoot() const
{
    return root;
}

const int TCHULib::getN() const
{
    return n;
}
const vector<int>& TCHULib::getPerm() const
{
    return perm;
}
double** TCHULib::getcovM()
{
    return covM;
}
double** TCHULib::getcorrM()
{
    return corrM;
}
Graph<int,int>* TCHULib::getDAG()
{
    return DAG;
}

Graph<Cherry,Separator>* TCHULib::getTCHJT()
{
    return tcjt;
}


///-----Properties

///-----Construction from 2 tCherry Junction Tree (Chow&Liu Tree)

void TCHULib::Getdonating_V_ariable(int CD, vector<int> &s, int CA)
{
    for (int vD : tcjt->GetVertexData(CD).nodes)
    {
        bool exists=false;
        for (int vS : tcjt->GetVertexData(CA).nodes)
            if(vS==vD)
            {
                exists=true;
                break;
            }
        if(!exists)
            s.push_back(vD);
    }
}
/*
void TCHULib::Get_D_ominatingVertex(int CA, vector<int> &s, int CD)
{
    for (int vA : tcjt->GetVertexData(CA).nodes)
    {
        bool exists=false;
        for (int vS : tcjt->GetVertexData(CD).nodes)
            if(vS==vA)
            {
                exists=true;
                break;
            }
        if(!exists)
            s.push_back(vA);
    }
}*/



void TCHULib::DFSCHLT(int vertex_id, vector<bool> & visited, Graph<int,int>& CHLT, bool root, int parentCluster)
{
    visited[vertex_id] = true;
    for (OutEdge<int> e : CHLT.GetVertex(vertex_id).GetOutgoingEdges())
    {
        int neighbor_id = e.GetDestID();
        if (!visited[neighbor_id])
        {
            Cherry cl;
            cl.nodes.push_back(vertex_id);
            cl.nodes.push_back(neighbor_id);
            cl.active=cl.exists=true;
            cl.weight=-0.5*log(1.0-corrM[vertex_id][neighbor_id]*corrM[vertex_id][neighbor_id]);
            weight+=cl.weight;
            int id=tcjt->AddVertex(cl);
            if(!root)
            {
                Separator s,s2;
                s.nodes.push_back(vertex_id);
                s2.nodes.push_back(vertex_id);
                Getdonating_V_ariable(id,s.v, parentCluster);
                Getdonating_V_ariable(parentCluster,s2.v,id);
                s.dVA=s2.v[0];
                s2.dVA=s.v[0];
                s.complete=s2.complete=true;
                s.weight=s2.weight=0.0;
                tcjt->AddEdge(parentCluster,id,s);
                tcjt->AddEdge(id,parentCluster,s2);
            }
            else
                parentCluster=id;
            root=false;
            DFSCHLT(neighbor_id,visited, CHLT,false,id);
        }
    }
}


TCHULib::TCHULib(Graph<int, int>& CHLT, int initNode, double** correlation, double** covariance, bool doDAG)
{
    k=2;
    n=CHLT.VertexCount();
    corrM=correlation;
    covM=covariance;
    root=0;
    tcjt= new Graph<Cherry,Separator>(true);
    DAG=new Graph<int,int>(true);
    vector<bool> visited(CHLT.VertexCount(),false);
    weight=0.0;
    DFSCHLT(initNode,visited,CHLT,true,-1);
    if(doDAG)
        buildDAG();
}


TCHULib::~TCHULib()
{
    delete tcjt;
    delete DAG;
}

///-----Construction from 2 tCherry Junction Tree (Chow&Liu Tree)

///-----Order Update


bool TCHULib::isValid(UTableRow& T)
{
    if(!(tcjt->GetVertex(T.CA).GetData().active) || !(tcjt->GetVertex(T.CD).GetData().exists))    ///If donor cluster does not exist or active cluster has already k+1 vars, continue
        return false;
    if((int)(tcjt->GetVertex(T.CA).GetData().nodes.size())>k && (int)(tcjt->GetVertex(T.CD).GetData().nodes.size())>k)  ///If both clusters size are k+1 then no potential update step exists
        return false;
    if(T.ucase==4)
        if((int)(tcjt->GetVertex(T.CA).GetData().nodes.size())>k || (int)(tcjt->GetVertex(T.CD).GetData().nodes.size())>k) ///If weight increase was calculated from case 4, this configuration has to be stiil in the tree
            return false;
    return true;
}


void TCHULib::increaseOrder()
{
    FillTable();
    /*priority_queue<UTableRow> ax;
    while(!pq.empty())
    {
        UTableRow t=pq.top();
        ax.push(t);
        pq.pop();
        cout<<t.CA<<" "<<t.CD<<" "<<" "<<t.v<<" "<<t.w<<endl;
    }
    while(!ax.empty())
    {
        UTableRow t=ax.top();
        pq.push(t);
        ax.pop();
    }*/

    int activeCherries=tcjt->VertexCount();
    bool rootUpdated=false;
    while(activeCherries>0)
    {
        UTableRow T=pq.top();
        pq.pop();
        if(isValid(T))
        {                                                           ///If maximum weight change is valid, add var v to active cluster
            tcjt->GetVertex(T.CA).GetData().nodes.push_back(T.v);
            tcjt->GetVertex(T.CA).GetData().active=false;
            tcjt->GetVertex(T.CA).GetData().weight=T.wCAn;
            activeCherries--;
            if(!rootUpdated)
            {
                root=T.CA;
                rootUpdated=true;
            }
            if(T.ucase==4)                                              ///If case 4 happens, drop donor cluster
            {
                for (OutEdge<Separator>& e : tcjt->GetVertex(T.CD).GetOutgoingEdges())       ///For each neighbor of the donor cluster
                    if(tcjt->GetVertex(e.GetDestID()).GetData().exists && e.GetDestID()!=T.CA)
                    {
                        e.GetData().complete=false;
                        tcjt->AddEdge(e.GetDestID(),T.CA,e.GetData(),0.0,false);  ///If case 4 happens, new separators are always of size k-1
                        tcjt->AddEdge(T.CA,e.GetDestID(),e.GetData(),0.0,false);   ///So keeping the old ones works fine
                    }
                    tcjt->GetVertex(T.CD).GetData().active=tcjt->GetVertex(T.CD).GetData().exists=false;
                    activeCherries--;
            }
            else if(T.ucase==2)
            {

                for (OutEdge<Separator>& e : tcjt->GetVertex(T.CA).GetOutgoingEdges())       ///Find the separators between CA and CD and add the variable v
                    if(e.GetDestID()==T.CD)
                    {
                        e.GetData().complete=false;
                        e.GetData().nodes.push_back(T.v);
                        e.GetData().weight=T.wSn;
                    }

                for (OutEdge<Separator>& e : tcjt->GetVertex(T.CD).GetOutgoingEdges())       ///Find the separator between CA and CD and add the variable v
                    if(e.GetDestID()==T.CA)
                    {
                        e.GetData().complete=false;
                        e.GetData().nodes.push_back(T.v);
                        e.GetData().weight=T.wSn;
                    }
            }

            for (OutEdge<Separator>& e : tcjt->GetVertex(T.CA).GetOutgoingEdges())       ///Add potential steps generated from adding the variable v to the neighbors of active cluster and deleting CD if happened
                if(tcjt->GetVertex(e.GetDestID()).GetData().exists && tcjt->GetVertex(e.GetDestID()).GetData().active)  ///If new or old neighbor of CA is still active, there is a potential step to add
                {
                    UTableRow newT;
                    vector<int> psa;
                    Getdonating_V_ariable(T.CA,psa,e.GetDestID());
                    for(int i=0;i<(int)psa.size();i++)
                    {
                        newT.v=psa[i];                ///Now variable v can can be donated from active cluster to its neighbors
                        newT.CA=e.GetDestID();
                        newT.ucase=2;
                        newT.CD=T.CA;             ///Active cluster is now the donor cluster
                        getW(newT,newT.CA,newT.CD,e.GetData().nodes,newT.v,newT.ucase);  ///Donor cluster (active cluster) is now size k+1, and its active neighbor cluster is size K, so case 2 always happen
                        pq.push(newT);
                    }
                }
        }
    }
    k++;
    ProcessBUDS();
    clean();
}


void TCHULib::ProcessBUDS()
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



void TCHULib::DFSGetBUDS(int vertex_id, vector<bool>& visited, vector<pair<int,pair<int,vector<int>>>>& buds)
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
}




void TCHULib::FillTable()
{
    vector<bool> visited(tcjt->VertexCount(),false);
    DFSTable(root,visited);
}


void TCHULib::DFSTable(int vertex_id,vector<bool> & visited)
{
    visited[vertex_id] = true;
    for (OutEdge<Separator> e : tcjt->GetVertex(vertex_id).GetOutgoingEdges())
    {
        int neighbor_id = e.GetDestID();
        if (!visited[neighbor_id])
        {
            UTableRow t,t2;
            t.CA=t2.CD=vertex_id;
            t.CD=t2.CA=neighbor_id;
            t.v=e.GetData().v[0];
            t2.v=e.GetData().dVA;
            t.ucase=t2.ucase=4;
            getW(t,vertex_id,neighbor_id,e.GetData().nodes,e.GetData().v[0],t.ucase);
            t2.w=t.w;
            t2.wCAn=t.wCAn;
            pq.push(t);
            pq.push(t2);
            DFSTable(neighbor_id,visited);
        }
    }
}


void TCHULib::getW(UTableRow &T, int CA, int CD, vector<int> &sep, int v, int ucase)
{
    if(ucase==2)
    {
        double dCAv,dCA, dSuv, dS;
        int sCA=tcjt->GetVertex(CA).GetData().nodes.size(), ssep=sep.size();
        double** aux=new double*[sCA+1];
        for(int c=0;c<=sCA;c++)
            aux[c]=new double[sCA+1];
        //separator
        for(int c=0;c<ssep;c++)
                for(int d=0;d<ssep;d++)
                    aux[c][d]=corrM[sep[c]][sep[d]];
        dS=determinant(aux,ssep);
        //dS=(2.0*3.14159265359*2.71828182846*covM[sep[0]][sep[0]]);

        //separator union v
        for(int c=0;c<ssep;c++)
            aux[ssep][c]=aux[c][ssep]=corrM[sep[c]][v];
        aux[ssep][ssep]=corrM[v][v];
        dSuv=determinant(aux,ssep+1);
        T.wSn=-0.5*log(dSuv);
        //Active Cluster CA
        for(int c=0;c<sCA;c++)
            for(int d=0;d<sCA;d++)
                aux[c][d]=corrM[tcjt->GetVertex(CA).GetData().nodes[c]][tcjt->GetVertex(CA).GetData().nodes[d]];
        dCA=determinant(aux,sCA);

        //Active Cluster CA union v
        for(int c=0;c<sCA;c++)
            aux[sCA][c]=aux[c][sCA]=corrM[tcjt->GetVertex(CA).GetData().nodes[c]][v];
        aux[sCA][sCA]=corrM[v][v];
        dCAv=determinant(aux,sCA+1);
        T.wCAn=-0.5*log(dCAv);
        for(int c=0;c<=sCA;c++)
            delete[] aux[c];
        delete[] aux;
        T.w=(dCA*dSuv)/(dCAv*dS);
    }
    else
    {
        double dCAv,dCA, dCD, dS;
        int sCA=tcjt->GetVertex(CA).GetData().nodes.size(),sCD=tcjt->GetVertex(CD).GetData().nodes.size(), ssep=sep.size();
        int ms=max(sCA,sCD);

        double** aux=new double*[ms+1];
        for(int c=0;c<=ms;c++)
            aux[c]=new double[ms+1];
        //separator
        for(int c=0;c<ssep;c++)
            for(int d=0;d<ssep;d++)
                aux[c][d]=corrM[sep[c]][sep[d]];
        dS=determinant(aux,ssep);

        //Donor Cluster CD
        for(int c=0;c<sCD;c++)
            for(int d=0;d<sCD;d++)
            {
                if(tcjt->GetVertex(CD).GetData().nodes[c]>=getN() || tcjt->GetVertex(CD).GetData().nodes[c]<0 || tcjt->GetVertex(CD).GetData().nodes[d]>=getN() || tcjt->GetVertex(CD).GetData().nodes[d]<0)
                {
                    FILE* dump;
                    dump=fopen("dump","wb");
                    for(int x=0;x<getN();x++)
                    {
                        for(int y=0;y<getN();y++)
                        {
                            printf("%.8lf\t",getcovM()[x][y]);
                            fwrite(&(getcovM()[x][y]),sizeof(double), 1, dump);
                        }
                        printf("\n");
                    }


                    fclose(dump);
                }
                aux[c][d]=corrM[tcjt->GetVertex(CD).GetData().nodes[c]][tcjt->GetVertex(CD).GetData().nodes[d]];
            }

        dCD=determinant(aux,sCD);

        //Active Cluster CA
        for(int c=0;c<sCA;c++)
            for(int d=0;d<sCA;d++)
                aux[c][d]=corrM[tcjt->GetVertex(CA).GetData().nodes[c]][tcjt->GetVertex(CA).GetData().nodes[d]];
        dCA=determinant(aux,sCA);

        //Active Cluster CA union v
        for(int c=0;c<sCA;c++)
            aux[sCA][c]=aux[c][sCA]=corrM[tcjt->GetVertex(CA).GetData().nodes[c]][v];
        aux[sCA][sCA]=corrM[v][v];
        dCAv=determinant(aux,sCA+1);
        T.wCAn=-0.5*log(dCAv);
        for(int c=0;c<=ms;c++)
            delete[] aux[c];
        delete[] aux;
        T.w=(dCA*dCD)/(dCAv*dS);
    }
}

void TCHULib::getWBUD(double &w, double &wS1, double &wS2, double &wC12, vector<int>& sep, int v1, int v2)
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
}


void TCHULib::clean()
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
}



///-----Order Update







///-----Building the DAG
void TCHULib::buildDAG()  ///Does a DFS in the cherry tree, creates a complete graph for the root, and adds for every child its dominating vertex
{
    vector<bool> visited(tcjt->VertexCount(),false);
    for(int c=0;c<n;c++)
        DAG->AddVertex(c);
    for(int c=0;c<k;c++)
    {
        perm.push_back(tcjt->GetVertex(root).GetData().nodes[c]);
        for(int d=0;d<c;d++)
            DAG->AddEdge(tcjt->GetVertex(root).GetData().nodes[d],tcjt->GetVertex(root).GetData().nodes[c],0,0.0,false);
    }
    //cout<<*DAG;
    DAGDFS(root,visited);
}



void TCHULib::DAGDFS(int vertex_id,vector<bool> & visited)
{
    visited[vertex_id] = true;
    for (OutEdge<Separator> e : tcjt->GetVertex(vertex_id).GetOutgoingEdges())
    {
        int neighbor_id = e.GetDestID();
        if (!visited[neighbor_id])
        {
            for(int origin : e.GetData().nodes)
                DAG->AddEdge(origin,e.GetData().v[0],0,0.0,false); ///add edges from current the separator
            perm.push_back(e.GetData().v[0]);
            DAGDFS(neighbor_id,visited);
       }
    }
}
///-----Building the DAG

/*






*/
ostream & operator<<(ostream & out, TCHULib & g)
{

    g.Print(out);
    return out;
}

void TCHULib::Print2File(int gen, int fun, int d, int k)
{
   //get the filename
    /*string path = "/home/tizo/rClassicEDA/graphs";
    int maxID = 0;
    DIR *dir;
    struct dirent *DirEntry;
    dir = opendir(path.c_str());

    while(DirEntry=readdir(dir))
    {
        string file = DirEntry->d_name;
        string strid = file.replace(file.begin(),file.begin()+1,"");
        char* p=0;
        strtol(strid.c_str(), &p, 10);
        if(!*p)
            if (maxID<atoi(strid.c_str()))
                maxID = atoi(strid.c_str());
    }
    closedir(dir);*/

    ofstream myfile;
    myfile.open ("/home/tizo/rClassicEDA/graphs/C" + to_string(d) + "_" + to_string(k) + "_" + to_string(fun) + "_" + to_string(gen));
    myfile<<"cluster,label,w_s"<<endl;
    for(int c=0;c<tcjt->VertexCount();c++)
    {
        if(tcjt->GetVertex(c).GetData().exists)
        {
            myfile<<"C_"<<tcjt->GetVertex(c).GetID()<<",{";
            for (int d : tcjt->GetVertex(c).GetData().nodes)
                myfile << d << " ";
            myfile<<"},"<<tcjt->GetVertex(c).GetData().weight<<"\n";
        }

    }
    myfile.close();

    myfile.open ("/home/tizo/rClassicEDA/graphs/S" + to_string(d) + "_" + to_string(k) + "_" + to_string(fun) + "_" + to_string(gen));
    myfile<<"from,to,label,w_s"<<endl;
    for(int c=0;c<tcjt->VertexCount();c++)
    {
        if(tcjt->GetVertex(c).GetData().exists)
        {
            Vertex<Cherry,Separator> v=tcjt->GetVertex(c);
            for (OutEdge<Separator> e : v.GetOutgoingEdges())
            {
                myfile << "C_" << v.GetID() << ",C_";
                myfile << tcjt->GetVertex(e.GetDestID()).GetID()<<",{";
                for (int nod : e.GetData().nodes)
                    myfile << " " << nod;
                myfile<<"},"<<e.GetData().weight<<endl;
            }
        }
    }
    myfile.close();
}

void TCHULib::Print(ostream & out)
{
    out << "W= "<<getWeight()<<endl<<"Clusters: \n";
    for(int c=0;c<tcjt->VertexCount();c++)
    {
        if(tcjt->GetVertex(c).GetData().exists)
        {
            out<<"# "<<tcjt->GetVertex(c).GetID()<<":( ";
            for (int d : tcjt->GetVertex(c).GetData().nodes)
                out << d << " ";
            out<<"), W= "<<tcjt->GetVertex(c).GetData().weight<<"\n";
        }

    }

    out << "\n\n";
    out << "Out Separators: \n";
    for(int c=0;c<tcjt->VertexCount();c++)
    {
        if(tcjt->GetVertex(c).GetData().exists)
        {
            Vertex<Cherry,Separator> v=tcjt->GetVertex(c);
            out << v.GetID() << "-> {";
            for (OutEdge<Separator> e : v.GetOutgoingEdges())
            {
                out << " " << tcjt->GetVertex(e.GetDestID()).GetID()<<" - C:"<<e.GetData().complete<<" - DV:"<<e.GetData().dVA<<" - V:";
                for(int nod : e.GetData().v)
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
    for(int c=0;c<tcjt->VertexCount();c++)
    {
        if(tcjt->GetVertex(c).GetData().exists)
        {
            Vertex<Cherry,Separator> v=tcjt->GetVertex(c);
            out << v.GetID() << "<-{";
            for (InEdge<Separator> e : v.GetIngoingEdges())
            {
                out << " " << tcjt->GetVertex(e.GetOriginID()).GetID()<<"- C:"<<e.GetData().complete<<" - DV:"<<e.GetData().dVA<<" - V:";
                for(int nod : e.GetData().v)
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
