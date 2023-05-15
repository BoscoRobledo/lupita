#include "Graph.h"
#include "Determinant.h"
#include "TCHLib.h"
#include <queue>
#include <iostream>

///-----Properties
const int TCHLib::getOrder() const
{
    return k;
}
const int TCHLib::getN() const
{
    return n;
}
const vector<int>& TCHLib::getPerm() const
{
    return perm;
}
double** TCHLib::getcovM()
{
    return covM;
}
double** TCHLib::getcorrM()
{
    return corrM;
}
Graph<int,int>* TCHLib::getDAG()
{
    return DAG;
}
///-----Properties

///-----Construction from Scratch

void TCHLib::ScFillTable()
{
    vector<bool> v(n);
    vector<int> nodes(k);                   ///Int vector to store the actual combination
    fill(v.end() - k, v.end(), true);       ///Boolean vector needed to generate the combinations
    double** aux=new double*[k];            ///Auxiliar Memory needed for the determinant calculation
        for(int c=0;c<k;c++)
            aux[c]=new double[k];

    do
    {
        int ci=0;
        for (int i = 0; i < n; ++i)         ///Generate each combination
            if (v[i])
            {
                nodes[ci]=i;
                ci++;
            }


        for(int c=0;c<k;c++)
            for(int d=0;d<k;d++)
                aux[c][d]=corrM[nodes[c]][nodes[d]];    ///Fill aux memory

        double wC=determinant(aux,k);                   ///Cluster determinant
        for(int i=0;i<k;i++)                            ///For each possible separator in the cluster
        {
            ScTableRow t;
            for(int j=0;j<k;j++)
            {
                if(i==j)
                    t.d=nodes[j];
                else
                    t.sep.push_back(nodes[j]);
            }
            if(k==2)
                t.w=1.0/(2.0*3.14159265359*2.71828182846*covM[t.sep[0]][t.sep[0]]*wC);      ///For the special case k==2 separator size is 1, and I(X)==h(X), so the separator variance is used.
            else
            {
                for(int c=0;c<(k-1);c++)
                    for(int d=0;d<(k-1);d++)
                        aux[c][d]=corrM[t.sep[c]][t.sep[d]];    ///Compute separator determinant
                        t.w=determinant(aux,k-1)/wC;
            }
            sc.push_back(t);
        }
    }
    while (next_permutation(v.begin(), v.end()));       ///For each combination
    for(int c=0;c<k;c++)
        delete[] aux[c];
    delete[] aux;
}


TCHLib::TCHLib(double** _corrM, double** _covM, int _n, int _k, bool doDAG)
{
    tcjt= new Graph<Cherry,Separator>(false);
    DAG= new Graph<int,int>(true);

    k=_k;
    n=_n;
    corrM=_corrM;
    covM=_covM;

    ScFillTable();

    sort(sc.begin(), sc.end());

    /*for(ScTableRow & t : sc)
    {
        std::cout<<t.d<<" (";
        for(int s : t.sep)
            cout<<s<<",";
        std::cout<<") "<<t.w<<endl;
    }*/

    int vertexCount=_n, tIndx=1;
    vector<bool> vAdded(vertexCount,false);
    Cherry root;
    for(int c=0;c<(k-1);c++)
    {
        root.nodes.push_back(sc[0].sep[c]);
        vAdded[sc[0].sep[c]]=true;
    }
    root.nodes.push_back(sc[0].d);
    vAdded[sc[0].d]=true;
    vertexCount-=k;
    int cherryIndx=tcjt->AddVertex(root);
    while(vertexCount>0)
    {
        if(!(vAdded[sc[tIndx].d]))
        {
            vector<bool> visited(tcjt->VertexCount(),false);
            if(addbyDFS(cherryIndx,visited,sc[tIndx],vAdded)>-1)
                vertexCount--;
        }
        tIndx++;
    }
    if(doDAG)
        buildDAG();
    //cout<<*DAG;
}

TCHLib::~TCHLib()
{
    delete tcjt;
    delete DAG;
}


int TCHLib::addbyDFS(int vertex_id,vector<bool> & visited, ScTableRow & t,vector<bool> & vAdded)
{
    visited[vertex_id] = true;
    int subset=true;
    for(int c=0;c<(k-1);c++)  ///Check if the actual separator is a subset of a cherry already present in the TCHJT
    {
        int d=0;
        for(;d<k;d++)
            if(tcjt->GetVertex(vertex_id).GetData().nodes[d]==t.sep[c])
                break;
        if(d==k)
        {
            subset=false;
            break;
        }
    }
    if(subset)
    {
        Cherry ch;
        Separator s;
        for(int c=0;c<(k-1);c++)
        {
            ch.nodes.push_back(t.sep[c]);
            s.nodes.push_back(t.sep[c]);
        }
        ch.nodes.push_back(t.d);
        s.v.push_back(t.d);
        vAdded[t.d]=true;
        int id=tcjt->AddVertex(ch);
        tcjt->AddEdge(vertex_id,id,s,0.0,false);
        return id;
    }
    else
    {
        for (OutEdge<Separator> e : tcjt->GetVertex(vertex_id).GetOutgoingEdges())
        {
            int neighbor_id = e.GetDestID();
            if (!visited[neighbor_id])
            {
                int id=addbyDFS(neighbor_id,visited,t,vAdded);
                if(id>=0)
                    return id;
            }

        }
        return -1;
    }
}
///-----Construction from Scratch

///-----Building the DAG

void TCHLib::buildDAG()  ///Does a DFS in the cherry tree, creates a complete graph for the root, and adds for every child its dominating vertex
{
    vector<bool> visited(tcjt->VertexCount(),false);
    for(int c=0;c<n;c++)
        DAG->AddVertex(c);
    for(int c=0;c<k;c++)
    {
        perm.push_back(tcjt->GetVertex(0).GetData().nodes[c]);
        for(int d=0;d<c;d++)
            DAG->AddEdge(tcjt->GetVertex(0).GetData().nodes[d],tcjt->GetVertex(0).GetData().nodes[c],0,0.0,false);
    }
    DAGDFS(0,visited);
}



void TCHLib::DAGDFS(int vertex_id,vector<bool> & visited)
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


ostream & operator<<(ostream & out, TCHLib & g)
{

    g.Print(out);
    return out;
}


void TCHLib::Print(ostream & out)
{
    out << "Cherrys = [";
    for(int c=0;c<tcjt->VertexCount();c++)
    {
        out<<"# "<<tcjt->GetVertex(c).GetID()<<":( ";
        for (int d : tcjt->GetVertex(c).GetData().nodes)
            out << d << " ";
        out<<"), ";
    }

    out << "]\n";
    out << "Out Separators = ";
    for(int c=0;c<tcjt->VertexCount();c++)
    {
        Vertex<Cherry,Separator> v=tcjt->GetVertex(c);
        out << "[" << v.GetID() << "->";
        for (OutEdge<Separator> e : v.GetOutgoingEdges())
        {
            out << " " << tcjt->GetVertex(e.GetDestID()).GetID()<<"(";
            for (int nod : e.GetData().nodes)
                out << " " << nod;
            out<<")";
        }
        out << "] ";
    }
    out << "\n";
    out << "In Separators = ";
    for(int c=0;c<tcjt->VertexCount();c++)
    {
        Vertex<Cherry,Separator> v=tcjt->GetVertex(c);
        out << "[" << v.GetID() << "<-";
        for (InEdge<Separator> e : v.GetIngoingEdges())
        {
            out << " " << tcjt->GetVertex(e.GetOriginID()).GetID()<<"(";
            for (int nod : e.GetData().nodes)
                out << " " << nod;
            out<<")";
        }
        out << "] ";
    }
    out << "\n";
}
