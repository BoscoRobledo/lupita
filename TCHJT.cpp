#include "Graph.hpp"
#include "TCHJT.hpp"
template <typename T>
const int TCHJT<T>::getOrder() const
{
    return k;
}

template <typename T>
const vector<vector<int> >& TCHJT<T>::getClusters() const
{
    return clusters;
}
template <typename T>
const vector<vector<int> >& TCHJT<T>::getSeparators() const
{
    return seps;
}

template <typename T>
TCHJT<T>::TCHJT(Graph<T>& CHLT, int initNode) : Graph<T>(false)
{
    k=2;
    for(int c=0;c<CHLT.VertexCount();c++)
        this->AddVertex(CHLT.GetVertex(c).GetData());
    for(int c=0;c<CHLT.VertexCount();c++)
        for (const OutEdge e : CHLT.GetVertex(c).GetOutgoingEdges())
            this->AddEdge(CHLT.GetVertex(c).GetID(), e.GetDestID(), e.GetCost(),false);


    vector<bool> visited(CHLT.VertexCount(), false);
    DFSCHLT(initNode, visited, CHLT,true);


}



template <typename T>
void TCHJT<T>::DFSCHLT(int vertex_id,vector<bool> & visited, Graph<T>& CHLT,bool root)
{
    visited[vertex_id] = true;
    for (const OutEdge e : CHLT.GetVertex(vertex_id).GetOutgoingEdges())
    {
        int neighbor_id = e.GetDestID();
        if (!visited[neighbor_id])
        {
            vector<int> cl;
            cl.push_back(vertex_id);
            cl.push_back(neighbor_id);
            clusters.push_back(cl);
            if(!root)
            {
                vector<int> sp;
                sp.push_back(vertex_id);
                seps.push_back(sp);
            }
            root=false;
            DFSCHLT(neighbor_id,visited, CHLT,false);
        }
    }
}



template <typename T>
ostream & operator<<(ostream & out, TCHJT<T> & g)
{
    //g.Print(out);
    g.prnt(out);
    return out;
}



template <typename T>
void TCHJT<T>::prnt(ostream & out) const
{
    out << endl<<"Clusters = ";
    for (const vector<int> c : clusters)
    {
        out<<"[";
        for (const int d : c)
            out << d << " ";
        out<<"] ";
    }

    out << endl<<"Separators = ";
    for (const vector<int> c : seps)
    {
        out<<"[";
        for (const int d : c)
            out << d << " ";
        out<<"] ";

    }

    out << "\n";
}
