#ifndef GRAPH_H
#define GRAPH_H
#include "Vertex.h"
#include "OutEdge.h"
#include "Edge.h"
#include "UFDS.h"
#include <utility>
#include <ostream>
#include <vector>
#include <cassert>
#include <queue>
#include <stack>
#include <algorithm>
using namespace std;
template <typename vD, typename eD>
class Graph
{
    private:
        vector<Vertex<vD,eD>> vertices;
        const bool directed;

        void DepthFirstSearchRecursive(int vertex_id, vector<int> & visit_order,vector<bool> & visited);
        void Print(ostream & out);
        void fillDAGDFSR(int vertex_id, vector<bool> &visited,Graph<vD,eD> &source);

    public:
        Graph(bool directed = false) : directed(directed) {}                            //Graph construction
        int AddVertex(vD value);
        bool EdgeExists(int start_id, int end_id, double cost);
        void AddEdge(int start_id, int end_id, eD edata, double cost = 0.0,bool repeat=true);



        int VertexCount() const;
        vD & GetVertexData(int vertex_id);                              //Graph querying
        Vertex<vD,eD> & GetVertex(int vertex_id);
        vector<int> GetAllVertexIDs() const;
        vector<Edge> getEdgeList(bool ordered = false);


        vector<int> DepthFirstSearch(int start_id);                           //Graph traversal



        Edge fillMSTFromGraph(Graph<vD,eD>& source);                                //Fill the graph from another one
        void fillDAG(int start_id,Graph<vD,eD> &source);


        template <typename U,typename V>
        friend ostream & operator<<(ostream & out, Graph<U,V> & g);

};


///-------------Graph Construction

template <typename vD, typename eD>
int Graph<vD, eD>::AddVertex(vD value)
{
    int id = VertexCount(); // id is the index into vertices array
    vertices.push_back(Vertex<vD, eD>(id, value));
    return id;
}

template <typename vD, typename eD>
bool Graph<vD, eD>::EdgeExists(int start_id, int end_id, double cost)
{
    for (OutEdge<eD> e : vertices[start_id].GetOutgoingEdges())
        if(e.GetDestID()==end_id && e.GetCost()==cost)
            return true;
    return false;
}
template <typename vD, typename eD>
void Graph<vD, eD>::AddEdge(int start_id, int end_id, eD edata, double cost, bool repeat)
{
    if(!(start_id >= 0 && start_id < VertexCount()) || !(end_id >= 0 && end_id < VertexCount()))
        start_id=end_id;
    assert(start_id >= 0 && start_id < VertexCount());
    assert(end_id >= 0 && end_id < VertexCount());
    if(!repeat)
        if(EdgeExists(start_id,end_id,cost))
           return;
    vertices[start_id].AddOutgoingEdge(end_id, cost, edata);
    vertices[end_id].AddIngoingEdge(start_id, cost, edata);
    if (!directed)
    {
        vertices[end_id].AddOutgoingEdge(start_id, cost, edata);
        vertices[start_id].AddIngoingEdge(end_id, cost, edata);
    }
}

///------Graph Construction



///-------------Graph Querying
template <typename vD, typename eD>
int Graph<vD, eD>::VertexCount() const
{
    return vertices.size();
}

template <typename vD, typename eD>
 vD & Graph<vD, eD>::GetVertexData(int vertex_id)
{
    return vertices[vertex_id].GetData();
}

template <typename vD, typename eD>
Vertex<vD, eD> & Graph<vD, eD>::GetVertex(int vertex_id)
{
    return vertices[vertex_id];
}

template <typename vD, typename eD>
vector<int> Graph<vD, eD>::GetAllVertexIDs() const
{
    vector<int> vertex_ids(VertexCount());
    for (size_t i = 0; i < vertex_ids.size(); ++i)
        vertex_ids[i] = i;
    return vertex_ids;
}


template <typename vD, typename eD>
vector<Edge> Graph<vD, eD>::getEdgeList(bool ordered)
{
    vector<Edge> lst;
    for (Vertex<vD, eD> v : vertices)
    for (OutEdge<eD> e : v.GetOutgoingEdges())
    {
        Edge ed(v.GetID(),e.GetDestID(),e.GetCost());
        lst.push_back(ed);
    }
    if(ordered)
        sort(lst.begin(),lst.end());
    return lst;
}

///---------------Graph Querying

///---------------Graph Traversal

template <typename vD, typename eD>
vector<int> Graph<vD, eD>::DepthFirstSearch(int start_id)
{
    vector<bool> visited(VertexCount(), false);
    vector<int> visit_order_recursive;
    DepthFirstSearchRecursive(start_id, visit_order_recursive, visited);
    return visit_order_recursive;
}

template <typename vD, typename eD>
void Graph<vD, eD>::DepthFirstSearchRecursive(int vertex_id, vector<int> & visit_order,vector<bool> & visited)
{
    visited[vertex_id] = true;
    visit_order.push_back(vertex_id);
    for (OutEdge<eD> e : vertices[vertex_id].GetOutgoingEdges())
    {
        int neighbor_id = e.GetDestID();
        if (!visited[neighbor_id])
            DepthFirstSearchRecursive(neighbor_id, visit_order, visited);
    }
}

///---------------Graph Traversal



///---------------Graph construction from another one

template <typename vD, typename eD>
void Graph<vD, eD>::fillDAG(int start_id,Graph<vD, eD> &source)
{
    vector<bool> visited(source.VertexCount(), false);
    for (Vertex<vD, eD> v : source.vertices)
        AddVertex(v.GetData());
    fillDAGDFSR(start_id, visited,source);
}

template <typename vD, typename eD>
void Graph<vD, eD>::fillDAGDFSR(int vertex_id, vector<bool> & visited,Graph<vD, eD> &source)
{
    visited[vertex_id] = true;
    for (const OutEdge<eD> e : source.vertices[vertex_id].GetOutgoingEdges())
    {
        int neighbor_id = e.GetDestID();
        if (!visited[neighbor_id])
        {
            AddEdge(vertex_id,neighbor_id,e.GetCost(),false);
            fillDAGDFSR(neighbor_id, visited,source);
        }
    }
}


template <typename vD, typename eD>
Edge Graph<vD, eD>::fillMSTFromGraph(Graph<vD, eD> &source)
{
    vector<Edge> l=source.getEdgeList(true);
    for (Vertex<vD, eD> v : source.vertices)
        AddVertex(v.GetData());
    UFDS UF(VertexCount());
    for (unsigned int i=0;i<l.size();i++)
    {
        Edge e = l[i];
        int cont=0;
        if (!UF.isSameSet(e.GetStartID(), e.GetDestID()))
        {
            AddEdge(e.GetStartID(),e.GetDestID(),e.GetCost(),false);
            UF.unionSet(e.GetStartID(), e.GetDestID());
            cont++;
            if(cont==(VertexCount()-1))
                return l[0];
        }
    }
    return l[0];
}

///---------------Graph construction from another one

///---------------Graph printing
template <typename vD, typename eD>
ostream & operator<<(ostream & out, Graph<vD, eD> & g)
{
    g.Print(out);
    return out;
}

template <typename vD, typename eD>
void Graph<vD, eD>::Print(ostream & out)
{
    out << "V = ";
    for (Vertex<vD, eD> v : vertices)
        out << v.GetData() << " ";
    out << "\n";
    out << "Out Edges = ";
    for (Vertex<vD, eD> v : vertices)
    {
        out << "[" << v.GetData() << ":";
        for (OutEdge<eD> e : v.GetOutgoingEdges())
        {
            out << " " << vertices[e.GetDestID()].GetData();
        }
        out << "] ";
    }
    out << "\n";
    out << " In edges = ";
    for (Vertex<vD, eD> v : vertices)
    {
        out << "[" << v.GetData() << ":";
        for (InEdge<eD> e : v.GetIngoingEdges())
        {
            out << " " << vertices[e.GetOriginID()].GetData();
        }
        out << "] ";
    }
    out << "\n";
}
#endif
