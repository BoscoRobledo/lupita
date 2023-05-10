#ifndef GRAPH_H
#define GRAPH_H
#include "Vertex.h"
#include "Edge.h"
#include "UFDS.h"
#include <utility>
#include <ostream>
#include <vector>
#include <cassert>
#include <queue>
#include <stack>
#include <algorithm>
#include <iostream>

using namespace std;
template <typename vD, typename eD>
class Graph
{
    private:
        vector<Vertex<vD,eD>> vertices;
        //vector<Edge<eD>> edges;

        const bool isDirected;
        const bool allowRepeatedEdges;

        void DFSRecursive(int vertex_id, vector<int> & visit_order, vector<bool> & visited);
        void Print(ostream & out);
        void FillDAGbyDFSRecursive(int vertex_id, vector<bool> &visited, Graph<vD,eD> &source);

    public:
        //Graph construction
        Graph(bool isDirected = false, bool allowRepeatedEdges = false) : isDirected(isDirected), allowRepeatedEdges(allowRepeatedEdges) {}
        int AddVertex(vD value);
        bool EdgeExists(int origin_id, int destination_id);
        void AddEdge(int origin_id, int destination_id, eD edata, bool repeat=true);

        //Graph querying
        int VertexCount() const;
        vD & GetVertexData(int vertex_id);
        Vertex<vD,eD> & GetVertex(int vertex_id);
        vector<int> GetAllVertexIDs() const;
        vector<Edge<eD>> GetEdgeList(bool ordered = false);

        //Graph traversal
        vector<int> DFS(int origin_id);

        //Fill the graph from another one
        Edge<eD> FillMWST(Graph<vD,eD>& source);
        void FillDAG(int origin_id,Graph<vD,eD> &source);

        //Graph output
        template <typename U,typename V>
        friend ostream & operator<<(ostream & out, Graph<U,V> & g);

};


///-------------Graph Construction

template <typename vD, typename eD>
int Graph<vD, eD>::AddVertex(vD value)
{
    int id = VertexCount();
    vertices.push_back(Vertex<vD, eD>(id, value));
    return id;
}

template <typename vD, typename eD>
bool Graph<vD, eD>::EdgeExists(int origin_id, int destination_id)
{   //ToDo: This can be done in O(1)
    for (Edge<eD> e : vertices[origin_id].GetOutgoingEdges())
        if(e.GetDestinationID()==destination_id)
            return true;
    return false;
}

template <typename vD, typename eD>
void Graph<vD, eD>::AddEdge(int origin_id, int destination_id, eD edata, bool repeat)
{
    assert(origin_id >= 0 && origin_id < VertexCount());
    assert(destination_id >= 0 && destination_id < VertexCount());
    if(!allowRepeatedEdges)
        if(EdgeExists(origin_id,destination_id))
           return;
    vertices[origin_id].AddOutgoingEdge(destination_id, edata);
    vertices[destination_id].AddIngoingEdge(origin_id, edata);
    //Edges are copied in reverse direction for undirected graphs
    if (!isDirected)
    {
        vertices[destination_id].AddOutgoingEdge(origin_id, edata);
        vertices[origin_id].AddIngoingEdge(destination_id, edata);
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
        vertex_ids[i] = vertices[i].GetID();
    return vertex_ids;
}


template <typename vD, typename eD>
vector<Edge<eD>> Graph<vD, eD>::GetEdgeList(bool ordered)
{
    vector<Edge<eD>> lst;
    for (Vertex<vD, eD> v : vertices)
      for (Edge<eD> e : v.GetOutgoingEdges())
      {
          lst.push_back(e);
      }
    if(ordered)
    {
        sort(lst.begin(),lst.end());
        reverse(lst.begin(),lst.end());
    }
    return lst;
}

///---------------Graph Querying

///---------------Graph Traversal

template <typename vD, typename eD>
vector<int> Graph<vD, eD>::DFS(int origin_id)
{
    vector<bool> visited(VertexCount(), false);
    vector<int> visit_order_recursive;
    DFSRecursive(origin_id, visit_order_recursive, visited);
    return visit_order_recursive;
}

template <typename vD, typename eD>
void Graph<vD, eD>::DFSRecursive(int vertex_id, vector<int> & visit_order,vector<bool> & visited)
{
    visited[vertex_id] = true;
    visit_order.push_back(vertex_id);
    for (Edge<eD> e : vertices[vertex_id].GetOutgoingEdges())
    {
        int neighbor_id = e.GetDestID();
        if (!visited[neighbor_id])
            DFSRecursive(neighbor_id, visit_order, visited);
    }
}

///---------------Graph Traversal



///---------------Graph construction from another one

template <typename vD, typename eD>
void Graph<vD, eD>::FillDAG(int origin_id,Graph<vD, eD> &source)
{
    vector<bool> visited(source.VertexCount(), false);
    for (Vertex<vD, eD> v : source.vertices)
        AddVertex(v.GetData());
    FillDAGbyDFSRecursive(origin_id, visited,source);
}

template <typename vD, typename eD>
void Graph<vD, eD>::FillDAGbyDFSRecursive(int vertex_id, vector<bool> & visited,Graph<vD, eD> &source)
{
    visited[vertex_id] = true;
    for (const Edge<eD> e : source.vertices[vertex_id].GetOutgoingEdges())
    {
        int neighbor_id = e.GetDestID();
        if (!visited[neighbor_id])
        {
            AddEdge(vertex_id,neighbor_id,false);
            FillDAGbyDFSRecursive(neighbor_id, visited,source);
        }
    }
}


template <typename vD, typename eD>
Edge<eD> Graph<vD, eD>::FillMWST(Graph<vD, eD> &source)
{
    vector<Edge<eD>> l = source.GetEdgeList(true);
    for (Vertex<vD, eD> v : source.vertices)
        AddVertex(v.GetData());
    UFDS UF(VertexCount());
    int cont=0;
    for (unsigned int i=0;i<l.size();i++)
    {
        Edge<eD> e = l[i];
        if (!UF.isSameSet(e.GetOriginID(), e.GetDestinationID()))
        {
            AddEdge(e.GetOriginID(),e.GetDestinationID(),false);
            UF.unionSet(e.GetOriginID(), e.GetDestinationID());
            cont++;
            if(cont==VertexCount()-1)
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
        for (Edge<eD> e : v.GetOutgoingEdges())
        {
            out << " " << vertices[e.GetDestinationID()].GetData();
        }
        out << "] ";
    }
    out << "\n";
    out << " In edges = ";
    for (Vertex<vD, eD> v : vertices)
    {
        out << "[" << v.GetData() << ":";
        for (Edge<eD> e : v.GetIngoingEdges())
        {
            out << " " << vertices[e.GetOriginID()].GetData();
        }
        out << "] ";
    }
    out << "\n";
}
#endif
