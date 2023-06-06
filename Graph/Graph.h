#ifndef GRAPH_H
#define GRAPH_H

#include <utility>
#include <ostream>
#include <vector>
#include <cassert>
#include <queue>
#include <stack>
#include <algorithm>
#include <iostream>

#include "Vertex.h"
#include "Edge.h"
#include "UFDS.h"

using namespace std;
template <typename vD, typename eD>

/** \brief Graph and traverse utilities.
  * \param type eD type for object stored in Edge.
  * \param type eD type for object stored in Edge.
  * \author Juan Bosco Robledo Muñoz, 2017
 */
class Graph
{
    public:
        // { Graph construction

        /** \brief Builds an empty graph.
         *
         * \param isDirected bool When false, adding an edge will result in both ways edges added.
         * \param allowRepeatedEdges bool When true adding more than an edge between the same origin and destination vertex is allowed.
         *
         */
        Graph(bool isDirected = false, bool allowRepeatedEdges = false) : isDirected(isDirected), allowRepeatedEdges(allowRepeatedEdges) {}

        /** \brief Add an unconnected vertex to the graph.
         *
         * \param value vD Value to store in vertex.
         * \return int Vertex index.
         *
         */
        int AddVertex(vD value);

        /** \brief Checks if edge between given vertexes exists.
         *
         * \param origin_id int origin vertex id
         * \param destination_id int destination vertex id
         * \return bool exists or not.
         *
         */
        bool EdgeExists(int origin_id, int destination_id);

        /** \brief Adds an edge between given vertexes. If isDirected = false, adds also an edge in opposite direction.
         *
         * \param origin_id int origin vertex id.
         * \param destination_id int destination vertex id.
         * \param edata eD Value to store in Edge.
         * \return void
         *
         */
        void AddEdge(int origin_id, int destination_id, eD edata);

        // }

        // { Graph querying

        /** \brief Gets number of vertexes
         *
         * \return int Number of vertexes
         *
         */
        int VertexCount() const { return vertices.size(); }

        /** \brief Gets data in vertex with given Id
         *
         * \param vertex_id int vertex id to retrieve data from
         * \return vD& vertex data
         *
         */
        vD & GetVertexData(int vertex_id) { return vertices[vertex_id].GetData(); }

        /** \brief Gets vertex by id
         *
         * \param vertex_id int vertex id
         * \return Vertex<vD,eD>& vertex
         *
         */
        Vertex<vD,eD> & GetVertex(int vertex_id) { return vertices[vertex_id]; }

        /** \brief Gets vector with vertexes id
         *
         * \return vector<int> vertexes id vector
         *
         */
        vector<int> GetAllVertexIDs() const;

        /** \brief Gets a list with all edges in graph
         *
         * \param ordered bool Edges ordered
         * \return vector<Edge<eD>> Edges List
         *
         */
        vector<Edge<eD>> GetEdgeList(bool ordered = false);

        // }

        // { Graph traversal

        /** \brief Do a DFS
         *
         * \param origin_id int Start vertex
         * \return vector<int> Vertexes ids visited in order
         *
         */
        vector<int> DFS(int origin_id);

        // }

        // { Graph construction from another one

        /** \brief Creates a MWST from the given graph
         *
         * \param source Graph<vD, eD>& Graph from which MWST is extracted
         * \return Edge<eD> Maximum Weight Edge
         *
         */
        Edge<eD> FillMWST(Graph<vD,eD>& source);

        /** \brief Creates a DAG from the given graph starting in given vertex by DFS
         *
         * \param origin_id int origin vertex id
         * \param source Graph<vD, eD>& Graph from which DAG is constructed
         * \return void
         *
         */
        void FillDAG(int origin_id,Graph<vD,eD> &source);

        // }

        // { Graph output

        template <typename U,typename V>
        /** \brief Operator for printing out graph as text
         *
         * \param out ostream& destination stream
         * \param g Graph<vD,eD>& Graph to print
         * \return ostream& destination stream
         *
         */
        friend ostream & operator<<(ostream & out, Graph<U,V> & g);

        // }

    private:

        /** \brief Vertexes set.
         *
         */
        vector<Vertex<vD,eD>> vertices;

        /** \brief Graph is directed.
         *
         */
        const bool isDirected;

        /** \brief Allow graph to have more than one directed edge between certain origin and destination vertexes.
         *
         */
        const bool allowRepeatedEdges;

        /** \brief Traverses the graph by DFS.
         *
         * \param vertex_id int Start vertex.
         * \param visit_order vector<int>& Visiting order.
         * \param visited vector<bool>& Vector for storing already visited vertexes.
         * \return void
         *
         */
        void DFSRecursive(int vertex_id, vector<int> & visit_order, vector<bool> & visited);

        /** \brief Make a Directed Acyclic Graph from an undirected one by a DFS.
         *
         * \param vertex_id int Start Vertex.
         * \param visited vector<bool>& Visiting order.
         * \param source Graph<vD, eD>& Result as DAG.
         * \return void
         *
         */
        void FillDAGbyDFSRecursive(int vertex_id, vector<bool> &visited, Graph<vD,eD> &source);

        /** \brief Prints the graph.
         *
         * \param out ostream& Stream in which the graph will be represented as text.
         * \return void
         *
         */
        void Print(ostream & out);


};


// { Graph Construction

template <typename vD, typename eD>
/** \brief Add an unconnected vertex to the graph.
 *
 * \param value vD Value to store in vertex.
 * \return int Vertex index.
 *
 */
int Graph<vD, eD>::AddVertex(vD value)
{
    int id = VertexCount();
    vertices.push_back(Vertex<vD, eD>(id, value));
    return id;
}

template <typename vD, typename eD>
/** \brief Checks if edge between given vertexes exists.
 *
 * \param origin_id int origin vertex id
 * \param destination_id int destination vertex id
 * \return bool exists or not.
 *
 */
bool Graph<vD, eD>::EdgeExists(int origin_id, int destination_id)
{   //ToDo: This can be done in O(1)
    for (Edge<eD> e : vertices[origin_id].GetOutgoingEdges())
        if(e.GetDestinationID()==destination_id)
            return true;
    return false;
}

template <typename vD, typename eD>
/** \brief Adds an edge between given vertexes. If isDirected = false, adds also an edge in opposite direction.
 *
 * \param origin_id int origin vertex id.
 * \param destination_id int destination vertex id.
 * \param edata eD Value to store in Edge.
 * \return void
 *
 */
void Graph<vD, eD>::AddEdge(int origin_id, int destination_id, eD edata)
{
    assert(origin_id >= 0 && origin_id < VertexCount());
    assert(destination_id >= 0 && destination_id < VertexCount());
    if(!allowRepeatedEdges)
        if(EdgeExists(origin_id, destination_id))
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

// }

// { Graph Querying


template <typename vD, typename eD>
/** \brief Gets vector with vertexes id
 *
 * \return vector<int> vertexes id vector
 *
 */
vector<int> Graph<vD, eD>::GetAllVertexIDs() const
{
    vector<int> vertex_ids(VertexCount());
    for (size_t i = 0; i < vertex_ids.size(); ++i)
        vertex_ids[i] = vertices[i].GetID();
    return vertex_ids;
}

template <typename vD, typename eD>
/** \brief Gets a list with all edges in graph
 *
 * \param ordered bool Edges ordered
 * \return vector<Edge<eD>> Edges List
 *
 */
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

// }

// { Graph Traversal

template <typename vD, typename eD>
/** \brief Do a DFS
 *
 * \param origin_id int Start vertex
 * \return vector<int> Vertexes ids visited in order
 *
 */
vector<int> Graph<vD, eD>::DFS(int origin_id)
{
    vector<bool> visited(VertexCount(), false);
    vector<int> visit_order_recursive;
    DFSRecursive(origin_id, visit_order_recursive, visited);
    return visit_order_recursive;
}

template <typename vD, typename eD>
/** \brief Traverses the graph by DFS.
 *
 * \param vertex_id int Start vertex.
 * \param visit_order vector<int>& Visiting order.
 * \param visited vector<bool>& Vector for storing already visited vertexes.
 * \return void
 *
 */
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

// }



// { Graph construction from another one

template <typename vD, typename eD>
/** \brief Creates a DAG from the given graph starting in given vertex by DFS
 *
 * \param origin_id int origin vertex id
 * \param source Graph<vD, eD>& Graph from which DAG is constructed
 * \return void
 *
 */
void Graph<vD, eD>::FillDAG(int origin_id,Graph<vD, eD> &source)
{
    vector<bool> visited(source.VertexCount(), false);
    for (Vertex<vD, eD> v : source.vertices)
        AddVertex(v.GetData());
    FillDAGbyDFSRecursive(origin_id, visited,source);
}

template <typename vD, typename eD>
/** \brief Traverses given graph by DFS and creates a DAG
 *
 * \param vertex_id int current vertex id
 * \param visited vector<bool> visited vertex list
 * \param source Graph<vD, eD>& Graph from which DAG is constructed
 * \return void
 *
 */
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
/** \brief Creates a MWST from the given graph
 *
 * \param source Graph<vD, eD>& Graph from which MWST is extracted
 * \return Edge<eD> Maximum Weight Edge
 *
 */
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

// }

// { Graph printing

template <typename vD, typename eD>
/** \brief Operator for printing out graph as text
 *
 * \param out ostream& destination stream
 * \param g Graph<vD,eD>& Graph to print
 * \return ostream& destination stream
 *
 */
ostream & operator<<(ostream & out, Graph<vD, eD> & g)
{
    g.Print(out);
    return out;
}

template <typename vD, typename eD>
/** \brief Prints the graph in a stream as text
 *
 * \param out ostream& stream destination
 * \return void
 *
 */
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

// }

#endif
