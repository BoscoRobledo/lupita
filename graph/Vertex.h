#ifndef VERTEX_HPP_INCLUDED
#define VERTEX_HPP_INCLUDED
#include "OutEdge.h"
#include "InEdge.h"
#include <utility>
#include <ostream>
#include <vector>
using namespace std;
template <typename vD, typename eD>
class Vertex
{
    public:
        Vertex(int id, vD value): id(id), data(value) {}
        void AddOutgoingEdge(int end_id, eD edata);
        void AddIngoingEdge(int origin_id, eD edata);
        const int GetID() const;

        vD & GetData();
        vector<OutEdge<eD>> & GetOutgoingEdges();
        vector<InEdge<eD>> & GetIngoingEdges();

    private:
        int id; // unique identifier
        vD data;
        vector<OutEdge<eD>> outgoing_edges;
        vector<InEdge<eD>> ingoing_edges;
};

template <typename vD, typename eD>
void Vertex<vD,eD>::AddOutgoingEdge(int end_id, eD edata)
{
    outgoing_edges.push_back(OutEdge<eD>(end_id, edata));
}

template <typename vD, typename eD>
void Vertex<vD,eD>::AddIngoingEdge(int origin_id, eD edata)
{
    ingoing_edges.push_back(InEdge<eD>(origin_id, edata));
}

template <typename vD, typename eD>
vD & Vertex<vD,eD>::GetData()
{
    return data;
}

template <typename vD, typename eD>
const int Vertex<vD,eD>::GetID() const
{
    return id;
}

template <typename vD, typename eD>
vector<OutEdge<eD>> & Vertex<vD,eD>::GetOutgoingEdges()
{
    return outgoing_edges;
}

template <typename vD, typename eD>
vector<InEdge<eD>> & Vertex<vD,eD>::GetIngoingEdges()
{
    return ingoing_edges;
}

#endif // VERTEX_HPP_INCLUDED
