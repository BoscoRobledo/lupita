#ifndef VERTEX_HPP_INCLUDED
#define VERTEX_HPP_INCLUDED
#include <utility>
#include <ostream>
#include <vector>
#include "Edge.h"

using namespace std;
template <typename vD, typename eD>
/** \brief Graph directed edge, containing object of template type vD
 * \param type vD type for object stored in Vertex.
 * \author Juan Bosco Robledo Muñoz, 2017
 */
class Vertex
{
    public:

        /** \brief Builds Vertex and sets its data
         *
         * \param id int Vertex id
         * \param value vD Data to store in the Vertex (typically description data)
         *
         */
        Vertex(int _id, vD _value): id(_id), data(_value) {}

        /** \brief Gets Vertex id
         *
         * \return const int Vertex id
         *
         */
        const int GetID() const { return id; }

        /** \brief Gets Data in Vertex
         *
         * \return vD& Reference to data (of type vD) in Vertex
         *
         */
        vD & GetData() { return data; }

        /** \brief Gets set of outgoing edges
         *
         * \return vector<Edge<eD>>& Reference to vector containing outgoing edges.
         *
         */
        vector<Edge<eD>> & GetOutgoingEdges() { return outgoing_edges; }

        /** \brief Gets set of ingoing edges
         *
         * \return vector<Edge<eD>>& Reference to vector containing ingoing edges.
         *
         */
        vector<Edge<eD>> & GetIngoingEdges() { return ingoing_edges; }

        /** \brief Adds an edge to a specific vertex.
         *
         * \param destination_id int destination vertex.
         * \param edata eD data in Edge
         * \return void
         *
         */
        void AddOutgoingEdge(int destination_id, eD edata) { outgoing_edges.push_back(Edge<eD>(id, destination_id, edata)); }

        /** \brief Add an edge from a specific vertex.
         *
         * \param origin_id int origin vertex
         * \param edata eD data in Edge
         * \return void
         *
         */
        void AddIngoingEdge(int origin_id, eD edata) { ingoing_edges.push_back(Edge<eD>(origin_id, id, edata)); }

    private:

        /** \brief Vertex id
         *
         */
        int id;

        /** \brief Vertex data (of template type vD)
         *
         */
        vD data;

        /** \brief Set of outgoing edges.
         *
         */
        vector<Edge<eD>> outgoing_edges;

        /** \brief Set of ingoing edges.
         *
         */
        vector<Edge<eD>> ingoing_edges;
};

#endif // VERTEX_HPP_INCLUDED
