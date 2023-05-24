#ifndef UFDS_HPP_INCLUDED
#define UFDS_HPP_INCLUDED
#include <vector>

using namespace std;

/** \brief Union Find Disjoint Set Algorithm. Code taken from https://sort-care.github.io/Disjoint-Sets/
 * \author Juan Bosco Robledo Muñoz, 2017
 */
class UFDS
{
     private:

        /** \brief Vertexes vector
         *
         */
        vector<int> p;

        /** \brief Rank aux Vector
         *
         */
        vector<int> rank;

    public:

        /** \brief Creates disjoint set data structure for algorithm
         *
         * \param N int size of data structure
         *
         */
        UFDS(int N);

        /** \brief Recursively find the set with vertex i
         *
         * \param i int vertex index
         * \return int set of vertex i
         *
         */
        int findSet(int i) { return (p[i] == i) ? i : ( p[i] = findSet(p[i])); }


        /** \brief Defines if vertexes i and j are in same set.
         *
         * \param i int vertex index
         * \param j int vertex index
         * \return bool true if vertexes i and j are in same set
         *
         */
        bool isSameSet(int i, int j) { return findSet(i) == findSet(j); }

        /** \brief Ranks union set for vertexes i and j.
         *
         * \param i int vertex index
         * \param j int vertex index
         * \return void
         *
         */
        void unionSet(int i, int j);
};

#endif // UFDS_HPP_INCLUDED
