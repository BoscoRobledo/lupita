#include "UFDS.h"

/** \brief Creates disjoint set data structure for algorithm
 *
 * \param N int size of data structure
 *
 */
UFDS::UFDS(int N)
{
    rank.assign(N, 0);
    p.assign(N, 0);
    for (int i=0;i<N;i++)
        p[i]=i;
}

/** \brief Ranks union set for vertexes i and j.
 *
 * \param i int vertex index
 * \param j int vertex index
 * \return void
 *
 */
void UFDS::unionSet(int i, int j)
{
    if (!isSameSet(i, j))
    {
        int x = findSet(i), y = findSet(j);
        if (rank[x] > rank[y])
            p[y] = x;
        else
        {
            p[x] = y;
            if (rank[x] == rank[y])
            rank[y]++;
        }
    }
}
