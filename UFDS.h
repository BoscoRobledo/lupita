#ifndef UFDS_HPP_INCLUDED
#define UFDS_HPP_INCLUDED
using namespace std;
class UFDS
{
    private:
        vector<int> p, rank;
    public:
        UFDS(int N)
        {
            rank.assign(N, 0);
            p.assign(N, 0);
            for (int i=0;i<N;i++)
                p[i]=i;
        }
    int findSet(int i)
    {
        return (p[i] == i) ? i : (p[i] = findSet(p[i]));
    }

    bool isSameSet(int i, int j)
    {
        return findSet(i) == findSet(j);
    }
    void unionSet(int i, int j)
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
};

#endif // UFDS_HPP_INCLUDED
