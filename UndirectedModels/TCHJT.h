#ifndef TCHJT_H
#define TCHJT_H

#include "Graph.hpp"

template <typename T>
class TCHJT : Graph<T>
{
    public:
        TCHJT(Graph<T> & CHLT, int initNode);
        virtual ~TCHJT() {}
        const int getOrder() const;
        const vector<vector<int> >& getClusters() const;
        const vector<vector<int> >& getSeparators() const;
        template <typename U>
        friend ostream & operator<<(ostream & out, TCHJT<U> & g);

    protected:


    private:
        void prnt(ostream & out) const;
        vector<vector<int> > clusters;
        vector<vector<int> > seps;
        void DFSCHLT(int vertex_id,vector<bool> & visited, Graph<T>& CHLT,bool root);
        int k;
};

#endif // TCHJT_H
