#ifndef TCHJT_HPP_INCLUDED
#define TCHJT_HPP_INCLUDED


#include "Graph.h"
#include <queue>
#include "Determinant.h"
#include <iostream>

class tChJT : public UndirectedModel
{
    public:
        tChJT(double** corrM, double** covM, int d, int k);
        virtual ~tChJT();
        void build();
        friend ostream & operator<<(ostream & out, tChJT & g);
    protected:

    private:
        Graph<Cluster, Separator>* tchjt;
        vector<ScTableRow> sc;
        void ScFillTable();
        //int addbyDFS(int vertex_id,vector<bool> & visited, ScTableRow & t,vector<bool> & vAdded);
        void Print(ostream & out);

};

#endif // TCHLIB_HPP_INCLUDED
