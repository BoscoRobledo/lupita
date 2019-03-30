#ifndef TCHLIB_HPP_INCLUDED
#define TCHLIB_HPP_INCLUDED


#include "Graph.h"
#include <queue>
#include "Determinant.h"
#include <iostream>


struct ScTableRow{
    double w;
    int d;
    vector<int> sep;
    inline bool operator<(const ScTableRow& comp) const
    {
        return w > comp.w;
    }
};

struct Cherry{
    vector<int> nodes;
    bool active;
    bool exists;
    double weight;
};

struct Separator{
    vector<int> nodes;
    vector<int> v;
    int dVA;
    bool complete;
    double weight;
};

class TCHLib
{
    public:
        TCHLib(double** corrM, double** covM, int n, int k, bool doDAG);  //Construction from scratch
        ~TCHLib();
        const int getOrder() const;
        const int getN() const;
        double** getcovM();
        double** getcorrM();
        const vector<int>& getPerm() const;
        Graph<int,int>* getDAG();
        void buildDAG();


        friend ostream & operator<<(ostream & out, TCHLib & g);

    protected:


    private:
        Graph<Cherry,Separator>* tcjt;
        Graph<int,int>* DAG;
        vector<int> perm;

        //Construction from scratch
        vector<ScTableRow> sc;
        void ScFillTable();
        int addbyDFS(int vertex_id,vector<bool> & visited, ScTableRow & t,vector<bool> & vAdded);

        //Generate the directed acyclic graph
        void DAGDFS(int vertex_id,vector<bool> & visited);

        void Print(ostream & out);
        int k,n;
        double** corrM;
        double** covM;
};

#endif // TCHLIB_HPP_INCLUDED
