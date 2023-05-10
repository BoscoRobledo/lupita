#ifndef TCHJT_HPP_INCLUDED
#define TCHJT_HPP_INCLUDED
#include "../../../graph/Graph.h"
#include "Cluster.h"
#include "Separator.h"
#include "PotentialUpdate.h"
#include "../UndirectedModel.h"
#include "../ChLT/ChLT.h"
#include <queue>
#include <iostream>

class tChJT : public UndirectedModel
{
    public:
        tChJT(double** corrM, double** covM, int d);
        tChJT(ChLT* chowLiuTree, double** corrM, double** covM, int d);
        virtual ~tChJT();
        void build(int _k);
        double getWeight();
        double getModelDifferentialEntropy();
        friend ostream & operator<<(ostream & out, tChJT & g);
        void increaseOrder();
        int getRootID();

    private:
        int rootCluster;
        void initVars();
        double weight, diffEntropy;
        //tCherry Junction Tree is declared directed, because donating (v) and destination (d) variables
        //in separator depends on donating and active cluster (ie edge direction)
        Graph<Cluster, Separator> tchjt;
        ChLT* baseModel;

        ///-----Construction from 2 tCherry Junction Tree (Chow&Liu Tree)
        void Getdonating_V_ariable(int CD, vector<int> &s, int CA);
        void DFSChLT(int vertex_id, vector<bool> & visited, bool root, int parentCluster);

        ///-----Order Update
        bool isValid(PotentialUpdate& T);
        void setW(PotentialUpdate& gamma,Edge<Separator>& e);
            ///--Priority queue creation
        priority_queue<PotentialUpdate> Gamma;
        void FillPriorityQueue();
        void DFSPQ(int vertex_id,vector<bool> & visited);


        //vector<ScTableRow> sc;
        //void ScFillTable();
        //int addbyDFS(int vertex_id,vector<bool> & visited, ScTableRow & t,vector<bool> & vAdded);
        void Print(ostream & out);

};

#endif //TCHJT_HPP_INCLUDED
