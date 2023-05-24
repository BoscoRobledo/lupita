#ifndef TCHJT_HPP_INCLUDED
#define TCHJT_HPP_INCLUDED
#include <queue>
#include <iostream>
#include "../../Graph/Graph.h"
#include "../UndirectedModel.h"
#include "../ChowLiuTree/ChowLiuTree.h"
#include "Cluster.h"
#include "Separator.h"
#include "PotentialUpdate.h"

/** \brief Class for fitting gaussian model in a k-order t-Cherry Junction Tree
 * \author Juan Bosco Robledo Muñoz - 2017
 */
class tCherryJunctionTree : public UndirectedModel
{
    public:

        /** \brief Sets data needed to fit the model.
         *
         * \param corrM double** Correlation Matrix
         * \param covM double** Covariance Matrix
         * \param d int Size of Matrices
         * \return ChowLiuTree(double** corrM, double** covM, int d):
         *
         */
        tCherryJunctionTree(double** corrM, double** covM, int d);

        /*
         * Sets data needed to fit the model.
         *  corrM: Data correlation matrix
         *  covM: Data covariance matrix.
         *  d: Size of Matrices
         *
         */
        tCherryJunctionTree(ChowLiuTree* chowLiuTree, double** corrM, double** covM, int d);
        virtual ~tCherryJunctionTree();
        void build(int _k);
        double getWeight();
        double getModelDifferentialEntropy();
        friend ostream & operator<<(ostream & out, tCherryJunctionTree & g);
        void increaseOrder();
        int getRootID();

    private:
        int rootCluster;
        void initVars();
        double weight, diffEntropy;
        //tCherry Junction Tree is declared directed, because donating (v) and destination (d) variables
        //in separator depends on donating and active cluster (ie edge direction)
        Graph<Cluster, Separator> tchjt;
        ChowLiuTree* baseModel;

        ///-----Construction from 2 tCherry Junction Tree (Chow&Liu Tree)
        void Getdonating_V_ariable(int CD, vector<int> &s, int CA);
        void DFSChowLiuTree(int vertex_id, vector<bool> & visited, bool root, int parentCluster);

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
