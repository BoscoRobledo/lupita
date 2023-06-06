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
         *
         */
        tCherryJunctionTree(double** corrM, double** covM, int d) : UndirectedModel(corrM,covM,d,0)
        { initVars(); baseModel=new ChowLiuTree(corrM,covM,d); }

        /** \brief Sets the base model and data needed to fit higher orders.
         *
         * \param chowLiuTree ChowLiuTree* ChowLiuTree as base model.
         * \param corrM double** Correlation Matrix.
         * \param covM double** Covariance Matrix.
         * \param d int Size of Matrices.
         *
         */
        tCherryJunctionTree(ChowLiuTree* chowLiuTree, double** corrM, double** covM, int d) : UndirectedModel(corrM,covM,d,0), baseModel(chowLiuTree)
        { initVars(); }

        /** \brief Frees all resources
         *
         */
        virtual ~tCherryJunctionTree() {}

        /** \brief Builds a k-tCherry Junction tree
         *
         * \param _k int Model order (max clique size in graph structure)
         *
         */
        void build(int _k);

        /** \brief Builds a k+1 order t-Cherry Junction Tree from a k order one
         *
         * \return void
         *
         */
        void increaseOrder();

        /** \brief Gets maximum weight vertex
         *
         * \return int vertex index.
         *
         */
        int getRootID();

        /** \brief Returns the weight of k-tCherry Junction tree
         *
         * \return double Total Weight
         *
         */
        double getWeight() { return weight; }

        /** \brief Returns models differential entropy
         *
         * \return double Differential Entropy
         *
         */
        double getModelDifferentialEntropy() { return diffEntropy; }

        /** \brief Prints out t-Cherry Junction Tree as text in given stream
         *
         * \param out ostream& Destination stream
         * \param g tCherryJunctionTree& Model to print
         * \return friend ostream& Destination stream
         *
         */
        friend ostream & operator<<(ostream & out, tCherryJunctionTree & g) { g.Print(out); return out; }

    private:

        /** \brief Graph for t-Cherry Junction Tree Structure. Graph is used directed through construction, because donating (v) and destination (d) variables in separator depends on donating and active cluster (i.e. edge direction)
         *
         */
        Graph<Cluster, Separator> tchjt;

        /** \brief Chow-Liu Tree used as base model (2nd order t-Cherry Junction Tree)
         *
         */
        ChowLiuTree* baseModel;

        /** \brief Initialize cumulative properties
         *
         * \return void
         *
         */
        void initVars(){ weight=0.0; diffEntropy=0.0; k=0; }

        /** \brief Total Weight of t-Cherry Junction Tree
         *
         */
        double weight;

        /** \brief Differential entropy of t-Cherry Junction Tree
         *
         */
        double diffEntropy;

        /** \brief Maximum Weight Cluster
         *
         */
        int rootCluster;

        /** \brief Constructs in structure a 2 order t-Cherry Junction Tree from Chow&Liu Tree
         *
         * \param vertex_id int Root vertex for DFS
         * \param visited vector<bool>& visited nodes vector
         * \param root bool States if current call is in root node
         * \param parentCluster int Parent cluster id
         * \return void
         *
         */
        void DFSChowLiuTree(int vertex_id, vector<bool> & visited, bool root, int parentCluster);


        ///-----Order Update

        /** \brief Gets a donating for a potential update
         *
         * \param CD int Donor cluster
         * \param s vector<int>&
         * \param CA int Active cluster
         * \return void
         *
         */
        void Getdonating_V_ariable(int CD, vector<int> &s, int CA);

        /** \brief Checks if a potential update is feasible.
         *
         * \param T PotentialUpdate& Potential update to check
         * \return bool true if valid
         *
         */
        bool isValid(PotentialUpdate& T);

        /** \brief Sets weight and diferential entropy changes for cluster and graph.
         *
         * \param gamma PotentialUpdate& Selected potential update.
         * \param e Edge<Separator>& Edge in where to apply the potential update
         * \return void
         *
         */
        void setW(PotentialUpdate& gamma,Edge<Separator>& e);

        /** \brief Potential Update priority queue
         *
         */
        priority_queue<PotentialUpdate> Gamma;

        /** \brief Fill PU priority queue from base model
         *
         * \return void
         *
         */
        void FillPriorityQueue();

        /** \brief Do a DFS in base Model to fill PU from base Model
         *
         * \param vertex_ id int root node for DFS
         * \param visited vector<bool>& visited nodes vector
         * \return void
         *
         */
        void DFSPQ(int vertex_id,vector<bool> & visited);


        /** \brief Prints t-Cherry Junction Tree as text representation in given stream
         *
         * \param out ostream& Destination stream
         * \return void
         *
         */
        void Print(ostream & out);


        //vector<ScTableRow> sc;
        //void ScFillTable();
        //int addbyDFS(int vertex_id,vector<bool> & visited, ScTableRow & t,vector<bool> & vAdded);

};

#endif //TCHJT_HPP_INCLUDED
