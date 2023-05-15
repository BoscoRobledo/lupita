#ifndef UNDIRECTED_MODEL_H
#define UNDIRECTED_MODEL_H

#include "../Graph/Graph.h"

/*
 * Class for fitting gaussian model in an undirected graph
 * Juan Bosco Robledo Muñoz - 2017
 */
class UndirectedModel
{
    public:
        UndirectedModel(double** _corrM, double** _covM, int _d, int _k);
        virtual void build(){}
        const int getOrder() const;
        const int getDimensionality() const;
        const int getMainVertex() const;
        double** getcovM();
        double** getcorrM();
        friend ostream & operator<<(ostream & out, UndirectedModel & g);
        // Structure of Graph that describes undirected model,
        // Vertex value: Variable index
        // Edge Value: correlation^2 as weight between vars in vertices connected
        Graph<int,double> structure;

    protected:
        // k: graph density order, d: number of vertices, mainV: Main Vertex (can be used as root for creating a DAG)
        int k,d,mainV;
        // Correlation Matrix from which the model will be fit
        double** corrM;
        // Covariance Matrix from which the model will be fit
        double** covM;
    private:
};

#endif // UNDIRECTED_MODEL_H
