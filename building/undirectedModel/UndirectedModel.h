#ifndef UNDIRECTED_MODEL_H
#define UNDIRECTED_MODEL_H

#include "../../graph/Graph.h"
class UndirectedModel
{
    public:
        UndirectedModel(double** _corrM, double** _covM, int _d, int _k);
        virtual void build(){}
        const int getOrder() const;
        const int getDimensionality() const;
        const int getBestVertex() const;
        double** getcovM();
        double** getcorrM();
        friend ostream & operator<<(ostream & out, UndirectedModel & g);
        Graph<int,double> structure;

    protected:
        int k,d,bestV;
        double** corrM;
        double** covM;
    private:


};

#endif // UNDIRECTED_MODEL_H
