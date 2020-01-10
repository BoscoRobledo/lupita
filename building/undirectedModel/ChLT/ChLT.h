#ifndef CHLT_HPP_INCLUDED
#define CHLT_HPP_INCLUDED

#include "../../../graph/Graph.h"
#include "../UndirectedModel.h"
class ChLT : public UndirectedModel
{
    public:
        ChLT(double** corrM, double** covM, int d, int k): UndirectedModel(corrM,covM,d,k){}
        void build()
        {
            Graph<int,int> aux(false);
            //build the complete graph from a full correlation matrix
            for(int c=0;c<d;c++)
              aux.AddVertex(c);
            for(int c=0;c<d;c++)
              for(int b=0;b<d;b++)
                  if(c!=b)
                      //Maximizing rho^2 is equivalent for maximizing gaussian MI
                      aux.AddEdge(c, b, b, corrM[c][b]*corrM[c][b], false);

            //Get MST from the auxiliar graph previously constructed.
            Edge mx=structure.fillMSTFromGraph(aux);
            if(covM[mx.GetStartID()][mx.GetStartID()]<=covM[mx.GetDestID()][mx.GetDestID()])
              bestV=mx.GetStartID();
            else
              bestV=mx.GetDestID();
        }

        friend ostream & operator<<(ostream & out, ChLT & g)
        {
            out<<g.structure;
            return out;
        }
};

#endif // CHLT_HPP_INCLUDED
