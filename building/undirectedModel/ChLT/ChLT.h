#ifndef CHLT_HPP_INCLUDED
#define CHLT_HPP_INCLUDED
#include "../../../graph/Graph.h"
#include "../UndirectedModel.h"
class ChLT : public UndirectedModel
{
    public:
        ChLT(double** corrM, double** covM, int d): UndirectedModel(corrM,covM,d,0){}
        void build()
        {
            Graph<int,double> aux(false);
            //build the complete graph from a full correlation matrix
            for(int c=0;c<d;c++)
              aux.AddVertex(c);
            for(int c=0;c<d;c++)
              for(int b=0;b<d;b++)
                  if(c!=b)
                      //Maximizing rho^2 is equivalent for maximizing gaussian MI
                      aux.AddEdge(c, b, corrM[c][b]*corrM[c][b], false);

            //Get MST from the auxiliar graph previously constructed.
            Edge<double> mx=structure.FillMWST(aux);
            if(covM[mx.GetOriginID()][mx.GetOriginID()]<=covM[mx.GetDestinationID()][mx.GetDestinationID()])
              bestV=mx.GetOriginID();
            else
              bestV=mx.GetDestinationID();
            k=2;
        }

        friend ostream & operator<<(ostream & out, ChLT & g)
        {

            out<<"Best Vertex="<<g.getBestVertex()<<endl<<g.structure;
            return out;
        }
};

#endif // CHLT_HPP_INCLUDED
