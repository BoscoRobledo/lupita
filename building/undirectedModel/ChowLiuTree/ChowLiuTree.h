#ifndef CHLT_HPP_INCLUDED
#define CHLT_HPP_INCLUDED
#include "../../../graph/Graph.h"
#include "../UndirectedModel.h"

/*
 * Class for fitting gaussian model in an Chow&Liu Tree
 * Juan Bosco Robledo Muñoz - 2017
 */
class ChowLiuTree : public UndirectedModel
{
    public:
        ChowLiuTree(double** corrM, double** covM, int d): UndirectedModel(corrM,covM,d,2){}

        /*
         * Builds the Chow&Liu Tree
         */
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

            //Set Main Vertex to min variance vertex in max correlation pair.
            if(covM[mx.GetOriginID()][mx.GetOriginID()]<=covM[mx.GetDestinationID()][mx.GetDestinationID()])
              mainV=mx.GetOriginID();
            else
              mainV=mx.GetDestinationID();
        }
};

#endif // CHLT_HPP_INCLUDED
