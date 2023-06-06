#include "../../Graph/Graph.h"
#include "../UndirectedModel.h"
#include "ChowLiuTree.h"


//{ Private members.

/** \brief Builds the Chow&Liu Tree
*
* \return void
*
*/
void ChowLiuTree::build()
{
    Graph<int,double> aux(false);

    //build the complete graph from a full correlation matrix
    for(int c=0;c<d;c++)
      aux.AddVertex(c);
    for(int c=0;c<d;c++)
      for(int b=0;b<d;b++)
          if(c!=b)
              //Maximizing rho^2 is equivalent for maximizing gaussian MI
              aux.AddEdge(c, b, corrM[c][b]*corrM[c][b]);

    //Get MST from the auxiliar graph previously constructed.
    Edge<double> maxCorrelationPair=structure.FillMWST(aux);

    //Set Main Vertex to min variance vertex in max correlation pair.
    if(covM[maxCorrelationPair.GetOriginID()][maxCorrelationPair.GetOriginID()]<=covM[maxCorrelationPair.GetDestinationID()][maxCorrelationPair.GetDestinationID()])
      mainV=maxCorrelationPair.GetOriginID();
    else
      mainV=maxCorrelationPair.GetDestinationID();
}

//}
