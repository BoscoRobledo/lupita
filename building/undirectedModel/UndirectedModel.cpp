#include "UndirectedModel.h"


/*
 * Gets graph density order.
 */
const int UndirectedModel::getOrder() const
{
    return k;
}

/*
 * Gets number of vars in model.
 */
const int UndirectedModel::getDimensionality() const
{
    return d;
}

/*
 * Gets covariance matrix.
 */
double** UndirectedModel::getcovM()
{
    return covM;
}

/*
 * Gets correlation matrix.
 */
double** UndirectedModel::getcorrM()
{
    return corrM;
}

/*
 * Constructor. Prepare all parameters required to build the undirected model.
 * _corrM: Correlation matrix
 * _covM: Covariance matrix
 * _d: Number of vars in model
 * _k: Graph density order
 */
UndirectedModel::UndirectedModel(double** _corrM, double** _covM, int _d, int _k)
{
    corrM=_corrM;
    covM=_covM;
    d=_d;
    k=_k;
}

/*
 * Prints a text representation of Undirected Model.
 * out: Destination stream
 * g: Model to print.
 * returns: reference to destination stream
 */
ostream & operator<<(ostream & out, UndirectedModel & g)
{
  out<<"Main Vertex="<<g.getMainVertex()<<endl<<g.structure;
  return out;
}

/*
 * Gets the graph main vertex.
 */
const int UndirectedModel::getMainVertex() const
{
    return mainV;
}
