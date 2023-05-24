#ifndef UNDIRECTED_MODEL_H
#define UNDIRECTED_MODEL_H

#include "../Graph/Graph.h"

/** \brief Class for fitting gaussian model in an undirected graph
 * \author Juan Bosco Robledo Muñoz - 2017
 */
class UndirectedModel
{
    public:

        /** \brief Constructor. Prepare all parameters required to build the undirected model.
         *
         * \param _corrM double** Correlation matrix
         * \param _covM double** Covariance matrix
         * \param _d int Number of vars in model
         * \param _k int Graph density order
         *
         */
        UndirectedModel(double** _corrM, double** _covM, int _d, int _k) : k(_k), d(_d), covM(_covM), corrM(_corrM) { build(); }

        /** \brief Gets undirected model order (maximum clusters # of vars).
         *
         * \return const int order
         *
         */
        const int getOrder() const { return k; };

        /** \brief Gets number of vars in model.
         *
         * \return const int order
         *
         */
        const int getDimensionality() const {return d; };

        /** \brief Gets the graph main vertex.
         *
         * \return const int Main Vertex index.
         *
         */
        const int getMainVertex() const {return mainV; };

        /** \brief Gets covariance matrix.
         *
         * \return double** Pointer to covariance matrix.
         *
         */
        double** getcovM() {return covM; };

        /** \brief Gets correlation matrix.
         *
         * \return double** Pointer to correlation matrix.
         *
         */
        double** getcorrM() {return corrM; };


        /** \brief Gets graph structure (Variable index as Vertex, Correlation^2 as Edge weight between vars in connected vertexes)
         *
         * \return Graph<int, double>& Graph structure.
         *
         */
        Graph<int, double>& getStructure() { return structure; }


        /** \brief Operator for printing a text representation of Undirected Model.
         *
         * \param out ostream& Destination stream
         * \param g UndirectedModel& Model to print.
         * \return ostream& reference to destination stream
         *
         */
        friend ostream & operator<<(ostream & out, UndirectedModel & g) { out<<"Main Vertex="<<g.getMainVertex()<<endl<<g.structure; return out; }

    protected:

        /** \brief Graph density order, i.e. Maximum Clique Size
         *
         */
        int k;

        /** \brief Number of variables
         *
         */
        int d;

        /** \brief Covariance Matrix from which the model will be fit
         *
         */
        double** covM;

        /** \brief Correlation Matrix from which the model will be fit
         *
         */
        double** corrM;

        /** \brief Main Vertex (typically a vertex with low variance or in a pair with high correlation. It can be used as root for creating a DAG)
         *
         */
        int mainV;

        /** \brief Structure of Graph that describes undirected model.
         *
         * \param int Variable index as Vertex
         * \param double Correlation^2 as Edge weight between vars in connected vertexes
         *
         */
        Graph<int,double> structure;

    private:

        /** \brief Virtual method for model fitting definition
         *
         * \return virtual void
         *
         */
        virtual void build(){}
};

#endif // UNDIRECTED_MODEL_H
