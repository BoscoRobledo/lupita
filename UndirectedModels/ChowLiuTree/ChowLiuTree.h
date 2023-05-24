#ifndef CHLT_HPP_INCLUDED
#define CHLT_HPP_INCLUDED
#include "../UndirectedModel.h"

/** \brief Class for fitting gaussian model in an Chow&Liu Tree
 * \author Juan Bosco Robledo Muñoz - 2017
 */
class ChowLiuTree : public UndirectedModel
{
    public:

        /** \brief Sets data and fits the model.
         *
         * \param corrM double** Correlation Matrix
         * \param covM double** Covariance Matrix
         * \param d int Size of Matrices
         *
         * \return Chow&Liu Tree built from Matrices using an MWST on a complete graph with weights given by correlation squared.
         *
         */
        ChowLiuTree(double** corrM, double** covM, int d): UndirectedModel(corrM,covM,d,2) { }

    private:

        /** \brief Builds the Chow&Liu Tree
         *
         * \return void
         *
         */
        void build();
};

#endif // CHLT_HPP_INCLUDED
