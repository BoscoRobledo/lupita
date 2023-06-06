#ifndef SC_SEPARATOR_H_INCLUDED
#define SC_SEPARATOR_H_INCLUDED

/** \brief Structure for separator representation
 */
struct Separator{

    /** \brief Nodes in separator
     *
     */
    vector<int> nodes;

    /** \brief Donating variables in separator
     *
     */
    vector<int> Xv;

    /** \brief Dominant variable in Active Cluster
     *
     */
    int Xu;

    /** \brief Separator size is k-1?
     *
     */
    bool complete;

    /** \brief Weight of separator
     *
     */
    double weight;

    /** \brief Separator variables differential entropy
     *
     */
    double diffEntropy;

    /** \brief Determinant of separator variables correlation matrix
     *
     */
    double corrMDet;

    /** \brief Determinant of separator variables covariance matrix
     *
     */
    double covMDet;

    /** \brief Cholesky decomposition of separator variables correlation matrix
     *
     */
    vector<double*> corrMCholDec;

    /** \brief Cholesky decomposition of separator variables covariance matrix
     *
     */
    vector<double*> covMCholDec;

    /** \brief Comparator by weights
     *
     */
    inline bool operator<(const Separator& comp) const
    {
        return weight < comp.weight;
    }
};
#endif
