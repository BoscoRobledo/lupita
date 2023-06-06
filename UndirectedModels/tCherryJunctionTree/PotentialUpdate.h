#ifndef POTENTIALUPDATE_H_INCLUDED
#define POTENTIALUPDATE_H_INCLUDED

/** \brief Struct for potential update representation
 */
struct PotentialUpdate{

    /** \brief Weight increment related
     *
     */
    double wInc;
    //store values for equation 10 terms,
    //useful for keeping clusters and separators weights updated

    /** \brief Weight of related separator union donating variable.
     *
     */
    double wSXv;

    /** \brief Weight of related separator union donating variable union active cluster dominant variable.
     *
     */
    double wSXvXu;

    /** \brief Weight of related separator union active cluster dominant variable.
     *
     */
    double wSXu;

    /** \brief Weight of related separator.
     *
     */
    double wS;

    /** \brief Differential Entropy decrement related.
     *
     */
    double hDec;

    /** \brief Differential Entropy of variables in related separator union donating variable.
     *
     */
    double hSXv;

    /** \brief Differential Entropy of variables in related separator union donating variable union active cluster dominant variable.
     *
     */
    double hSXvXu;

    /** \brief Differential Entropy of variables in related separator union active cluster dominant variable.
     *
     */
    double hSXu;

    /** \brief Differential Entropy of variables in related separator.
     *
     */
    double hS;

    /** \brief Donating variable index.
     *
     */
    int Xv;

    /** \brief Dominant variable in active cluster index.
     *
     */
    int Xu;

    /** \brief Active cluster index.
     *
     */
    int CA;

    /** \brief Donor Cluster index.
     *
     */
    int CD;

    /** \brief Potential Update case (in terms of cluster sizes)
     *
     */
    short updateCase;

    /** \brief Stores covariance matrix Cholesky decompositions for determinant calculations in weight increment for each potential update in order to set them to separators when potential update is selected.
    *
    */
    vector<double*> covMCholDec;

    /** \brief Stores correlation matrix Cholesky decompositions for determinant calculations in weight increment for each potential update in order to set them to separators when potential update is selected.
     *
     */
    vector<double*> corrMCholDec;

    /** \brief Comparation of potential updates by weight increment
     *
     */
    inline bool operator<(const PotentialUpdate& comp) const
    {
        return wInc < comp.wInc;
    }
};
#endif // POTENTIALUPDATE_H_INCLUDED
