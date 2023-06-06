#ifndef CLUSTER_H_INCLUDED
#define CLUSTER_H_INCLUDED

/** \brief Structure for Cluster representation.
 */
struct Cluster{

    /** \brief Nodes in cluster
     *
     */
    vector<int> nodes;

    /** \brief Cluster is updated?
     *
     */
    bool updated;

    /** \brief Cluster was deleted when included in a cluster that contains already all its variables
     *
     */
    bool deleted;

    /** \brief Cluster weight
     *
     */
    double weight;

    /** \brief Cluster differential entropy
     *
     */
    double diffEntropy;

    /*vector<double*> corrMCholDec;
    vector<double*> covMCholDec;*/
};
#endif
