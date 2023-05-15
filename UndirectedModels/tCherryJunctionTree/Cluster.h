#ifndef CLUSTER_H_INCLUDED
#define CLUSTER_H_INCLUDED
struct Cluster{
    vector<int> nodes;
    bool updated;
    bool deleted;
    double weight;
    double diffEntropy;
    /*vector<double*> corrMCholDec;
    vector<double*> covMCholDec;*/
};
#endif
