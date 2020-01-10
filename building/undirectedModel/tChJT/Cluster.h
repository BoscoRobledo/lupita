#ifndef SC_CLUSTER_H_INCLUDED
#define SC_CLUSTER_H_INCLUDED
struct Cluster{
    vector<int> nodes;
    bool active;
    bool exists;
    double weight;
};
#endif
