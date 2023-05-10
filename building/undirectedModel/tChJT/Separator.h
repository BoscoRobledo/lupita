#ifndef SC_SEPARATOR_H_INCLUDED
#define SC_SEPARATOR_H_INCLUDED
struct Separator{
    vector<int> nodes;
    vector<int> Xv;
    int Xu;
    bool complete;
    double weight;
    double diffEntropy;
    double corrMDet;
    double covMDet;
    vector<double*> corrMCholDec;
    vector<double*> covMCholDec;
    inline bool operator<(const Separator& comp) const
    {
        return weight < comp.weight;
    }
};
#endif
