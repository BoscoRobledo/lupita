#ifndef POTENTIALUPDATE_H_INCLUDED
#define POTENTIALUPDATE_H_INCLUDED

struct PotentialUpdate{
    //Total correlation changes
    double wInc;
    //store values for equation 10 terms,
    //useful for keeping clusters and separators weights updated
    double wSXv;
    double wSXvXu;
    double wSXu;
    double wS;

    //Differential entropy changes
    double hDec;
    double hSXv;
    double hSXvXu;
    double hSXu;
    double hS;

    int Xv;
    int Xu;
    int CA;
    int CD;
    short updateCase;

    //Stores a copy of Cholesky decompositions for determinant calculations
    //while evaluating function 10 in order to avoid recalculations once a potential
    //step is selected
    vector<double*> covMCholDec;
    vector<double*> corrMCholDec;

    inline bool operator<(const PotentialUpdate& comp) const
    {
        return wInc < comp.wInc;
    }
};
#endif // POTENTIALUPDATE_H_INCLUDED
