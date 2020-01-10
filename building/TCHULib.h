#ifndef TCHULIB_H_INCLUDED
#define TCHULIB_H_INCLUDED
#include "Graph.h"
#include "TCHLib.h"
struct ScTableRow{
    double w;
    int d;
    vector<int> sep;
    inline bool operator<(const ScTableRow& comp) const
    {
        return w > comp.w;
    }
};

struct UTableRow{
    double w;
    double wCAn;
    double wSn;
    int v;
    int ucase;
    int CA;
    int CD;
    inline bool operator<(const UTableRow& comp) const
    {
        return w < comp.w;
    }
};

class TCHULib
{
    public:
        TCHULib(Graph<int, int> & CHLT, int initNode, double** corrM, double** covM, bool doDAG);
        ~TCHULib();
        double getWeight();
        const int getOrder() const;
        const int getN() const;
        double** getcovM();
        double** getcorrM();
        const int getRoot() const;
        const vector<int>& getPerm() const;
        Graph<int,int>* getDAG();
        Graph<Cherry,Separator>* getTCHJT();
        void buildDAG();

        void increaseOrder();
        friend ostream & operator<<(ostream & out, TCHULib & g);
        void Print2File(int gen, int fun, int d, int k);





    protected:


    private:
        int root;
        Graph<Cherry,Separator>* tcjt;
        Graph<Cherry,Separator>* aux;
        Graph<int,int>* DAG;
        vector<int> perm;
        //construction of the potential steps table
        priority_queue<UTableRow> pq;

        double getW(vector<int> vars);
        void FillTable();
        bool isValid(UTableRow& T);
        void Getdonating_V_ariable(int CD, vector<int> &s, int CA);
        void Get_D_ominatingVertex(int CA, vector<int> &s, int CD);
        void DFSTable(int vertex_id,vector<bool> & visited);
        void getW(UTableRow &T,int CA, int CD, vector<int> &sep,int v, int ucase);
        void DFSGetBUDS(int vertex_id, vector<bool>& visited, vector<pair<int,pair<int,vector<int>>>>& buds);
        void ProcessBUDS();
        void getWBUD(double &w, double &wS1, double &wS2, double &wC12, vector<int>& sep, int v1, int v2);
        void clean();
        //Construction from a Chow-Liu Tree

        void DFSCHLT(int vertex_id,vector<bool> & visited, Graph<int, int>& CHLT,bool root, int parentCluster);


        //Generate the directed acyclic graph



        void DAGDFS(int vertex_id,vector<bool> & visited);
        void Print(ostream & out);
        int k,n;
        double** corrM;
        double** covM;
        double weight;
};

#endif // TCHULIB_H_INCLUDED
