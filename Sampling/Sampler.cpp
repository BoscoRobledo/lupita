#include "TCHLib.h"
#include "TCHULib.h"
#include "Graph.h"
#include "PrecMatrix.h"
#include "Simulate.h"
#include <iostream>
#include <random>
#include <chrono>
#include <Eigen/Dense>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>


void getSamplesMarginal(double* mu, double* var, double ** sample, int D, int n)
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    std::normal_distribution<double> distribution (0.0,1.0);

    for(int c=0;c<D;c++)
        for(int i=0;i<n;i++)            ///simulate n instances
            sample[i][c]=mu[c]+sqrt(var[c])*distribution(generator);
}


void getSamplesChowLiu(double* mu, double** covM, double** corrM, double ** sample, int D, int n)
{
    Graph<int,int> G,CHLT;
    for(int c=0;c<D;c++)
        G.AddVertex(c);
    for(int c=0;c<D;c++)
        for(int d=0;d<D;d++)
            if(c!=d)
                G.AddEdge(c,d,d,corrM[c][d]*corrM[c][d],false);  ///build the complete graph from a full correlation matrix

    Edge minEdge=CHLT.fillMSTFromGraph(G);                           ///Get MST from the graph previously constructed. Undirected graph

    int root=0;
    if(covM[minEdge.GetStartID()][minEdge.GetStartID()]<=covM[minEdge.GetDestID()][minEdge.GetDestID()])  ///Pick a root and turn the graph to a DAG
        root=minEdge.GetStartID();
    else
        root=minEdge.GetDestID();
    Graph<int,int> DAGC(true);
    DAGC.fillDAG(root,CHLT);
    vector<int> AO= DAGC.DepthFirstSearch(root);       ///Do a DFS to get the ancestral ordering of the tree
    //for(int i: AO)
    //    cout<<i<<" ";
    //cout<<endl;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();      ///Standard normal distributed pseudorandom numbers generator
    std::default_random_engine generator (seed);
    std::normal_distribution<double> distribution (0.0,1.0);

    for(int c: AO)            ///Calculate for each variable conditional variance and mean
    {

        vector<InEdge<int>> parentset=DAGC.GetVertex(c).GetIngoingEdges();
        int pssize=parentset.size();
        if(pssize==0)  ///Root node. Unconditional mean and variance
        {
            for(int i=0;i<n;i++)            ///simulate n instances
                sample[i][c]=mu[c]+sqrt(covM[c][c])*distribution(generator);
        }
        else
        {
            int parent=parentset[0].GetOriginID();
            double var_c=covM[c][c]-(covM[c][parent]*covM[c][parent])/covM[parent][parent];
            for(int i=0;i<n;i++)        ///simulate n instances
            {
                double mu_c=mu[c]+(covM[c][parent]*(sample[i][parent]-mu[parent]))/covM[parent][parent];  ///conditional mean
                sample[i][c]=mu_c+sqrt(var_c)*distribution(generator);
            }
        }
    }
}

/** \brief Get samples from a Gaussian Network built from a t-cherry junction tree.
 *
 * \param mu double* Unconditional means vector
 * \param covM double** Covariance Matrix to build the Cherry tree from
 * \param corrM double** Correlation Matrix to build the Cherry tree from
 * \param sample double** Memory to fill with samples
 * \param D int Sample dimensionality
 * \param n int Sample size
 * \param k int Desired order of the Cherry tree
 * \return int Weight of the Cherry tree modeled form population
 *
 */
double getSamplesTChJT(double* mu, double** covM, double** corrM, double ** sample, int D, int n, int k, int evals, int func)
{
    //TCHLib chT(corrM,covM,D,k,true);
    int gen=evals/n;
    Graph<int,int> G,CHLT;
    for(int c=0;c<D;c++)
        G.AddVertex(c);
    for(int c=0;c<D;c++)
        for(int d=0;d<D;d++)
            if(c!=d)
                G.AddEdge(c,d,d,corrM[c][d]*corrM[c][d],false);  //build the complete graph from a full correlation matrix
    Edge mx=CHLT.fillMSTFromGraph(G);
    TCHULib chT(CHLT,mx.GetStartID(),corrM,covM,false);
    for(int i=3;i<=k;i++)   ///Getting the junction tree distribution using the order update algorithm
        chT.increaseOrder();
    if(gen % 7==0)
        chT.Print2File(gen,func,D,k);
    chT.buildDAG();
    double res=chT.getWeight();
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();      ///Standard normal distributed pseudorandom numbers generator
    std::default_random_engine generator (seed);
    std::normal_distribution<double> distribution (0.0,1.0);

    double* condVar= new double[D];        ///Conditional variances are not dependent on the value taken by the parents. So they are calculated once per node
    double** beta=new double*[D];           ///Regression coeffs beta are calculated once per node, in order to use them to compute conditional means
        for(int d=0;d<D;d++)
            beta[d]=new double[k-1];

    double* sigmaxY=new double[k-1];       ///Aux memory for computing conditional means and variances
    double* Y_x=new double[k-1];
    double** psW=new double*[k-1];         ///Memory for the covariance matrix of the parent sets
    for(int d=0;d<(k-1);d++)
        psW[d]=new double[k-1];


    vector<int> AO=chT.getPerm();
    //for(int i: AO)
     //   cout<<i<<" ";
    //cout<<endl;
    for(int c: AO)            ///Calculate for each variable conditional variances and beta coeffs from its parent set
    {

        vector<InEdge<int>> parentset=chT.getDAG()->GetVertex(c).GetIngoingEdges();
        int pssize=parentset.size();
        if(pssize==0)  ///Root node. Unconditional mean and variance
        {
            for(int i=0;i<n;i++)            ///simulate n instances
            {
                sample[i][c]=mu[c]+sqrt(chT.getcovM()[c][c])*distribution(generator);
                if(sample[i][c]!=sample[i][c])
                {
                    FILE* dump;
                    dump=fopen("dump2","wb");
                    for(int x=0;x<chT.getN();x++)
                    {
                        for(int y=0;y<chT.getN();y++)
                        {
                            printf("%.8lf\t",chT.getcovM()[x][y]);
                            fwrite(&(chT.getcovM()[x][y]),sizeof(double), 1, dump);
                        }
                        printf("\n");
                    }
                }
            }
        }
        else
        {
            vector<int> AOPS;
            for(int d : AO)       ///Non-efficient. Sort ingoing edges according to Ancestral ordering. Solution. Ensure ingoing edge is sorted when creating it.
            {
                for(InEdge<int>& e : parentset)
                    if(d==e.GetOriginID())
                    {
                        AOPS.push_back(d);
                        break;
                    }
                if((int)AOPS.size()==pssize)
                    break;
            }
            if(pssize==k)
                cout<<chT;
            recursiveW(psW,chT.getcovM(),pssize,AOPS.data()); ///Get parent set precision matrix

            for(int d=0;d<pssize;d++)
                sigmaxY[d]=chT.getcovM()[AOPS[d]][c];      ///covariance vector between the variable i and its parent set
            vectorSqMatrixProduct(sigmaxY,psW,pssize,beta[c]);   ///computing beta coeffs vector
            condVar[c]=chT.getcovM()[c][c]-vectorDotProduct(beta[c],sigmaxY,pssize); ///conditional variance per node

            for(int i=0;i<n;i++)        ///simulate n instances
            {
                for(int j=0;j<pssize;j++)
                    Y_x[j]=sample[i][AOPS[j]]-mu[AOPS[j]];      ///Values taken by parents minus parents mean
                double mu_c=mu[c]+vectorDotProduct(beta[c],Y_x,pssize);  ///conditional mean
                if(condVar[c]<0.0)
                    condVar[c]*=-1.0;
                sample[i][c]=mu_c+sqrt(condVar[c])*distribution(generator);
                if(sample[i][c]!=sample[i][c])
                {
                    FILE* dump;
                    dump=fopen("dumpSimulate","wb");
                    for(int x=0;x<chT.getN();x++)
                    {
                        for(int y=0;y<chT.getN();y++)
                        {
                            printf("%.8lf\t",chT.getcovM()[x][y]);
                            fwrite(&(chT.getcovM()[x][y]),sizeof(double), 1, dump);
                        }
                        printf("\n");
                    }
                }
            }
        }
    }

    delete[] condVar;
    for(int d=0;d<D;d++)
        delete []beta[d];
    delete [] beta;
    delete [] sigmaxY;
    delete [] Y_x;
    for(int d=0;d<(k-1);d++)
        delete [] psW[d];
    delete [] psW;
    return res;
}

void getSamplesMVN(double* mu, double** covM,  double ** sample, int D, int n)
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    Eigen::internal::scalar_normal_dist_op<double> randN; // Gaussian functor
    Eigen::internal::scalar_normal_dist_op<double>::rng.seed(seed); // Seed the rng
    Eigen::VectorXd mean(D);
    Eigen::MatrixXd covar(D,D);
    for(int c=0;c<D;c++)
    {
        mean(c)=mu[c];
        for(int d=0;d<D;d++)
            covar(c,d) = covM[c][d];
    }
    Eigen::MatrixXd normTransform(D,D);
    Eigen::LLT<Eigen::MatrixXd> cholSolver(covar);

    if (cholSolver.info()==Eigen::Success)
        normTransform = cholSolver.matrixL();
    else
    {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar);
        normTransform = eigenSolver.eigenvectors() * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
    }

    Eigen::MatrixXd samples = (normTransform * Eigen::MatrixXd::NullaryExpr(D,n,randN)).colwise()+ mean;
    for(int c=0;c<D;c++)
    {
        for(int d=0;d<n;d++)
        {
            if(samples(c,d)==-NAN || samples(c,d)==NAN)
                sample[c][d]=0.0;
            else
                sample[d][c]=samples(c,d);
        }
    }
}
