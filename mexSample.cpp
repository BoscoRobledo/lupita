#include "mex.h"
#include <iostream>
#include "Simulate.h"

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    int D; //dimension of the problem
    int order;   //order of statistics used in the tree probability approx (1. All variables are independent, 2. Chow & Liu Tree, 3. 3-t-CherryJunctionTree )
    int n;  //number of samples to generate
    double *corrM;               // Population Correlation Matrix
    double *covM;               // Population Covariance Matrix
    double *mu;                 //population mean
    double **bdcorrM=(double**)malloc(D*sizeof(double*));  //pointer arrays to handle matrices as bidimensional arrays
    double **bdcovM=(double**)malloc(D*sizeof(double*));
    double **bdsamples=(double**)malloc(n*sizeof(double*));

    D = (int)mxGetScalar(prhs[0]);
    n = (int)mxGetScalar(prhs[1]);
    order = (int)mxGetScalar(prhs[2]);
    mu = mxGetPr(prhs[3]);
    covM = mxGetPr(prhs[4]);
    corrM = mxGetPr(prhs[5]);



    plhs[0] = mxCreateDoubleMatrix(n,D,mxREAL); //memory for precision matrix
    for(int c=0;c<D;c++)
    {
        bdcorrM[c]=corrM+c*D;
        bdcovM[c]=covM+c*D;
    }
    for(int c=0;c<n;c++)
        bdsamples[c]=mxGetPr(plhs[0])+c*D;

    if(order==1)
    {
        double* var=(double*)malloc(D*sizeof(double));
        for(int c=0;c<D;c++)
            var[c]=bdcovM[c][c];
        getSamplesMarginal(mu,var,bdsamples,D,n);
        free(var);
    }
    else if(order==-2)
        getSamplesChowLiu(mu,bdcovM,bdcorrM,bdsamples,D,n);
    else if(order==D)
        getSamplesMVN(mu,bdcovM,bdsamples,D,n);
    else
        getSamplesTChJT(mu,bdcovM,bdcorrM,bdsamples,D,n,order);
}




/* check for proper number of arguments
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Problem dimension, approx order and correlation matrix are required.");
    }
    if( !mxIsDouble(prhs[0]) ||
         mxIsComplex(prhs[0]) ||
         mxGetNumberOfElements(prhs[0])!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","First argument must be problem dimension.");
    }
    if( !mxIsDouble(prhs[1]) ||
         mxIsComplex(prhs[1]) ||
         mxGetNumberOfElements(prhs[1])!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","Second argument must be order approx.");
    }

    if( !mxIsDouble(prhs[2]) ||
         mxIsComplex(prhs[2])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Third argument must be population correlation matrix");
    }


    D = (int)mxGetScalar(prhs[0]);
    order = (int)mxGetScalar(prhs[1]);
    covM = mxGetPr(prhs[2]);
    corrM = mxGetPr(prhs[3]);
    if(nlhs!=order*2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","Returns the precision matrix for each order<=second argument order.");
    }
    Wfun* wf=new Wfun[order];
    wf[0]=DW;
    wf[1]=WChL;
    for(int c=0;c<order;c++)
    {
        plhs[c*2] = mxCreateDoubleMatrix(D,D,mxREAL); //memory for precision matrix
        plhs[c*2+1] = mxCreateDoubleMatrix(1,D,mxREAL); //memory for permutation indices
        W[c] = mxGetPr(plhs[c*2]);
        perm[c] = mxGetPr(plhs[c*2+1]);
        wf[c](W[c],covM,corrM,perm[c],D);

    }*/
