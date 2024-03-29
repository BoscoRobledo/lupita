#include "mex.h"
#include <iostream>
#include "PrecMatrix.h"
/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    int D; //dimension of the problem
    int order;   //order of statistics used in the tree probability approx (1. All variables are independent, 2. Chow & Liu Tree, 3. 3-t-CherryJunctionTree )
    double *corrM;               /* Population Correlation Matrix */
    double *covM;               /* Population Correlation Matrix */
    double *W;               /* Population Correlation Matrix */
    double *perm;              /* Population permutation indexes*/
    double **bdcorrM=(double**)malloc(D*sizeof(double*));
    double **bdcovM=(double**)malloc(D*sizeof(double*));
    double **bdW=(double**)malloc(D*sizeof(double*));

    D = (int)mxGetScalar(prhs[0]);
    order = (int)mxGetScalar(prhs[1]);
    covM = mxGetPr(prhs[2]);
    corrM = mxGetPr(prhs[3]);
    plhs[0] = mxCreateDoubleMatrix(D,D,mxREAL); //memory for precision matrix
    plhs[1] = mxCreateDoubleMatrix(D,1,mxREAL); //memory for permutation indexes
    W = mxGetPr(plhs[0]);
    perm = mxGetPr(plhs[1]);
    /*mexPrintf("bdW Base pointer: %p\n",bdW);
    mexPrintf("W Base pointer: %p\n",W);
    mexPrintf("bdcovM Base pointer: %p\n",bdcovM);
    mexPrintf("covM Base pointer: %p\n",covM);
    mexPrintf("bdcorrM Base pointer: %p\n",bdcorrM);
    mexPrintf("corrM Base pointer: %p\n",corrM);*/
    for(int c=0;c<D;c++)
    {
        bdcorrM[c]=&corrM[c*D];  //Note: Matlab stores matrices in a column-wise manner. Being covariance, correlation and precision matrix symmetric, we can handle it in a row-wise manner.
        bdcovM[c]=&covM[c*D];
        bdW[c]=&W[c*D];
        /*mexPrintf("bdW pointer %d: %p\n",c,bdW[c]);
        mexPrintf("bdcovM pointer %d: %p\n",c,bdcovM[c]);
        mexPrintf("bdcorrM pointer %d: %p\n",c,bdcorrM[c]);*/
    }
    /*mexPrintf("bdW Base pointer: %p\n",bdW);
    mexPrintf("W Base pointer: %p\n",W);
    mexPrintf("bdcovM Base pointer: %p\n",bdcovM);
    mexPrintf("covM Base pointer: %p\n",covM);
    mexPrintf("bdcorrM Base pointer: %p\n",bdcorrM);
    mexPrintf("corrM Base pointer: %p\n",corrM);*/
    if(order==-2)
        WChL(bdW,bdcovM,bdcorrM,perm,D);
    else
        getW(bdW,bdcovM,bdcorrM,perm,D,order);
    /*mexPrintf("bdW Base pointer: %p\n",bdW);
    mexPrintf("W Base pointer: %p\n",W);
    mexPrintf("bdcovM Base pointer: %p\n",bdcovM);
    mexPrintf("covM Base pointer: %p\n",covM);
    mexPrintf("bdcorrM Base pointer: %p\n",bdcorrM);
    mexPrintf("corrM Base pointer: %p\n",corrM);*/
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
