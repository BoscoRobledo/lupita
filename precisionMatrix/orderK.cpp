





    /*int nSamples=500000;
    double** sampchl=new double*[nSamples];         ///Memory for the covariance matrix of the parent sets
    double** sampcht=new double*[nSamples];
    for(int d=0;d<nSamples;d++)
    {
        sampchl[d]=new double[D];
        sampcht[d]=new double[D];
    }
    double* mu=new double[D];
    getSamplesChowLiu(mu,covM,corrM,sampchl,D,nSamples);
    FILE* outchl;
    outchl=fopen("resCHL","w");
    for(int c=0;c<nSamples;c++)
    {
        for(int d=0;d<(D-1);d++)
            fprintf(outchl,"%.10lf\t",sampchl[c][d]);
        fprintf(outchl,"%.10lf\n",sampchl[c][D-1]);
    }
    fclose(outchl);

    getSamplesTChJT(mu,covM,corrM,sampcht,D,nSamples,2);
    FILE* outcht;
    outcht=fopen("resCHT","w");
    for(int c=0;c<nSamples;c++)
    {
        for(int d=0;d<(D-1);d++)
            fprintf(outcht,"%.10lf\t",sampcht[c][d]);
        fprintf(outcht,"%.10lf\n",sampcht[c][D-1]);
    }
    fclose(outcht);




    ch3.increaseOrder();
    ch3.increaseOrder();

    cout<<endl;

    cout<<endl;

    int* p= new int[D];
    TCHLib ch3(corrM,covM,D,2);

    for(int c=0;c<D;c++)
    {
        p[c]=ch3.getPerm()[c];
        cout<<p[c]<<" ";
    }

    cout<<endl<<ch3;

    getWfromTChJT(ch3,W,p);
    for(int c=0;c<D;c++)
        perm[c]=(double)p[c];
    delete []p;*/
    //TCHLib tl(CHLT,0,cM,cV);
    //tl.increaseOrder();
    //cout<<tl;

    /*cout<<endl<<"Beta:    ";
    for(int c=0;c<D;c++)
        cout<<setprecision(3)<<setw(7)<<beta[c]<<" ";
    cout<<endl<<"CondCov: ";
    for(int c=0;c<D;c++)
        cout<<setprecision(3)<<setw(7)<<cv[c]<<" ";*/
