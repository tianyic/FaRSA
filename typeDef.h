//
//  typeDef.h
//  fPlusl1cpp
//
//  Created by chentianyi on 1/12/16.
//  Copyright (c) 2016 ChenTianyi. All rights reserved.
//

#ifndef fPlusl1cpp_typeDef_h
#define fPlusl1cpp_typeDef_h

struct problem{
    int m, n; // m is the number of samples, n is the number of features
    double *y;
    struct Feature_node **x;
    double bias;
};

// sparse data (index, value)
struct Feature_node{
    int index;
    double value;
};


struct Parameter{
    
    int max_iter;
    int max_time;
    int maxCG_iter;
    int maxCD_iter;
    bool tryCG;
    bool tryCD;
    bool tryCrossOver;
    int printlevel;
    std::string printfile;
    int printevery;
    int checkResEvery;
    double Gamma;
    double eta_r;
    double eta;
    double xi;
    double tol_absolute;
    double tol_relative;
    double betaFrac;
    double phiFrac;
    double CrossOverTol;
    double maxback;
    double fracViol;
    double nVmaxAllow;
    int termType;

    
};

struct OutputCG{
    double *d;
    std::string subprobFlag;
    double res;
    int subits;
    double HvProds;
    double res_target;
    int nV;
    int nVmax;
};

#endif
