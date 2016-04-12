//
//  SubspaceSolvers.h
//  fPlusl1cpp
//
//  Created by chentianyi on 1/16/16.
//  Copyright (c) 2016 ChenTianyi. All rights reserved.
//

#ifndef fPlusl1cpp_SubspaceSolvers_h
#define fPlusl1cpp_SubspaceSolvers_h

#include "typeDef.h"
#include <iostream>
#include "Utils.h"
#include <map>

typedef std::map<int, int> IdxMap;
template <class T> static inline void swap(T& x, T& y) { T t=x; x=y; y=t; }

struct OutputCG CGsolver(struct Parameter parms, double *X, LogRegCost fun, int *Snn, int *Sp, int *Sn, int *Snn_p,int nS, double *gradF_S){
    
    struct OutputCG outputCG;
    double eta_r = parms.eta_r;
    double fracViol = parms.fracViol;
    int nVmaxAllow = parms.nVmaxAllow;
    int maxCG_iter = parms.maxCG_iter;
    
    int n = fun.num_features;
    
    int nVmax = std::max(nVmaxAllow, int(fracViol*nS));
    int counter = 1;
    bool term = false;
    double *r0 = Malloc(double, nS);
    
    int i;
//    for (i=0; i< nS; i++) {
//        std::cout<<gradF_S[i]<<" "<<Snn[i]<<" "<<Snn_p[i]<<std::endl;
//    }
    
    Copy(gradF_S, r0, 1.0, nS);
    double normr0 = Norm(r0, nS);
    double normr = normr0;
    double *r = Malloc(double, nS);
    Copy(r0, r, 1.0, nS);
    double *p = Malloc(double, nS);
    Copy(r0, p, -1.0, nS);
    
    int max_loop = std::min(maxCG_iter, nS);
    //std::cout<<max_loop<<std::endl;
    double res_target = std::max(std::min(eta_r,pow(normr0, 1.5))*normr0,1e-10);
    //std::cout<<"res_target:"<<res_target<<std::endl;
    //std::cout<<nS<<std::endl;
    double *Hp = Malloc(double, nS);
    double pTHp;
    double alphaCG;
    double betaCG;
    double *d = Malloc(double, nS);
    double normd;
    double normr_old;
    double res;
    int nV;
    
    while (1) {
        
        fun.HvS(p, Hp, nS, Snn, Snn_p);
//        for (i = 0; i<nS; i++) {
//            std::cout<<Hp[i]<<" "<<p[i]<<std::endl;
//        }
        pTHp = dotProduct(p, Hp,nS);
        alphaCG = normr*normr/pTHp;
        //std::cout<<"alphaCG:"<<alphaCG<<std::endl;
        axPlusv(alphaCG, p, d, nS);
        normd = Norm(d, nS);
        axPlusv(alphaCG, Hp, r, nS);
        normr_old = normr;
        normr = Norm(r, nS);
        
        // Violating sets
        nV = 0;
        for (i = 0; i < nS; i++) {
            if (X[Snn_p[i]]+d[i]<0 && X[Snn_p[i]]>0) {
                nV++;
            }else if (X[Snn_p[i]]+d[i]>0 && X[Snn_p[i]]<0){
                nV++;
            }
            //std::cout<<Snn_p[i]<<std::endl;
        }
        
        // Check for termination of CG
        if (normr <= res_target) {
            //std::cout<<"CGtol"<<std::endl;
            term = true;
        }else if (nV > nVmax){
            term = true;
            //std::cout<<"CGvio"<<std::endl;
        }else if (normd >= 1e3){
            //std::cout<<"CGbig"<<std::endl;
            term = true;
        }else if (counter > max_loop){
            //std::cout<<"CGmax"<<std::endl;
            term = true;
        }
        
        if (term) {
            res = normr;
            //std::cout<<"res:"<<res<<std::endl;
            break;
        }
        
        betaCG = normr*normr/(normr_old*normr_old);
        axPlusby(-1.0, r, betaCG, p, nS);
//        for(int i=0;i<nS;i++){
//            std::cout<<p[i]<<std::endl;
//        }
//        std::cout<<"========="<<std::endl;
//        for(int i=0;i<nS;i++){
//            std::cout<<d[i]<<std::endl;
//        }
//        std::cout<<"========="<<std::endl;
        counter++;

    }
    //set outputCG
    outputCG.d = d;
    
//    for(int i=0;i<nS;i++){
//        std::cout<<d[i]<<std::endl;
//    }
//
    return(outputCG);
}

struct OutputCD CDsolver(struct Parameter parms, double *X, LogRegCost fun, int *Snn, int *Sp, int *Sn, int *Snn_p,int nS, double *gradF_S){
    
    struct OutputCD outputCD;
    double fracViol = parms.fracViol;
    int nVmaxAllow = parms.nVmaxAllow;
    double eta_r = parms.eta_r;
    int n = fun.num_features;
    int m = fun.num_samples;
    int maxCD_iter = parms.maxCD_iter;
    int *index = Malloc(int, nS);
    Copy_int(Snn_p, index, 1, nS);
    
    int nVmax = std::max(nVmaxAllow, int(fracViol*nS));
    
    int i, j, idx, ind;
    double *d = Malloc(double, nS);
    int  iter_CD = 0;
    bool term = false;
    
    double *r0 = Malloc(double, nS);
    Copy(gradF_S, r0, 1.0, nS);
    double normr0 = Norm(r0, nS);
    double res_target = std::max(std::min(eta_r,pow(normr0, 1.5))*normr0, 1e-10);
    
    double *Hdiag = Malloc(double, n);
    Feature_node *x;
    double *xTd = Malloc(double, n);
    int maxInnerLoop_Iter =  nS;
    int iter_innerloop = 0;
    double Hjj;
    double Gradj;
    double z;
    double normd;
    int nV;
    double decent;
    double sufficient_decent;
    IdxMap idxMap;
    
    // set Hdiag
    for (i = 0; i < nS; i++) {
        idx = Snn_p[i];
        Hdiag[idx] = 1e-8;
        idxMap[idx] = i;
        double tmp = 0.0;
        x = fun.X[idx];
        while (x->index != -1) {
            ind = x->index-1;
            Hdiag[idx] += x->value*x->value*fun.diag[ind];
            
        }
        Hdiag[idx] = Hdiag[idx]/(double)fun.num_samples;
    }
    

    
    
    while (true) {
        
        
        for (i = 0; i < m; i++) {
            xTd[i] = 0.0;
        }
        
        while (true) {
            if (iter_innerloop >= maxInnerLoop_Iter) {
                break;
            }
            // shuffle indexes
            for (i = 0; i < nS; i++) {
                j = i + rand()&(nS - i);
                swap(index[i], index[j]);
            }
            
            for (i = 0; i < nS; i++) {
                j = index[i];
                Hjj = Hdiag[j];
                x = fun.X[j];
                Gradj = gradF_S[j];
                while (x->index != -1) {
                    ind = x->index-1;
                    Gradj += x->value*fun.diag[ind]*xTd[ind];
                    x++;
                }
                
                z = -Gradj/Hjj;
                
                if (fabs(z) < 1e-12) {
                    continue;
                }
                
                // insert set d[]+= d[]+z
                // p code, find j in [1:nS]
                d[idxMap.at(j)] += z;
                z = std::min(std::max(z , -10.0), 10.0);
                SparseOperator::axPlusv(z, x, xTd);
                
            }
            iter_innerloop++;
        }
        
        normd = Norm(d, nS);
        nV = 0;
            
        for (i = 0; i < nS; i++) {
            if (X[Snn_p[i]]+d[Snn_p[i]] < 0 && X[Snn_p[i]] > 0) {
                nV ++;
            }else if (X[Snn_p[i]]+d[Snn_p[i]] > 0 && X[Snn_p[i]] < 0){
                nV ++;
            }
        }
            
        decent = dotProduct(gradF_S, d, nS);
        
        if (decent <= sufficient_decent && nV > nVmax) {
            term = true;
        }else if (decent <= sufficient_decent && normd >= 1e3){
            term = true;
        }else if (iter_CD >= maxCD_iter){
            term = true;
        }
        if (term) {
            break;
        }
        
        iter_CD++;
        
    }
    outputCD.d = d;
    outputCD.DirDec = decent;
    
    return(outputCD);
}


#endif
