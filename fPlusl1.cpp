//
//  main.cpp
//  fPlusl1cpp
//
//  Created by chentianyi on 1/12/16.
//  Copyright (c) 2016 ChenTianyi. All rights reserved.
//

#include <iostream>
#include <vector>
#include <math.h>
#include "typeDef.h"
#include "DataStream.h"
#include "LogRegCost.h"
#include "Parms.h"
#include "SubspaceSolvers.h"
#include "Utils.h"

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
#define INF HUGE_VAL


struct Parameter parms;
//template <class T> static inline void swap(T& x, T& y) { T t=x; x=y; y=t; }
template <class T> static inline void pointswap(T* x, T* y){T t=*x; *x=*y;*y=t;}



int main(int argc, const char * argv[]) {
    // insert code here...
    char input_file_name[1024];
    strcpy(input_file_name, "/Users/chentianyi/Documents/optimization research/numeric/fPlusl1cpp/fPlusl1cpp/RawData/real-sim");
    
    void print(double *name[],int n);
    void setBetaPhi(double *gradf,double *betaphi, double *absbetaphi, int *Sz, int *Sp, int *Sn, int *Snn, int n, double *X, double lambda,double &normBeta, double &normPhi, int &numz, int &nump, int &numn,int &numnn);
    void selectBeta(double frac, double *betaphi, double *absbetaphi,int *Sz, int numz, int numzz,int &nS);
    void selectPhi(double frac, double *betaphi, double *absbetaphi, int *Snn, int numnn, int &nS);
    double oneNorm(double *x, int n);
    void project(double *X, double *d, int *Snn_p, int nS, bool &isSameOrthant,double alpha);
    LogRegCost fun;
    
    // load file
    readProblem(input_file_name,fun);
    
    std::cout<<"num of samples:"<<fun.num_samples<<std::endl;
    std::cout<<"num of features:"<<fun.num_features<<std::endl;
    
    /****** set Parameters ******/
    setParms(parms);
    
    int m = fun.num_samples;
    int n = fun.num_features;
    int max_iter = parms.max_iter;
    int max_time = parms.max_time;
    int printlevel = parms.printlevel;
    std::string printfile = parms.printfile;
    int printevery = parms.printevery;
    double Gamma = parms.Gamma;
    double eta_r = parms.eta_r;
    double eta = parms.eta;
    double xi = parms.xi;
    double tol_absolute = parms.tol_absolute;
    double tol_relative = parms.tol_relative;
    double betaFrac = parms.betaFrac;
    double phiFrac = parms.phiFrac;
    int maxCG_iter = parms.maxCG_iter;
    int maxCD_iter = parms.maxCD_iter;
    double CrossOverTol = parms.CrossOverTol;
    double maxback = parms.maxback;
    double fracViol = parms.fracViol;
    int nVmaxAllow = parms.nVmaxAllow;
    int checkResEvery = parms.checkResEvery;
    bool tryCG = parms.tryCG;
    bool tryCD = parms.tryCD;
    bool tryCrossOver = parms.tryCrossOver;
    int termType = parms.termType;
    
    double lambda = 1.0/fun.num_samples;
    /***************************/
    int i;
    // convert labels into 1, -1, since some RawDatasets don't use 1, -1 as their default labels
    double *y =  Malloc(double, m);
    for ( i = 0; i < m; i++) {
        if (fun.y[i] > 0) {
            y[i] = 1.0;
        }else{
            y[i] = -1.0;
        }
        //std::cout<<y[i]<<std::endl;
    }
    fun.y = y;
    
    
    
    // Initialize X
    double *X = Malloc(double, n);
    for (i = 0; i < n; i++) {
        X[i] = 0.0;
    }
    
    fun.setExpterm(X);
    //std::cout<<fun.func(X)<<std::endl;
    double *gradf = Malloc(double, n); // gradient of f
    fun.grad(X, gradf);
    int iter = 0;

    double normBeta; // calculate 1-norm
    double normPhi;  // calculate 1-norm
    double normBeta0;
    double normPhi0;
    double norm0;
    
    int numz;
    int nump;
    int numn;
    int numnn;

    int j;
    double alpha;
    double F_old;
    double F;
    
    std::cout<<"lambda:"<<lambda<<std::endl;

//    for (int i = 0; i < n; i++) {
//        std::cout<<i<<":"<<gradf[i]<<std::endl;
//    }
    
    clock_t startTime = clock();
    
    while (1) {
        
        double *betaphi= Malloc(double, n);  // a "mixture" of beta, and phi
        double *absbetaphi = Malloc(double, n); // absolute betaphi
        int *Sz = Malloc(int , n);
        int *Sn = Malloc(int , n);
        int *Sp = Malloc(int , n);
        int *Snn = Malloc(int, n);
        struct OutputCG outputCG;
        struct OutputCD outputCD;

        
//        for(i = 0;i<n;i++){
//            std::cout<<X[i]<<std::endl;
//        }
        
        setBetaPhi(gradf, betaphi,absbetaphi, Sz, Sp, Sn, Snn, n, X, lambda, normBeta, normPhi, numz, nump, numn, numnn);
        if (iter == 0) {
            normBeta0 = normBeta;
            normPhi0 = normPhi;
            norm0 = std::max(normBeta0, normPhi0);
            std::cout<<"norm0:"<<norm0<<std::endl;
        }
        
        //std::cout<<numz<<std::endl;
        /*****/

        
        // Termination
        if (termType == 1) {
            //std::cout<<"error:"<<std::max(normBeta, normPhi)<<std::endl;
            if (std::max(normBeta, normPhi) <= std::max(tol_absolute, tol_relative*norm0)) {
                std::cout<<"Optimal solution has been founded."<<" normbeta:"<<normBeta<<" normPhi:"<<normPhi<<std::endl;
                break;
            }
        }
        
        iter++;
        if (iter > max_iter) {
            std::cout<<"Maximum iteration has been reached"<<std::endl;
            break;
        }
//        std::cout<<" normbeta:"<<normBeta<<" normPhi:"<<normPhi<<std::endl;
//        for (int i = 0; i < n; i++) {
//            std::cout<<i<<":"<<betaphi[i]<<" "<<X[i]<<std::endl;
//        }
        
        if (normBeta <= Gamma*normPhi) {
            //std::cout<<"PHIPHI"<<std::endl;
            int nS = 0;
            int i;
            
            // need to modify, Snn_p should be with length nS
            int *Snn_p = Malloc(int, numnn);
            Copy_int(Snn, Snn_p, 1, numnn);
            //std::cout<<"numnn:"<<numnn<<std::endl;
//            for (i=0; i< numnn; i++) {
//                std::cout<<Snn_p[i]<<std::endl;
//            }
            
            selectPhi(phiFrac, betaphi, absbetaphi, Snn, numnn, nS);
            
            double *d = Malloc(double, nS);
            
            if (std::max(normBeta, normPhi) <= CrossOverTol && tryCrossOver) {
                tryCG = true;
                tryCD = false;
            }
            
//            for (i = 0; i<n; i++) {
//                std::cout<<i<<" "<<betaphi[i]<<" "<<absbetaphi[i]<<std::endl;
//            }

            
            double *gradF_S = Malloc(double, nS);

            for (i = 0 ; i < nS; i++) {
                if (X[Snn_p[i]] > 0) {
                    gradF_S[i] = gradf[Snn_p[i]] + lambda;
                }else if (X[Snn_p[i]] < 0){
                    gradF_S[i] = gradf[Snn_p[i]] - lambda;
                }
                //std::cout<<gradF_S[i]<<" "<<" "<<gradf[Snn_p[i]]<<" "<<Snn[i]<<" "<<Snn_p[i]<<" "<<lambda<<"  gradF_S"<<std::endl;
            }
            
            double DirDer = 0.0;
            if (tryCG) {
                outputCG = CGsolver(parms, X, fun, Snn, Sp, Sn,Snn_p, nS, gradF_S);
                d = outputCG.d;
                for (i = 0; i < nS; i++) {
                    DirDer += gradF_S[i]*d[i];
                }
            }
            
            if (tryCD) {
                outputCD = CDsolver(parms, X, fun, Snn, Sp, Sn, Snn_p, nS, gradF_S);
                d = outputCD.d;
                DirDer = outputCD.DirDec;
            }
            
            
            // line search
            j = 0;
            alpha = 1.0;
            F_old =fun.func()+lambda*oneNorm(X, n);

            
            double *tmp = Malloc(double,n);
            Copy(X,tmp,1.0,n);
            
            bool isSameOrthant = true;
            //std::cout<<"DirDer:"<<DirDer<<std::endl;
            
//            for (i=0; i<nS; i++) {
//                std::cout<<Snn_p[i]<<" "<<Snn[i]<<std::endl;
//            }
            
            while (1) {
                
                // may need improvement
                project(X, d, Snn_p, nS, isSameOrthant,alpha);
                
                fun.setExpterm(X);
                F =fun.func()+lambda*oneNorm(X, n);

                if (F - F_old <= eta * alpha * DirDer && isSameOrthant) {
                    std::cout<<"iter:"<<iter-1<<" "<<F<<" Phi descent step"<<std::endl;
                    break;
                }else if (F < F_old && !isSameOrthant){
                    std::cout<<"iter:"<<iter-1<<" "<<F<<" Phi sign change step"<<std::endl;
                    break;
                }
                
                if (j > maxback) {
                    break;
                }
                j++;
                alpha *= xi;
                
                for (i = 0; i < nS; i++) {
                    X[Snn_p[i]] = tmp[Snn_p[i]];
                }
                
            }
            //std::cout<<"j:"<<j<<" logf:"<<fun.func()<<std::endl;
        }else if (normBeta > Gamma*normPhi){
            
            int nS;
            int i;
            selectBeta(betaFrac, betaphi, absbetaphi, Sz, numz,n, nS);
            
            // scale Beta
            double scaleFactor = 0.1;
            double normCBeta = 0.0;
            for (i = 0; i < nS; i++) {
                normCBeta += betaphi[Sz[numz-i-1]]*betaphi[Sz[numz-i-1]];
                //std::cout<<normCBeta<<" ";
            }
            //std::cout<<std::endl;
            normCBeta = sqrt(normCBeta);
            //std::cout<<nS<<std::endl;
            //std::cout<<betaphi[Sz[numz-1]]<<std::endl;
            //std::cout<<"normCBeta"<<normCBeta<<std::endl;
            //std::cout<<scaleFactor<<std::endl;
            
            //std::cout<<"nS"<<nS<<" "<<numz<<std::endl;
            if (normCBeta > scaleFactor) {
                for (i = 0; i < nS; i++) {
                    betaphi[Sz[numz-i-1]] = -scaleFactor*betaphi[Sz[numz-i-1]]/normCBeta;
                }
            }else{
                for (i = 0; i < nS; i++) {
                    betaphi[Sz[numz-i-1]] = -betaphi[Sz[numz-i-1]];
                }
            }
            
//            for (i = 0; i < n; i++) {
//                std::cout<<"scalebetaphi"<<betaphi[i]<<std::endl;
//            }
            
            //std::cout<<betaphi[numz-1]<<std::endl;
            
            // set DirDer
            double DirDer = 0.0;
            for (i = 0; i < nS; i++) {
                if (betaphi[Sz[numz-i-1]]>0) {
                    DirDer += (gradf[Sz[numz-i-1]]+lambda)*betaphi[Sz[numz-i-1]];
                }else{
                    DirDer += (gradf[Sz[numz-i-1]]-lambda)*betaphi[Sz[numz-i-1]];
                }
            }
            //std::cout<<"DirDer:"<<DirDer<<std::endl;
            
            // line search
            j = 0; // line search num
            F_old =fun.func()+lambda*oneNorm(X, n);
            
            alpha = 1.0;

            double *tmp = Malloc(double,n);
            
            while (1) {

                // may need improvement
                for (i = 0; i < nS; i++) {
                    tmp[Sz[numz-i-1]] = X[Sz[numz-i-1]];
                    X[Sz[numz-i-1]] += betaphi[Sz[numz-i-1]]*alpha;
                    //std::cout<<Sz[numz-i-1]<<":"<<X[Sz[numz-i-1]]<<std::endl;
                }
                fun.setExpterm(X);
                F =fun.func()+lambda*oneNorm(X, n);
                //std::cout<<F<<std::endl;
                if (F - F_old <= eta * alpha * DirDer) {
                    std::cout<<"iter:"<<iter-1<<" "<<F<<" Beta step"<<std::endl;
                    break;
                }
                if (j >= maxback) {
                    break;
                }
                
                j++;
                alpha *= xi;
                for (i = 0; i < nS; i++) {
                    X[Sz[numz-i-1]] = tmp[Sz[numz-i-1]];
                }
                
            }
        }
        fun.grad(X, gradf);
        // need to reset betaphi, absbetaphi, Sz, Sp, Sn, Snn
        
        
    }
    
    std::cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds." << std::endl;
    
    return 0;
    
}

void setBetaPhi(double *gradf,double *betaphi, double *absbetaphi, int *Sz, int *Sp, int *Sn, int *Snn,int n, double *X, double lambda, double &normBeta, double &normPhi, int &numz, int &nump, int &numn, int &numnn){
    
    int i;
    double tmp;
    
    normBeta = 0.0;
    normPhi = 0.0;
    numz = 0;
    nump = 0;
    numn = 0;
    numnn = 0;
    
    for (i = 0; i < n; i++) {
        if (X[i] == 0) {
            if (gradf[i] < -lambda) {
                betaphi[i] = gradf[i] + lambda;
                Sz[numz] = i;
                numz++;
                normBeta += betaphi[i]*betaphi[i];
                absbetaphi[i] = -betaphi[i];
                continue;
            }else if (gradf[i] > lambda){
                betaphi[i] = gradf[i] - lambda;
                Sz[numz] = i;
                numz++;
                normBeta += betaphi[i]*betaphi[i];
                absbetaphi[i] = betaphi[i];
                continue;
            }
            
        }else if (X[i] > 0){

            betaphi[i] = gradf[i] + lambda;
            if(gradf[i]+lambda > 0){
                betaphi[i] = std::min(betaphi[i], std::max(X[i],gradf[i]-lambda));
            }
            absbetaphi[i] = betaphi[i]>=0?betaphi[i]:-betaphi[i];
            normPhi += betaphi[i]*betaphi[i];
            Sp[nump] = i;
            nump++;
            Snn[numnn] = i;
            numnn++;
            continue;
        }else if (X[i] < 0){
            betaphi[i] = gradf[i] - lambda;
            if (gradf[i]-lambda < 0) {
                betaphi[i] = std::max(betaphi[i],std::min(X[i],gradf[i]+lambda));
            }
            absbetaphi[i] = betaphi[i]>=0?betaphi[i]:-betaphi[i];
            normPhi += betaphi[i]*betaphi[i];
            Sn[numn] = i;
            numn++;
            Snn[numnn] = i;
            numnn++;
            continue;

        }
        

    }

    normBeta =  sqrt(normBeta);
    normPhi = sqrt(normPhi);
    //TO DO: new phi
}

void selectBeta(double frac, double *betaphi, double *absbetaphi,int *Sz, int numz,int n, int &nS){
    void bubbleSort(double *betaphi, double *absbetaphi, int num, int *S, int nS);
    if (frac < 1) {
        nS = (int)floor(n*frac);
        nS = std::max(nS, 1);
        nS = std::min(nS, numz);
        bubbleSort(betaphi, absbetaphi,numz, Sz, nS);
    }
}

void selectPhi(double frac, double *betaphi, double *absbetaphi, int *Snn, int numnn, int &nS){
    void bubbleSort(double *betaphi, double *absbetaphi, int num, int *S, int nS);
    if (frac < 1){
        nS = (int)floor(numnn*frac);
        nS = std::max(nS,1000);
        nS = std::min(nS, numnn);
        bubbleSort(betaphi, absbetaphi, numnn, Snn, nS);
        
    }
}


void bubbleSort(double *betaphi, double *absbetaphi, int num, int *S, int nS){
    
    
    int i , j;
    int *tmpS = Malloc(int, num);
    for (i = 0; i < num; i++) {
        tmpS[i] = S[i];
    }
    //std::cout<<numz<<std::endl;
    for (i = 0; i < nS; i++) {
        for (j = 0; j < num-1-i; j++) {
            if (absbetaphi[tmpS[j]] > absbetaphi[tmpS[j+1]]) {
                swap(absbetaphi[tmpS[j]],absbetaphi[tmpS[j+1]]);
                //swap(betaphi[tmpS[j]],betaphi[tmpS[j+1]]);
                swap(S[j], S[j+1]);
            }
        }
    }
    
}

void project(double *X, double *d, int *Snn_p, int nS, bool &isSameOrthant,double alpha){
    int i;
    for (i = 0; i < nS; i++) {
        if (X[Snn_p[i]] > 0) {
            X[Snn_p[i]] += alpha*d[i];
            if (X[Snn_p[i]] <= 0) {
                X[Snn_p[i]] = 0.0;
                isSameOrthant = false;
                continue;
            }
            continue;
        }else if(X[Snn_p[i]] < 0){
            X[Snn_p[i]] += alpha*d[i];
            if (X[Snn_p[i]] >= 0) {
                X[Snn_p[i]] = 0.0;
                isSameOrthant = false;
                continue;
            }
            continue;
        }
    }
}

void print(double *name,int n)
{
    int i;
    for(i=0;i<n;i++)
        std::cout<<name[i]<<std::endl;
}