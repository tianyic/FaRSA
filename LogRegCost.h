//
//  logCost.h
//  fPlusl1cpp
//
//  Created by chentianyi on 1/13/16.
//  Copyright (c) 2016 ChenTianyi. All rights reserved.
//

#ifndef fPlusl1cpp_logCost_h
#define fPlusl1cpp_logCost_h

#include "SparseOperator.h"
#include <math.h>
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

// setExpterm is contained in func, setSigmoid, and setDiag is contained in grad to decrease loop
class LogRegCost{
    
public:
    int num_features;       //number of features
    int num_samples;        //number of sample
    double *y;              //label
    struct Feature_node **X; // data matrix
    double bias=0.0;
    double lambda;
    double *expTerm;
    
    LogRegCost(){}
    
    
    double func(){
        
        double fvalue = 0.0;
        double *y = this->y;
        int m = this->num_samples;
        int n = this->num_features;
        
        
        // calculate logistic cost
        for (int i=0; i < m; i++) {
            
            double yexpTerm = y[i]*this->expTerm[i];
            
            if(yexpTerm >= 0){
                fvalue +=  log(1 + exp(-yexpTerm));
            }else{
                fvalue +=  (-yexpTerm + log(1+exp(yexpTerm)));
            }
        }
        
        return(fvalue/this->num_samples);
        
    }
    
    
    
    void grad(double *w, double *g){
        
        double *y = this->y;
        int m = this->num_samples;
        int n = this->num_features;
        
        this->sigmoid = Malloc(double, m);
        this->diag = Malloc(double, m);
        double *tmp = new double[m];

        for (int i = 0; i < m; i++) {
            sigmoid[i] = 1/(1+exp(-y[i]*expTerm[i]));
            //std::cout<<i<<":"<<sigmoid[i]<<std::endl;
            diag[i] = sigmoid[i]*(1-sigmoid[i]);
            tmp[i] = (sigmoid[i]-1)*y[i]/(double)this->num_samples;
        }
        XTv(tmp, g);
        
        delete [] tmp;
        
    }

    void Hv(double *s,double *Hs){
        
        int m = this->num_samples;
        int n = this->num_features;
        double *tmp = new double[m];
        Feature_node **X = this->X;
        
        for (int i = 0; i < n; i++) {
            Hs[i] = 0;
        }
        
        for (int i = 0; i < m; i++) {
            Feature_node * const Xi=X[i];
            tmp[i] = SparseOperator::dot(s, Xi);
            tmp[i] = this->diag[i]*tmp[i];
            SparseOperator::axPlusbv(tmp[i], 1.0, Xi, Hs);
        }
        delete [] tmp;
    }
    
    void HvS(double *s, double *Hs, int nS, int *Snn, int *Snn_p){
        
        int m = this->num_samples;
        double num_sam = (double) m;
        double *tmp = Malloc(double, this->num_samples);
        Feature_node **X = this->X;
        int i,j;
        
//        for (i = 0 ; i<nS; i++) {
//            std::cout<<Snn_p[i]<<":"<<s[i]<<std::endl;
//        }
        
        Feature_node *x = X[0];
        
        for (i = 0; i < nS; i++) {
            Hs[i] = 0.0;
        }
        
        for (i = 0; i < m; i++) {
            Feature_node * const Xi = X[i];
            tmp[i] = SparseOperator::dot1(s, Xi, Snn, Snn_p, nS);
            //std::cout<<tmp[i]<<std::endl;
            tmp[i] =  this->diag[i]*tmp[i]/num_sam;
            // need to modify
            //SparseOperator::axPlusv(tmp[i], Xi, Hs);
            SparseOperator::axPlusvS(tmp[i], Xi, Hs, Snn, Snn_p, nS);
        }
        delete [] tmp;
        
    }
    
    void Xv(double *v, double *Xv){
        
        int m = this->num_samples;
        Feature_node **X = this->X;
        for(int i=0;i < m; i++){
            Xv[i] = SparseOperator::dot(v, X[i]);
        }
        
    }
    
    void XTv(double *v,double *XTv){
        
        int m = this->num_samples;
        int n = this->num_features;
        Feature_node **X = this->X;
        
        for (int i = 0; i < n; i++) {
            XTv[i] = 0.0;
        }
        
        for (int i = 0; i < m; i++) {
            //std::cout<<v[i]<<std::endl;
            SparseOperator::axPlusv(v[i], X[i], XTv);
        }
        
    }
    
    void setExpterm(double *w){
        // set expTerm
        this->expTerm = Malloc(double, this->num_samples);
        Xv(w, this->expTerm);
    }
    
private:
    
    double *g;
    double *sigmoid;
    double *diag;
    double *tmp;
    
};


#endif
