//
//  SparseOperator.h
//  fPlusl1cpp
//
//  Created by chentianyi on 1/13/16.
//  Copyright (c) 2016 ChenTianyi. All rights reserved.
//

#ifndef fPlusl1cpp_SparseOperator_h
#define fPlusl1cpp_SparseOperator_h

#include "typeDef.h"

class SparseOperator{

public:
    
    // dot product
    static double dot(const double *v, const Feature_node *x){
        double vTx=0.0;
        while (x->index != -1) {
            //std::cout<<x->index<<":"<<x->value<<std::endl;
            vTx += v[x->index-1]*x->value;
            x++;
        }
        return(vTx);
    }
    
    // v=ax+bv
    static void axPlusbv(const double a, const double b, const Feature_node *x, double *v){
        while (x->index != -1) {
            v[x->index-1] = b*v[x->index-1] + a*x->value;
            x++;
        }
    }
    
    static void axPlusv(const double a, const Feature_node *x, double *v)
    {
        while(x->index != -1)
        {
            v[x->index-1]+=a*x->value;
            //std::cout<<x->index-1<<" "<<a<<" "<<x->value<<" "<<v[x->index-1]<<std::endl;
            x++;
        }
        //std::cout<<v[x->index-2]<<std::endl;
    }
    
    static void axPlusvS(const double a, const Feature_node *x, double *v, const int *Sn, const int *Sn_p, int nS){
        int i = 0;
        int pointer = 0;
        while (x->index != -1) {
            for ( i = pointer; i < nS; i++) {
                if (Sn_p[i] == x->index-1) {
                    v[i] += a*x->value;
                    pointer++;
                    break;
                }else if(Sn_p[i] < x->index-1){
                    pointer++;
                    continue;
                }else{
                    break;
                }
            }
            x++;
        }
    }
    
    static double dot1(const double *s, const Feature_node *x, const int *Sn, const int *Sn_p,int nS){
        
        int i;
        double value = 0.0;
        int pointer = 0;
        while (x->index != -1) {
            for (i = pointer; i < nS; i++) {
                if (Sn_p[i] == x->index-1) {
                    value += x->value*s[i];
                    pointer++;
                    break;
                }else if (Sn_p[i]< x->index-1){
                    pointer++;
                    continue;
                }else{
                    break;
                }
            }
            x++;
        }
        
        return(value);
    }
};

#endif
