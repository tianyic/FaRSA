//
//  Parms.h
//  fPlusl1cpp
//
//  Created by chentianyi on 1/14/16.
//  Copyright (c) 2016 ChenTianyi. All rights reserved.
//

#ifndef fPlusl1cpp_Parms_h
#define fPlusl1cpp_Parms_h

#include "typeDef.h"

void setParms(struct Parameter &parms){

    parms.max_iter = 2;
    parms.max_time = 3600;
    parms.printlevel = 2;
    parms.printfile = 'testing';
    parms.printevery = 20;
    parms.Gamma = 1.0;
    parms.eta_r = 0.1;
    parms.eta = 0.01;
    parms.xi = 0.5;
    parms.tol_absolute = 1.0e-6;
    parms.tol_relative = 1.0e-6;
    parms.betaFrac = 0.1;
    parms.phiFrac = 0.2;
    parms.maxCG_iter = INFINITY;
    parms.maxCD_iter = 1000;
    parms.CrossOverTol = 0.1;
    parms.maxback = 100;
    parms.fracViol = 0.1;
    parms.nVmaxAllow = 1000;
    parms.checkResEvery = 1;
    parms.tryCG = true;
    parms.tryCD = false;
    parms.tryCrossOver = false;
    parms.termType = 1;
}

#endif
