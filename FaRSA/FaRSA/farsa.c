//
//  farsa.c
//  FaRSA
//
//  Created by chentianyi on 10/12/16.
//  Copyright (c) 2016 ChenTianyi. All rights reserved.
//
#include "farsa.h"


void farsa( int argc, char **argv, struct Input_FaRSA *input_farsa ){
    
    /* get arguments from parser */
    parse_command_line( argc, argv );
    
    /* load setting from profile */
    load_setting_profile();
    
    fprintf(stderr, "Problem file: %s\n", param.prob_name  );
    fprintf(stderr, "Format %s\n", param.data_format );
    /* objective function is not personalized, then read problem */
    
    fprintf(stderr, "param.loss_type %d\n", param.loss_type );
    if( param.loss_type != 3 ){
        
        read_problem();
        
        m = problem.num_samples;
        n = problem.num_features;
        y = Malloc( double, m );
        dm = (double) m;
        
    }else if( param.loss_type == 3 ){
        if( n == 0 ){
            fprintf(stderr, "Current n is 0. Please set n.\n" );
            exit(1);
        }
        m = 1;
        dm = (double) m;
    }else{
        fprintf(stderr, "Loss type is not successfully loaded\n" );
        exit(1);
    }
    
    fprintf(stderr, "n : %d, m : %d\n", n, m );
    
    int i, j;
    
    iter = 0;
    beta_iter = 0;
    phi_iter = 0;
    
    rsType = 1;
    norm_beta = 0.0;
    norm_phi = 0.0;
    norm_phi0 = 0.0;
    norm_beta0 = 0.0;
    
    beta = Malloc( double, n );
    phi = Malloc( double, n );
    
    grad_f = Malloc( double, n );
    grad_F = Malloc( double, n );
    
    S = Malloc( int, n + 1 );
    Sz = Malloc( int, n );
    Sn = Malloc( int, n );
    Sp = Malloc( int, n );
    x = Malloc( double, n );
    Gamma = param.Gamma;
    tol_absolute = param.tol_absolute;
    tol_relative = param.tol_relative;
    
    d = Malloc( double, n );
    TRradiusBeta = param.TRradiusBeta;
    TRradiusPhi = param.TRradiusPhi;
    alpha = 1.0;
    x_linesearch = Malloc( double, n );
    dirDer = 0.0;
    F_old = 0.0;  /* objective function value for last step */
    F = 0.0; /* objective function value on current step */
    step = Malloc( double, n ); /* difference between two steps for trust region update */
    norm_step = 0.0;
    maxback = param.maxback;
    xi = param.xi;
    sameOrthant = TRUE;
    
    lambda = param.lambdaCoeff / (double) problem.num_samples;
    
    fprintf(stderr, "111111\n" );
    /* if loss function has optimized version, then directly call them
     if not, run the generic routine. */
    clock_t begin = clock();
    if( param.loss_type == 1 ){
        fprintf(stdout, "Logistic loss plus l1 regularizer...\n" );
        // problem = *prob;
        logistic_loss( param );
    }else if( param.loss_type == 2 ){
        fprintf(stdout, "Lasso Loss plus l1 regularizer...\n" );
    }else if( param.loss_type == 3 ){
        fprintf(stdout, "Personalized loss function...\n" );
        // problem = *prob;
        generic_loss( param, input_farsa );
    }
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    fprintf( stdout, "runtime: %fs\n", time_spent );
    fprintf( stderr, "Problem: %s, runtime: %fs\n", param.prob_name, time_spent );
    
}

/****************** functions that can be used for any function including special optimized loss function ***/
struct projResult project( double *x, double *d, double alpha ){
    
    struct projResult projresult;
    double *x_proj_line;
    x_proj_line = Malloc( double, n );
    int i;
    int same_sign;
    same_sign = TRUE;
    
    for( i = 0; i < n; i++ ){
        x_proj_line[i] = x[i] + alpha * d[i];
        if( Sp[i] == 1 && x_proj_line[i] < 0 ){
            x_proj_line[i] = 0.0;
            same_sign = FALSE;
        }else if( Sn[i] == 1 && x_proj_line[i] > 0 ){
            x_proj_line[i] = 0.0;
            same_sign = FALSE;
        }
    }
    
    projresult.project_vector = x_proj_line;
    projresult.same_sign = same_sign;
    
    return projresult;
}

void setBetaPhi(){
    
    int i;
    int *p;
    p = S;
    
    norm_beta = 0.0;
    norm_phi = 0.0;
    beta = Malloc( double, n );
    phi = Malloc( double, n );
    
    for ( i = 0; i < n; i++ ){
        
        if( x[i] == 0 ){
            Sz[i] = 1;
            Sp[i] = 0;
            Sn[i] = 0;
            if ( grad_f[i] < -lambda ){
                beta[i] = grad_f[i] + lambda;
            }else if( grad_f[i] > lambda ){
                beta[i] = grad_f[i] - lambda;
            }
            // use l2 norm at first
            norm_beta += beta[i] * beta[i];
            phi[i] = 0.0;
        }else if( x[i] > 0 ){
            phi[i] = grad_f[i] + lambda;
            if( phi[i] > 0 ){
                phi[i] = MIN( phi[i], MAX( x[i], grad_f[i] - lambda ) );
            }
            norm_phi += phi[i] * phi[i];
            Sp[i] = 1;
            Sz[i] = 0;
            Sn[i] = 0;
            grad_F[i] = grad_f[i] + lambda;
            if( phi[i] != 0.0 ){
                nS++;
                *p = i;
                p++;
            }
            
        }else{
            phi[i] = grad_f[i] - lambda;
            if( phi[i] < 0 ){
                phi[i] = MAX( phi[i], MIN( x[i], grad_f[i] + lambda ) );
            }
            norm_phi += phi[i] * phi[i];
            Sn[i] = 1;
            Sz[i] = 0;
            Sp[i] = 0;
            grad_F[i] = grad_f[i] - lambda;
            if( phi[i] != 0.0 ){
                nS++;
                *p = i;
                p++;
            }
        }
        
    }
    *p = -1;
    p = S;
    
    norm_phi = sqrt( norm_phi );
    norm_beta = sqrt( norm_beta );
    
}

/*************************************************************/

/***********  Personalized Loss Section Start ****************/

void generic_loss(
                  const struct Parameter param,
                  struct Input_FaRSA *input_farsa
                  ){
    int i, j;
    int iter = 0;
    F = input_farsa->func(x);
    // fprintf(stderr, "%f\n", F);
    grad_f = input_farsa->grad_f( x );
    
    while(1){
        
        
        setBetaPhi();
        fprintf(stderr, "Iter: %d\n", iter );
        fprintf(stderr, "norm_beta: %f, norm_phi: %f\n", norm_beta, norm_phi );
        
        if( iter == 0 ){
            norm_phi0 = norm_phi;
            norm_beta0 = norm_beta;
            ttol = MAX( tol_absolute, tol_relative * MAX( norm_beta0, norm_phi0 ) );
        }
        
        // termination
        if( MAX( norm_beta, norm_phi ) <= ttol ){
            fprintf(stdout, "Optimal solution has been found.\n" );
            fprintf(stdout, "Objective function value: %f\n", F );
            fprintf(stdout, "Iteration: %d\n", iter);
            break;
        }
        
        // break;
        if( iter >= param.max_iter ){
            fprintf(stdout, "Maximum iteration has been reached.\n" );
            break;
        }
        
        if( norm_beta <= Gamma * norm_phi ){
            
            iter_type = 1;
            struct InputCG input_cg;
            input_cg.n = n;
            input_cg.nS = nS;
            input_cg.x = x;
            input_cg.S = S;
            input_cg.Sp = Sp;
            input_cg.Sn = Sn;
            input_cg.eta_r = param.eta_r;
            input_cg.grad_F = grad_F;
            input_cg.nVmaxAllow = ceil( MAX((double)param.nVmaxAllow,input_cg.fracViol*(double)nS) );
            input_cg.TRradius = TRradiusPhi;
            input_cg.maxCG_iter = param.max_CG_iter;
            input_cg.fracViol = param.fracViol;
            input_cg.rsType = rsType;
            input_cg.hessVecProd = input_farsa->hessVec;
            
            struct OutputCG output_cg = CGsolver( input_cg );
            
            d = output_cg.d_full;
            
            if( output_cg.nV == 0 && norm_beta < tol_absolute ){
                rsType = 2;
            }else{
                rsType = 1;
            }
            
            // Perform a line search
            j = 0;
            alpha = 1.0;
            dirDer = output_cg.dirDer;
            double suf_descent_1 = param.eta * dirDer;
            F_old = F;
            struct projResult projresult = project( x, d, alpha );
            x_linesearch = projresult.project_vector;
            sameOrthant = projresult.same_sign;
            F = input_farsa->func(x_linesearch);
            
            while(1){
                if( sameOrthant && F - F_old <=  suf_descent_1 * alpha ){
                    // set x, and grad_f
                    for( i = 0; i < n; i++ ){
                        step[i] = x_linesearch[i] - x[i];
                        norm_step += step[i] * step[i];
                        x[i] = x_linesearch[i];
                        // grad_f[i] = -1.0 / dm * dot_r( ysigmoid, problem.cX[i] );
                    }
                    norm_step = sqrt( norm_step );
                    break;
                }
                
                if( !sameOrthant && F < F_old ){
                    // set x, and grad_f
                    for( i = 0; i < n; i++ ){
                        step[i] = x_linesearch[i] - x[i];
                        norm_step += step[i] * step[i];
                        x[i] = x_linesearch[i];
                        // grad_f[i] = -1.0 / dm * dot_r( ysigmoid, problem.cX[i] );
                    }
                    norm_step = sqrt( norm_step );
                    break;
                }
                
                if ( j > maxback ){
                    // set x, and grad_f
                    for( i = 0; i < n; i++ ){
                        step[i] = x_linesearch[i] - x[i];
                        norm_step += step[i] * step[i];
                        x[i] = x_linesearch[i];
                        // grad_f[i] = -1.0 / dm * dot_r( ysigmoid, problem.cX[i] );
                    }
                    norm_step = sqrt( norm_step );
                    break;
                }
                
                alpha *= xi;
                projresult = project( x, d, alpha );
                x_linesearch = projresult.project_vector;
                sameOrthant = projresult.same_sign;
                F = input_farsa->func(x_linesearch);
                j++;
                
            }
            
            grad_f = input_farsa->grad_f( x_linesearch );
            
            
        }else{
            
            iter_type = 0;
            double TRbeta_norm_beta = -( TRradiusBeta / norm_beta );
            // fprintf(stderr, "%f\n", TRbeta_norm_beta );
            for ( i = 0; i < n; i++ ){
                d[i] = TRbeta_norm_beta * beta[i];
                x_linesearch[i] = x[i] + d[i];
                if( d[i] > 0 ){
                    grad_F[i] = grad_f[i] + lambda;
                }else if( d[i] < 0 ){
                    grad_F[i] = grad_f[i] - lambda;
                }else{
                    grad_F[i] = 0.0;
                }
            }
            
            double normd = TRradiusBeta;
            // Perform line search to get updated x.
            j = 0;
            alpha = 1.0;
            dirDer = dot_n5( grad_F, d, n );
            
            F_old = F;
            
            F = input_farsa->func(x_linesearch);
            
            double suf_descent_1 = param.eta * dirDer;
            double tmp = 0.0; // used for update grad_f
            
            // perform line search
            while( 1 ){
                if( F - F_old <= alpha * suf_descent_1 ){
                    // set x, and grad_f
                    for( i = 0; i < n; i++ ){
                        step[i] = x_linesearch[i] - x[i];
                        norm_step += step[i] * step[i];
                        x[i] = x_linesearch[i];
                        // grad_f[i] = -1.0 / dm * dot_r( ysigmoid, problem.cX[i] );
                    }
                    norm_step = sqrt( norm_step );
                    break;
                }
                if ( j > maxback ){
                    // set x, and grad_f
                    for( i = 0; i < n; i++ ){
                        step[i] = x_linesearch[i] - x[i];
                        norm_step += step[i] * step[i];
                        x[i] = x_linesearch[i];
                        // grad_f[i] = -1.0 / dm * dot_r( ysigmoid, problem.cX[i] );
                    }
                    norm_step = sqrt( norm_step );
                    break;
                }
                j++;
                alpha *= xi;
                
                // set new trial step
                for( i = 0 ; i < n; i++ ){
                    x_linesearch[i] = x[i] + alpha * d[i];
                }
                
                F = input_farsa->func(x_linesearch);
            }
            
            grad_f = input_farsa->grad_f( x_linesearch );
            
            // fprintf(stderr, "%f\n", F);
            
        }
        
        d = Malloc( double, n );
        nS = 0;
        // Set trust-region radius for next iteration.
        if( iter_type == 1 ){
            TRradiusPhi = MAX( 1E-3, MIN( 1E3, 10 * norm_step ) );
        }else{
            TRradiusBeta = MAX( 1E-5, MIN( 1.0, norm_step ) );
        }
        norm_step = 0.0;
        iter++;
    }
    
    
}


/***********  Personalized Loss Section End   ****************/

/***********  Logistic Regression Loss Section Start *********/

void logistic_loss( const struct Parameter param ){
    
    int i, j;
    
    double *sigmoid;
    double *ysigmoid;
    double *expterm;
    double *wTx;
    
    /* help_vector for logistic loss is sigmoid vector */
    help_vector = Malloc( double, m );
    ysigmoid = Malloc( double, m );
    expterm = Malloc( double, m );
    wTx = Malloc( double, m );
    
    fprintf(stderr, "Line 428 \n" );
    
    /* Set labels and Initialize some stuff */
    for ( i = 0; i < m; i++ ){
        if( problem.y[i] > 0 )
            y[i] = 1.0;
        else
            y[i] = -1.0;
        expterm[i] = 1.0;
        help_vector[i] = 0.5;
        ysigmoid[i] = y[i] * help_vector[i];
    }
    
    fprintf(stderr, "Line 441 \n" );
    
    // Continue initialize
    for( i = 0 ; i < n; i++ ){
        x[i] = 0.0;
        grad_f[i] = -1.0  / dm * dot_r( ysigmoid, problem.cX[i] );
    }
    
    F = logistic_func( x, expterm );
    
    fprintf(stderr, "Line 451 \n" );
    
    // fprintf(stderr, "%f\n", F);
    
    while(1){
        
        setBetaPhi();
        // fprintf(stderr, "norm beta:%f norm phi:%f, gamma: %f\n", norm_beta, norm_phi, Gamma );
        
        
        if( iter == 0 ){
            norm_phi0 = norm_phi;
            norm_beta0 = norm_beta;
            ttol = MAX( tol_absolute, tol_relative * MAX( norm_beta0, norm_phi0 ) );
        }
        
        // termination
        if( MAX( norm_beta, norm_phi ) <= ttol ){
            fprintf(stdout, "Optimal solution has been found.\n" );
            fprintf(stdout, "Objective function value: %f\n", F );
            fprintf(stdout, "Iteration: %d\n", iter);
            break;
        }
        
        // break;
        if( iter >= param.max_iter ){
            fprintf(stdout, "Maximum iteration has been reached.\n" );
            break;
        }
        
        if( norm_beta <= Gamma * norm_phi ){
            
            iter_type = 1;
            struct InputCG input_cg;
            input_cg.n = n;
            input_cg.nS = nS;
            input_cg.x = x;
            input_cg.S = S;
            input_cg.Sp = Sp;
            input_cg.Sn = Sn;
            input_cg.eta_r = param.eta_r;
            input_cg.grad_F = grad_F;
            input_cg.nVmaxAllow = ceil( MAX((double)param.nVmaxAllow,input_cg.fracViol*(double)nS) );
            input_cg.TRradius = TRradiusPhi;
            input_cg.maxCG_iter = param.max_CG_iter;
            input_cg.fracViol = param.fracViol;
            input_cg.rsType = rsType;
            input_cg.hessVecProd = &logistic_hessVecProd;
            
            logistic_setXF();
            
            struct OutputCG output_cg = CGsolver( input_cg );
            
            d = output_cg.d_full;
            
            if( output_cg.nV == 0 && norm_beta < tol_absolute ){
                rsType = 2;
            }else{
                rsType = 1;
            }
            
            // Perform a line search
            j = 0;
            alpha = 1.0;
            dirDer = output_cg.dirDer;
            double suf_descent_1 = param.eta * dirDer;
            F_old = F;
            struct projResult projresult = project( x, d, alpha );
            x_linesearch = projresult.project_vector;
            sameOrthant = projresult.same_sign;
            logistic_setExpTerm( x_linesearch, wTx, help_vector, ysigmoid, expterm );
            F = logistic_func( x_linesearch, expterm );
            
            while(1){
                if( sameOrthant && F - F_old <=  suf_descent_1 * alpha ){
                    // set x, and grad_f
                    for( i = 0; i < n; i++ ){
                        step[i] = x_linesearch[i] - x[i];
                        norm_step += step[i] * step[i];
                        x[i] = x_linesearch[i];
                        grad_f[i] = -1.0 / dm * dot_r( ysigmoid, problem.cX[i] );
                    }
                    norm_step = sqrt( norm_step );
                    break;
                }
                
                if( !sameOrthant && F < F_old ){
                    // set x, and grad_f
                    for( i = 0; i < n; i++ ){
                        step[i] = x_linesearch[i] - x[i];
                        norm_step += step[i] * step[i];
                        x[i] = x_linesearch[i];
                        grad_f[i] = -1.0 / dm * dot_r( ysigmoid, problem.cX[i] );
                    }
                    norm_step = sqrt( norm_step );
                    break;
                }
                
                if ( j > maxback ){
                    // set x, and grad_f
                    for( i = 0; i < n; i++ ){
                        step[i] = x_linesearch[i] - x[i];
                        norm_step += step[i] * step[i];
                        x[i] = x_linesearch[i];
                        grad_f[i] = -1.0 / dm * dot_r( ysigmoid, problem.cX[i] );
                    }
                    norm_step = sqrt( norm_step );
                    break;
                }
                
                alpha *= xi;
                projresult = project( x, d, alpha );
                x_linesearch = projresult.project_vector;
                sameOrthant = projresult.same_sign;
                logistic_setExpTerm( x_linesearch, wTx, help_vector, ysigmoid, expterm );
                F = logistic_func( x_linesearch, expterm );
                j++;
                
            }
            // fprintf(stdout, "CG F: %f\n", F);
            
            
        }else{
            iter_type = 0;
            double TRbeta_norm_beta = -( TRradiusBeta / norm_beta );
            for ( i = 0; i < n; i++ ){
                d[i] = TRbeta_norm_beta * beta[i];
                x_linesearch[i] = x[i] + d[i];
                if( d[i] > 0 ){
                    grad_F[i] = grad_f[i] + lambda;
                }else if( d[i] < 0 ){
                    grad_F[i] = grad_f[i] - lambda;
                }else{
                    grad_F[i] = 0.0;
                }
            }
            
            double normd = TRradiusBeta;
            // Perform line search to get updated x.
            j = 0;
            alpha = 1.0;
            dirDer = dot_n5( grad_F, d, n );
            
            F_old = F;
            logistic_setExpTerm( x_linesearch, wTx, help_vector, ysigmoid, expterm );
            F = logistic_func( x_linesearch, expterm );
            
            double suf_descent_1 = param.eta * dirDer;
            double tmp = 0.0; // used for update grad_f
            
            // perform line search
            while( 1 ){
                if( F - F_old <= alpha * suf_descent_1 ){
                    // set x, and grad_f
                    for( i = 0; i < n; i++ ){
                        step[i] = x_linesearch[i] - x[i];
                        norm_step += step[i] * step[i];
                        x[i] = x_linesearch[i];
                        grad_f[i] = -1.0 / dm * dot_r( ysigmoid, problem.cX[i] );
                    }
                    norm_step = sqrt( norm_step );
                    break;
                }
                if ( j > maxback ){
                    // set x, and grad_f
                    for( i = 0; i < n; i++ ){
                        step[i] = x_linesearch[i] - x[i];
                        norm_step += step[i] * step[i];
                        x[i] = x_linesearch[i];
                        grad_f[i] = -1.0 / dm * dot_r( ysigmoid, problem.cX[i] );
                    }
                    norm_step = sqrt( norm_step );
                    break;
                }
                j++;
                alpha *= xi;
                
                // set new trial step
                for( i = 0 ; i < n; i++ ){
                    x_linesearch[i] = x[i] + alpha * d[i];
                }
                
                logistic_setExpTerm( x_linesearch, wTx, help_vector, ysigmoid, expterm );
                F = logistic_func( x_linesearch, expterm );
            }
            
            // fprintf(stdout, "%f\n", F );
        }
        
        
        
        d = Malloc( double, n );
        nS = 0;
        // Set trust-region radius for next iteration.
        if( iter_type == 1 ){
            TRradiusPhi = MAX( 1E-3, MIN( 1E3, 10 * norm_step ) );
        }else{
            TRradiusBeta = MAX( 1E-5, MIN( 1.0, norm_step ) );
        }
        norm_step = 0.0;
        iter++;
        
    }
}

double logistic_func( double *x, double *expterm ){
    
    double funcValue = 0.0;
    for (int i = 0; i < m; i++ ){
        funcValue += log( 1.0 + expterm[i] );
    }
    
    funcValue /= (double) m;
    
    funcValue += lambda * l1_n5( x, n );
    
    return funcValue;
}

void logistic_setExpTerm( double *x, double *wTx, double *sigmoid, double *ysigmoid, double *expterm ){
    for( int i = 0; i < m; i++ ){
        wTx[i] = dot_c( x, problem.rX[i] );
        
        expterm[i] = exp( -1.0 * problem.y[i] * wTx[i] );
        help_vector[i] = expterm[i] / ( 1.0 + expterm[i] );
        ysigmoid[i] = problem.y[i] * help_vector[i];
    }
}

double *logistic_hessVecProd( double *v ){
    
    double *hv = Malloc( double, nS );
    
    int i;
    
    double *tmp_vector = Malloc( double, m );
    for( i = 0; i < m; i++ ){
        tmp_vector[i] = help_vector[i]*(1-help_vector[i])*dot_c(v,problem.rXS[i]);
    }
    
    for( i = 0; i < nS; i++ ){
        hv[i] = 1.0 / (double) m * dot_r( tmp_vector, problem.cXS[i] ) + 1E-8 * v[i];
    }
    
    free (tmp_vector);
    return hv;
}

void logistic_setXF(){
    int i;
    int nnz = 0;
    // set column format data
    problem.cXS = Malloc( struct rcNode *, nS );
    
    
    for(  i = 0; i < nS; i++ ){
        problem.cXS[i] = &X_col[col_ptr_head[S[i]]];
    }
    
    problem.rXS = Malloc( struct rcNode *, m );
    
    int *row_ptr_sp = Malloc( int, m + 1 );
    
    
    for( i = 0; i < m + 1; i++ ){
        row_ptr_sp[i] = 0;
    }
    
    for( i = 0; i < nS; i++ ){
        struct rcNode *x = problem.cXS[i];
        while( x->r_index != -1 ){
            nnz++;
            row_ptr_sp[x->r_index+1 ]++;
            x++;
        }
    }
    
    for( i = 1; i < m + 1; i++ ){
        row_ptr_sp[i] += row_ptr_sp[i-1] + 1;
    }
    
    X_row_S = Malloc( struct rcNode, nnz + m );
    
    for( i = 0; i < m; i++ ){
        problem.rXS[i] = &X_row_S[row_ptr_sp[i]];
    }
    
    for( i = 0; i < nS; i++ ){
        struct rcNode *x = problem.cXS[i];
        while( x->r_index != -1 ){
            int idx = x->r_index; // without - 1
            X_row_S[row_ptr_sp[idx]].c_index = i; // starts from 0, (1)
            X_row_S[row_ptr_sp[idx]].r_index = idx;
            X_row_S[row_ptr_sp[idx]].value = x->value;
            row_ptr_sp[idx]++;
            x++;
        }
    }
    
    for( i = 0; i < m; i++ ){
        X_row_S[row_ptr_sp[i]].c_index = -1;
        X_row_S[row_ptr_sp[i]].r_index = -1;
    }
    
}

/* **********  Logistic Regression Loss Section End ******** */


/* **********  Least Square Loss Section Start ******** */



/* **********  Least Square Loss Section End ******** */

/* Sparse Operators Start */

double dot_r( const double *v, const struct rcNode *x ){
    double vTx = 0.0;
    while( x->r_index != -1 ){
        vTx += v[x->r_index] * x->value;
        x++;
    }
    return vTx;
}


double dot_c( const double *v, const struct rcNode *x ){
    double vTx = 0.0;
    while( x->c_index != -1 ){
        vTx += v[x->c_index] * x->value;
        x++;
    }
    return vTx;
}


double norm2_sq_r( const struct rcNode *x ){
    double val = 0.0;
    while( x->r_index != -1 ){
        val =+ x->value * x->value;
        x++;
    }
    return val;
}

double norm2_sq_c( const struct rcNode *x ){
    double val = 0.0;
    while( x->c_index != -1 ){
        val =+ x->value * x->value;
        x++;
    }
    return val;
}

double dot_S_r( const double *v, const struct rcNode *x, int *S ){
    
    double vTx = 0.0;
    
    while( x->r_index != -1 && *S != -1 ){
        if( x->r_index == *S ){
            vTx += *v * x->value;
            v++;
            x++;
            S++;
        }else if( x->r_index > *S ){
            S++;
            v++;
        }else if( x->r_index < *S ){
            x++;
        }
    }
    
    return vTx;
}

double dot_S_c( const double *v, const struct rcNode *x, int *S ){
    
    double vTx = 0.0;
    
    while( x->c_index != -1 && *S != -1 ){
        if( x->c_index == *S ){
            vTx += *v * x->value;
            v++;
            x++;
            S++;
        }else if( x->c_index > *S ){
            S++;
            v++;
        }else if( x->c_index < *S ){
            x++;
        }
    }
    
    return vTx;
}

/* Sparse Operators End */

/* n5 operators start */
double dot_n5( double *v1, double *v2, int n ){
    int i, n5;
    double result = 0.0;
    
    if( n <= 0 ) return result;
    
    n5 = n % 5;
    
    for( i = 0; i < n5; i++ ){
        result += v1[i] * v2[i];
    }
    
    for( ; i < n; i += 5 ){
        result += v1[i]*v2[i] + v1[i+1]*v2[i+1] + v1[i+2]*v2[i+2] + v1[i+3]*v2[i+3] + v1[i+4]*v2[i+4];
    }
    return result;
}

double l1_n5( double *v, int n ){
    int i, n5;
    double result = 0.0;
    
    if( n <= 0 ) return result;
    
    n5 = n % 5;
    
    for( i = 0; i < n5; i++ ){
        result += fabs(v[i]);
    }
    
    for( ; i < n; i += 5 ){
        result += fabs(v[i]) + fabs(v[i+1]) + fabs(v[i+2]) + fabs(v[i+3]) + fabs(v[i+4]);
    }
    return result;
}

/* n5 operators end */

/* print tools */
void print( double *vector, int n ){
    
    for( int i = 0; i < n; i++ ){
        fprintf(stdout, "%f ", vector[i] );
    }
    fprintf(stdout, "\n" );
}

void printInt( int *vector, int n ){
    
    for( int i = 0; i < n; i++ ){
        fprintf(stdout, "%d ", vector[i] );
    }
    fprintf(stdout, "\n" );
}

/* CG reduced-space solver */
struct OutputCG CGsolver( struct InputCG input_cg ){
    
    struct OutputCG output_cg;
    
    int i;
    int nV = 0;
    int iter_cg = 0;
    int *S = input_cg.S;
    int nS = input_cg.nS;
    int n = input_cg.n;
    int nVmaxAllow = input_cg.nVmaxAllow;
    int maxCG_iter = input_cg.maxCG_iter;
    int max_loop = MIN( maxCG_iter, nS );
    double *x = input_cg.x;
    double *grad_F = input_cg.grad_F;
    double TRradius = input_cg.TRradius;
    double *d_full = Malloc( double, n );
    double *d = Malloc( double, nS );
    double *p = Malloc( double, nS );
    double *r = Malloc( double, nS );
    double *Hp = Malloc( double, nS );
    double normr0 = 0.0;
    double normr;
    double normrSq;
    double res_target = 0.0;
    double pTHp = 0.0;
    double alphaCG = 0.0;
    double betaCG = 0.0;
    double normd = 0.0;
    double norm_old = 0.0;
    
    
    for( i = 0; i < nS; i++ ){
        r[i] = grad_F[ S[i] ];
        p[i] = -grad_F[ S[i] ];
        normr0 += p[i] * p[i];
    }
    
    normrSq = normr0;
    normr0 = sqrt( normr0 );
    normr = normr0;
    
    if( input_cg.rsType == 1 ){
        res_target = MAX( input_cg.eta_r * normr0, 1E-12 );
    }else{
        res_target = MAX( MIN( input_cg.eta_r, normr0 ) * normr0, 1E-12 );
    }
    
    while(1){
        /* Compute next linear CG iterate. */
        Hp = input_cg.hessVecProd(p);
        
        pTHp = dot_n5( p, Hp, nS );
        alphaCG = normrSq / pTHp;
        norm_old = normr;
        normr = 0.0;
        for( i = 0; i < nS; i++ ){
            d[i] += alphaCG * p[i];
            normd += d[i] * d[i];
            d_full[S[i]] = d[i];
            r[i] = r[i] + alphaCG * Hp[i];
            normr += r[i] * r[i];
        }
        
        normrSq = normr;
        normd = sqrt( normd );
        normr = sqrt( normr );
        
        /* calculate violating sets */
        for( i = 0; i < n; i++ ){
            if( x[i] > 0 && x[i] + d_full[i] < 0 ){
                nV++;
            }else if( x[i] < 0 && x[i] + d_full[i] > 0){
                nV++;
            }
        }
        
        /* Check for termination of CG. */
        if( normr <= res_target ){
            /* fprintf(stdout, "CG teriminate : CG residual\n" ); */
            break;
        }else if( nV > nVmaxAllow ){
            /* fprintf(stdout, "CG teriminate : CG violate\n" ); */
            break;
        }else if( normd >= TRradius ){
            /* fprintf(stdout, "CG teriminate : CG big\n" ); */
            break;
        }else if( iter_cg > max_loop ){
            /* fprintf(stdout, "CG teriminate : CG maximum\n" ); */
            break;
        }
        
        betaCG = normrSq / ( norm_old * norm_old );
        
        for( i = 0; i < nS; i++ ){
            p[i] = -r[i] + betaCG * p[i];
        }
        nV = 0;
        normd = 0.0;
        iter_cg++;
    }
    
    /* calculate dirder */
    output_cg.dirDer = 0.0;
    for( i = 0; i < nS; i++ ){
        output_cg.dirDer += d[i] * grad_F[S[i]];
    }
    output_cg.d_full = d_full;
    output_cg.nV = nV;
    output_cg.iter = iter_cg;
    return output_cg;
}

/************** Read parses and settings ************/
void parse_command_line( int argc, char **argv ){
    
    int i;
    param.max_iter = 1000;
    param.printlevel = 2;
    param.Gamma = 1.0;
    param.eta_r = 1E-1;
    param.eta = 1E-2;
    param.xi = 0.5;
    param.tol_absolute = 1E-6;
    param.tol_relative = 1E-6;
    param.betaFrac = 1.0;
    param.phiFrac = 1.0;
    param.tryCG = 0;
    param.tryCD = 1;
    param.tryCrossOver = 0;
    param.max_CG_iter = 1000;
    param.max_CD_iter = 1000;
    param.crossOverTol = 1E-1;
    param.TRradiusPhi = 1E3;
    param.TRradiusBeta = 1E-1;
    param.maxback = 100;
    param.fracViol = 0.1;
    param.nVmaxAllow = 1000;
    param.checkResEvery = 1;
    param.termType = 2;
    param.lambdaCoeff = 1.0;
    param.prob_name = Malloc( char, max_line_size );
    param.scaling = 0;
    param.loss_type = 1;
    param.profile_file = Malloc(char, max_line_size);
    param.data_format = Malloc( char, max_line_size );
    
    for ( i = 1; i < argc; i++ ){
        if( argv[i][0] != '-' ) break;
        if( ++i >= argc ){
            perror("Invalid arguments.");
            exit(0);
        }
        switch( argv[i-1][1] ){
            case 'p':
                param.profile_file = Malloc( char, 1024 );
                strcat( param.profile_file, argv[i] );
                break;
        }
    }
}

/* read setting from profile */
void load_setting_profile(){
    
    char *attribute;
    char *value;
    FILE *fp;
    strcat( param.profile_file, "/Users/chentianyi/Documents/FaRSA/FaRSA/FaRSA/FaRSA.profile" );
    fp = fopen( param.profile_file, "r" );
    
    if( fp == NULL ){
        fprintf( stderr, "Can not open profile file...\n" );
        perror("profile_file: fopen");
        exit(1);
    }
    
    max_line_size = 1024;
    line = Malloc( char, max_line_size );
    
    while( get_line(fp) != NULL ){
        attribute = strtok( line, ":");
        
        value = strtok( NULL, " \t\n");
        
        /* If current value is null, then continue */
        if( value == NULL ){
            continue;
        }
        
        fprintf(stderr, "%s %s\n", attribute, value);
        
        if( !strcmp( attribute, "maximum iteration" ) ){
            param.max_iter = atoi(value);
            continue;
        }else if( !strcmp( attribute, "lambda coefficient" ) ){
            param.lambdaCoeff = atof(value);
            continue;
        }else if( !strcmp( attribute, "Gamma" ) ){
            param.Gamma = atof(value);
            continue;
        }else if( !strcmp( attribute, "eta" ) ){
            param.eta = atof(value);
            continue;
        }else if( !strcmp( attribute, "eta_r" ) ){
            param.eta_r = atof(value);
            continue;
        }else if( !strcmp( attribute, "xi" ) ){
            param.xi = atof(value);
            continue;
        }else if( !strcmp( attribute, "xi" ) ){
            param.xi = atof(value);
            continue;
        }else if( !strcmp( attribute, "tolerance" ) ){
            param.tol_relative = atof(value);
            param.tol_absolute = atof(value);
            continue;
        }else if( !strcmp( attribute, "betaFrac" ) ){
            param.betaFrac = atof(value);
            param.betaFrac = param.betaFrac > 1.0 ? 1.0 : param.betaFrac;
            continue;
        }else if( !strcmp( attribute, "phiFrac" ) ){
            param.phiFrac = atof(value);
            param.phiFrac = param.phiFrac > 1.0 ? 1.0 : param.phiFrac;
            continue;
        }else if( !strcmp( attribute, "tryCG" ) ){
            param.tryCG = atoi(value);
            continue;
        }else if( !strcmp( attribute, "tryCD" ) ){
            param.tryCD = atoi(value);
            continue;
        }else if( !strcmp( attribute, "tryCrossOver" ) ){
            param.tryCrossOver = atoi(value);
            continue;
        }else if( !strcmp( attribute, "max_CG_iter") ){
            param.max_CG_iter = atoi(value);
            continue;
        }else if( !strcmp( attribute, "max_CD_iter") ){
            param.max_CD_iter = atoi(value);
            continue;
        }else if( !strcmp( attribute, "crossOverTol") ){
            param.crossOverTol = atof(value);
            continue;
        }else if( !strcmp( attribute, "TRradiusPhi") ){
            param.TRradiusPhi = atof(value);
            continue;
        }else if( !strcmp( attribute, "TRradiusBeta") ){
            param.TRradiusBeta = atof(value);
            continue;
        }else if( !strcmp( attribute, "maxback") ){
            param.maxback = atof(value);
            continue;
        }else if( !strcmp( attribute, "fracViol") ){
            param.fracViol = atof(value);
            param.fracViol = param.fracViol > 1.0 ? 1.0 : param.fracViol;
            continue;
        }else if( !strcmp( attribute, "nVmaxAllow") ){
            param.nVmaxAllow = atoi(value);
            continue;
        }else if( !strcmp( attribute, "checkResEvery") ){
            param.checkResEvery = atoi(value);
            continue;
        }else if( !strcmp( attribute, "termType") ){
            param.termType = atoi(value);
            continue;
        }else if( !strcmp( attribute, "objective function type") ){
            param.loss_type = atoi(value);
            continue;
        }else if( !strcmp( attribute, "data file") ){
            strcat( param.prob_name, value );
            fprintf(stderr, "Data file %s\n", param.prob_name );
            continue;
        }
        // else if( !strcmp( attribute, "format") ){
        // 	strcat( param.data_format, value );
        // 	fprintf(stderr, "Data format %s\n", param.data_format );
        // 	continue;
        // }
        
    }
}

/* function for selecting suitable read data */
void read_problem(){
    if( param.scaling == 0 ){
        read_sparse_problem();
    }else if( param.scaling == 1 ){
        read_scaling_problem();
    }else if( param.scaling == 2 ){
        read_scaling_problem();
    }
}

/* get line function */
char* get_line( FILE *input ){
    
    int len;
    
    if( fgets( line, max_line_size, input ) == NULL )
        return NULL;
    
    /* if current max_line_size doesn't find the \n, then enlarge
     the max_line_size */
    while( strrchr( line, '\n' ) == NULL ){
        max_line_size *= 2;
        line = (char *) realloc( line, max_line_size );
        len = (int) strlen(line);
        if( fgets( line+len, max_line_size-len, input ) == NULL )
            break;
    }
    return line;
    
}

void read_sparse_problem(){
    
    int max_index;
    char *label, *idx, *val;
    size_t entries;
    
    FILE *fp = fopen( param.prob_name, "r" );
    
    if( fp == NULL ){
        fprintf( stdout, "Can't open input file %s\n", param.prob_name );
        // exit(0) behave like return 0 in main() function, exit(1) behave like return 1 
        exit(1);
    }
    
    max_line_size = 1024;
    problem.num_samples = 0;
    entries = 0;
    
    fprintf(stderr, "Line 1245 \n" );
    
    line = Malloc( char, max_line_size );
    
    
    while( get_line(fp) != NULL ){
        
        label = strtok( line, " \t\n");
        // features
        while(1){
            
            idx = strtok( NULL, ":" );
            val = strtok( NULL, " \t" );
            
            if( val == NULL || idx == NULL )
                break;
            
            entries++;      
        }
        problem.num_samples++; 		
    }
    rewind( fp );
    
    
    fprintf(stderr, "Line 1269 \n" );
    
    problem.y = Malloc( double, problem.num_samples );
    problem.rX = Malloc( struct rcNode *, problem.num_samples ); 
    
    // sort by samples
    X = Malloc( struct rcNode, entries + problem.num_samples );
    
    max_index = -1;
    
    int j = 0;
    
    fprintf(stderr, "Line 1280 \n" );
    
    for ( int i = 0; i < problem.num_samples; i++ ){
        
        get_line(fp);
        
        problem.rX[i] = &X[j];
        
        /* strtok: The C library function char *strtok(char *str, const char *delim) breaks 
         string str into a series of tokens using the delimiter delim. */
        label = strtok( line, " \t\n");
        
        if( label == NULL )
            continue;
        
        
        problem.y[i] = atof(label);
        
        while(1){
            /* The first call to strtok() returns a pointer to the first token in the string 
             pointed to by s1. Subsequent calls to strtok() must pass a NULL pointer as the 
             first argument, in order to get the next token in the string. */
            idx = strtok( NULL, ":" );
            val = strtok( NULL, " \t" );
            
            if( val == NULL || idx == NULL )
                break;
            
            X[j].c_index = atoi(idx) - 1; /* c_index is index of feature */
            X[j].r_index = i; /* r_index is index of sample */
            X[j].value = atof(val);
            
            if( max_index < X[j].c_index ){
                max_index = X[j].c_index;
            }
            j++;
            
        }	
        
        X[j].r_index = -1;
        X[j++].c_index = -1;
        
    }
    problem.num_features = max_index + 1;
    
    fprintf(stderr, "Line 1326 \n" );
    
    transpose();
    
    fprintf(stdout, "Dataset Description:\n" );
    fprintf(stdout, "Number of samples : %d\n", problem.num_samples );
    fprintf(stdout, "Number of features: %d\n", problem.num_features );
    
    fclose(fp);
    
}


void transpose(){
    
    int i;
    
    int m = problem.num_samples;
    int tmp_n = problem.num_features;
    
    // fprintf(stderr, "m: %d, n: %d\n", m, n);
    // fprintf(stderr, "Line 1346\n" );
    int nnz;
    nnz = 0;
    
    fprintf(stderr, "Line 1352\n" );
    
    col_ptr = Malloc( int, tmp_n + 1 );
    
    fprintf(stderr, "Line 1357\n" );
    
    col_ptr_head = Malloc( int, n + 1 );
    
    problem.cX = Malloc( struct rcNode *, problem.num_features );
    
    for( i = 0; i < n + 1; i++ ){
        col_ptr[i] = 0;
        col_ptr_head[i] = 0;
    }
    
    fprintf(stderr, "Line 1358\n" );
    for( i = 0; i < m; i++ ){
        struct rcNode *x = problem.rX[i];
        while( x->c_index != -1 ){
            nnz++;
            col_ptr[x->c_index+1 ]++;
            col_ptr_head[x->c_index+1]++;
            x++;
        }
    }
    
    
    for( i = 1; i < n + 1; i++ ){
        col_ptr[i] += col_ptr[i-1] + 1;
        col_ptr_head[i] += col_ptr_head[i-1] + 1;
    }
    
    
    X_col = Malloc( struct rcNode, nnz + n );
    
    for( i = 0; i < n; i++ ){
        problem.cX[i] = &X_col[col_ptr[i]];
    }
    
    for( i = 0; i < m; i++ ){
        struct rcNode *x = problem.rX[i];
        while( x->c_index != -1 ){
            int idx = x->c_index; // without - 1
            X_col[col_ptr[idx]].r_index = i; // starts from 0, (1)
            X_col[col_ptr[idx]].c_index = idx;
            X_col[col_ptr[idx]].value = x->value;
            col_ptr[idx]++;
            x++;
        }
    }
    
    for( i = 0; i < n; i++ ){
        X_col[col_ptr[i]].r_index = -1;
        X_col[col_ptr[i]].c_index = -1;
    }
    
}

void read_scaling_problem(){
    
    int max_index;
    char *label, *idx, *val;
    size_t entries;
    struct rcNode *X_tmp;
    struct rcNode **tmp_rX;
    
    FILE *fp = fopen( param.prob_name, "r" );
    
    if( fp == NULL ){
        fprintf( stdout, "Can't open input file %s\n", param.prob_name );
        exit(1);
    }
    
    max_line_size = 1024;
    problem.num_samples = 0;
    entries = 0;
    
    line = Malloc( char, max_line_size );
    
    max_index = -1;
    
    while( get_line(fp) != NULL ){
        
        label = strtok( line, " \t\n");
        /* features */
        while(1){
            
            idx = strtok( NULL, ":" );
            val = strtok( NULL, " \t" );
            
            if( val == NULL || idx == NULL )
                break;
            
            if( max_index < atoi(idx) - 1  ){
                max_index = atoi(idx) - 1;
            }
            
            entries++;      
        }
        problem.num_samples++; 		
    }
    rewind( fp );
    
    problem.num_features = max_index + 1;
    
    problem.y = Malloc( double, problem.num_samples );
    
    X_tmp = Malloc( struct rcNode, entries + problem.num_samples );
    
    tmp_rX = Malloc( struct rcNode *, problem.num_samples );
    
    fprintf(stdout, "Dataset Description:\n" );
    fprintf(stdout, "Number of samples : %d\n", problem.num_samples );
    fprintf(stdout, "Number of features: %d\n", problem.num_features );
    
    mean_column = Malloc( double, problem.num_features );
    min_column = Malloc( double, problem.num_features );
    max_column = Malloc( double, problem.num_features );
    
    int max_entries = problem.num_features * problem.num_samples;
    
    fprintf(stderr, "%d\n", max_entries);
    
    if( param.scaling == 1 ){
        entries = max_entries;
    }else if( param.scaling == 2 ){
        entries = 100;
    }
    
    fprintf( stdout, "Collecting statistic...\n" );
    
    int j = 0;
    int i;
    
    for ( i = 0; i < problem.num_samples; i++ ){
        
        get_line(fp);
        
        tmp_rX[i] = &X_tmp[j];
        /* strtok: The C library function char *strtok(char *str, const char *delim) breaks 
         string str into a series of tokens using the delimiter delim. */
        label = strtok( line, " \t\n");
        
        if( label == NULL )
            continue;
        
        problem.y[i] = atof(label);
        
        while(1){
            /* The first call to strtok() returns a pointer to the first token in the string 
             pointed to by s1. Subsequent calls to strtok() must pass a NULL pointer as the 
             first argument, in order to get the next token in the string. */
            idx = strtok( NULL, ":" );
            val = strtok( NULL, " \t" );
            
            if( val == NULL || idx == NULL )
                break;
            
            int index = atoi(idx) - 1;
            double value = atof(val);
            
            if( min_column[index] > value ){
                min_column[index] = value;
            }	
            
            if( max_column[index] < value ){
                max_column[index] = value;
            }
            
            mean_column[index] += value;
            
            X_tmp[j].c_index = index; /* c_index is index of feature */
            X_tmp[j].r_index = i; /* r_index is index of sample */
            X_tmp[j].value = value;
            
            j++;
            
        }	
        
        X_tmp[j].c_index = -1;
        X_tmp[j++].r_index = -1;
        
    }
    
    rewind( fp );
    
    for( i = 0; i < problem.num_features; i++ ){
        mean_column[i] /= (double) problem.num_samples;
    }
    
    if( param.scaling == 1 ){
        X = Malloc( struct rcNode, max_entries + problem.num_samples );
        j = 0;
        for( i = 0; i < problem.num_samples; i++ ){
            int fea_idx = 0;
            problem.rX[i] = &X[j];
            for( fea_idx = 0; fea_idx < problem.num_features; fea_idx++ ){
                X[j+fea_idx].value = ( X[j+fea_idx].value - mean_column[fea_idx] );
            }
            
            while( tmp_rX[i]->c_index != -1 ){
                tmp_rX[i]++;
            }
            
        }
    }
    
    /* select sparse or dense data format for scaling type 2 */
    if( param.scaling == 2 ){
        int count2 = 0;
        for( i = 0; i < problem.num_samples; i++ ){
            while( tmp_rX[i]->c_index != -1 ){
                if( tmp_rX[i]->value == min_column[i] ){
                    count2++;
                }
                tmp_rX[i]++;
            }
        }
        fprintf(stdout, "Scaling 2 needs to be done\n" );
        fprintf(stderr, "count2: %d\n", count2);
    }
    
    fclose( fp );
}




