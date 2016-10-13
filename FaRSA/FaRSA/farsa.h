//
//  farsa.h
//  FaRSA
//
//  Created by chentianyi on 10/12/16.
//  Copyright (c) 2016 ChenTianyi. All rights reserved.
//

#ifndef __FaRSA__farsa__
#define __FaRSA__farsa__


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif


/* Sparse data node that contains the row index, column index and the value of data point*/
struct rcNode
{
    int r_index; /* r_index is the index of row that is index of sample */
    int c_index; /* c_index is the index of column that is index of feature */
    double value;
};

/* Problem structure to save information */
struct Problem
{
    int num_features;
    int num_samples;
    double *y;
    struct rcNode **rX;
    struct rcNode **cX;
    struct rcNode **cXS;
    struct rcNode **rXS;
};

/* Parameter structure */
struct Parameter
{
    int max_iter;
    int printlevel;
    double Gamma;
    double eta_r;
    double eta;
    double xi;
    double tol_absolute;
    double tol_relative;
    double betaFrac;
    double phiFrac;
    int tryCG;
    int tryCD;
    int tryCrossOver;
    int max_CG_iter;
    int max_CD_iter;
    int scaling;
    double crossOverTol;
    double TRradiusPhi;
    double TRradiusBeta;
    int maxback;
    double fracViol;
    int nVmaxAllow;
    int checkResEvery;
    int termType;
    double lambda;
    char *prob_name;
    double lambdaCoeff;
    double ttol;
    int loss_type;
    char *profile_file;
    char *data_format;
};

/* Input argument for farsa */
struct Input_FaRSA{
    double ( *func )( double * );
    double *( *grad_f )( double * );
    double *( *hessVec )( double * );
};

/* Projection result */
struct projResult{
    double *project_vector;
    int same_sign;
};

/* Subspace solver part starts */
struct InputCG{
    double *x;
    int *S;
    int *Sp;
    int *Sn;
    int nS;
    int nVmaxAllow;
    double *grad_F;
    double eta_r;
    double rsType;
    double fracViol;
    double TRradius;
    double HvProds;
    int maxCG_iter;
    int n;
    double *( *hessVecProd )( double *);
};

struct OutputCG{
    double *d_full;
    double dirDer;
    int nV;
    int iter;
    
};

struct InputCD{
    double *x;
    int *S;
    int *Sp;
    int *Sn;
    int nS;
    int nVmaxAllow;
    double *grad_F;
    double eta_r;
    double rsType;
    double fracViol;
    double TRradius;
    double HvProds;
    int maxCG_iter;
    int n;
    double *( *hessVecProd )( double * );
};

struct OutputCD{
    double *d_full;
    double dirDer;
    int nV;
    int iter;
    
};

/* Subspace solver part ends */

/* Global variables */
int iter;
int beta_iter;
int phi_iter;
int m;
int n;
int rsType; /* Target Residual Type in phi step */
int iter_type; /* 1 is phi step, 0, is beta step */
int *S;
int *Sz;
int *Sn;
int *Sp;
int n5, nS; /* n5 use for accelerate calculation process */
double *beta;
double *phi;
double *grad_f;
double *grad_F;
double norm_beta;
double norm_phi;
double norm_beta0;
double norm_phi0;
double lambda;
double Gamma;
double *x;
double ttol;
double tol_relative;
double tol_absolute;
double *help_vector;
double *y;

double dm;
double *d;
double TRradiusBeta;
double TRradiusPhi;
double alpha;
double *x_linesearch;
double dirDer;
double F_old;  /* objective function value for last step */
double F; /* objective function value on current step */
double *step; /* difference between two steps for trust region update */
double norm_step;
int maxback;
double xi;
int sameOrthant;

struct rcNode *X;
struct rcNode *X_col;
struct rcNode *X_row_S;

int *col_ptr;
int *col_ptr_head;
double *mean_column;
double *min_column;
double *max_column;

static int max_line_size = 0;
static char *line = NULL;


/* global problem for simplifying argument in coding */
struct Problem problem;

/* global parameter structure for setting information */
struct Parameter param;

/* Load information from parser arguments */
void parse_command_line( int argc, char **argv );

/* Load information from profile */
void load_setting_profile();

char* get_line( FILE *input );

/* select the correct way to read problem */
void read_problem();

/* read data by the libsvm format */
void read_sparse_problem();

/* transpose libsvm format data */
void transpose();

/* read dense type data */
void read_dense_problem();

/* read libsvm format data with scaling techniques */
void read_scaling_problem();


void farsa( int argc, char **argv, struct Input_FaRSA *input_farsa );

/* reduced space CG solver */
struct OutputCG CGsolver( struct InputCG input_cg );

/* reduced space CD solver */
struct OutputCD CDsolver( struct InputCD input_cg );

/* Projected gradient descent */
struct projResult project( double *x, double *d, double alpha );


void setBetaPhi();


void logistic_loss( const struct Parameter param );

/* calculate logistic loss function value */
double logistic_func( double *x, double *expterm );

/* update some intermediate variables of logistic loss */
void logistic_setExpTerm(  double *x, double *wTx, double *sigmoid, double *ysigmoid, double *expterm );

/* logistic loss hess vector product */
double *logistic_hessVecProd( double *v );

void logistic_setXF();

void generic_loss( const struct Parameter param, struct Input_FaRSA *input_farsa );


/* Sparse operators start */
/* Operator 1: the inner product of one vector and a sparse vector with the same dimension via row */
double dot_r( const double *v, const struct rcNode *x );

/*  Operator 2: the inner product of one vector and a sparse vector with the same dimension via column */
double dot_c( const double *v, const struct rcNode *x );

/*  Operator 3: l2 norm square by row */
double norm2_sq_r( const struct rcNode *x );

/*  Operator 4: l2 norm square by column */
double norm2_sq_c( const struct rcNode *x );

/*  Operator 5: inner product of v and subvector of the sparse vector given index S by row. */
double dot_S_r( const double *v, const struct rcNode *x, int *S );

/*  Operator 6: inner product of v and subvector of the sparse vector given index S by column. */
double dot_S_c( const double *v, const struct rcNode *x, int *S );

/* Sparse operators end */


/* n5 operators start */
/* Operator 7: n5 type inner product */
double dot_n5( double *v1, double *v2, int n );

/* Operator 8: n5 type l1 norm */
double l1_n5( double *v, int n );

/* n5 operators end */

/* print tools */
void print( double *vector, int n );

void printInt( int *vector, int n );


#endif /* defined(__FaRSA__farsa__) */
