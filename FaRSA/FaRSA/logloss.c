//
//  logloss.c
//  FaRSA
//
//  Created by chentianyi on 10/12/16.
//  Copyright (c) 2016 ChenTianyi. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include "farsa.h"

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define len(array) ( sizeof(array) / sizeof( array[0]) )
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif


int main( int argc, char **argv ){
    
    
    // parse_command_line( argc, argv, &param );
    
    // read_problem( param, &logloss );
    
    // param.loss_type = "logistic_loss";
    /* For testing */
    farsa( argc, argv, NULL );
    
    return 1;
    
}

