//
//  DataStream.h
//  fPlusl1cpp
//
//  Created by chentianyi on 1/13/16.
//  Copyright (c) 2016 ChenTianyi. All rights reserved.
//

#ifndef fPlusl1cpp_DataStream_h
#define fPlusl1cpp_DataStream_h

#include "LogRegCost.h"
#include "typeDef.h"

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
#define INF HUGE_VAL

static char *line = NULL;
static int max_line_len;

struct Feature_node *X_space;

static char* readLine(FILE *input){
    int len;
    
    if(fgets(line,max_line_len,input) == NULL)
        return NULL;
    
    while(strrchr(line,'\n') == NULL){
        max_line_len *= 2;
        line = (char *) realloc(line,max_line_len);
        len = (int) strlen(line);
        if(fgets(line+len,max_line_len-len,input) == NULL)
            break;
    }
    return line;
}

void readProblem(const char *filename, LogRegCost &fun){
    int max_index, sample_max_index; // use for find the num_features
    int elements;

    FILE *fp = fopen(filename,"r");
    char *endptr; // use for functions: strtod, strtol
    char *index, *value, *label;
    
    if(fp == NULL){
        fprintf(stderr,"can't open input file %s\n",filename);
        exit(1);
    }
    
    fun.num_samples=0;
    max_line_len = 1024;
    line = Malloc(char,max_line_len);
    
    while(readLine(fp)!=NULL){
        char *p = strtok(line," \t"); // label
        
        // features
        while(1){
            p = strtok(NULL," \t");
            
            // check '\n' as ' ' after the last feature
            if(p == NULL || *p == '\n'){
                break;
            }

            elements++;
        }
        elements++; // for bias term
        //std::cout<<elements<<std::endl;
        fun.num_samples++;
    }

    
    rewind(fp);
    
    fun.y = Malloc(double, fun.num_samples);
    fun.X = Malloc(struct Feature_node *, fun.num_samples);
    //std::cout<<elements<<std::endl;
    //std::cout<<elements+fun.num_samples<<std::endl;
    X_space = Malloc(struct Feature_node, elements+fun.num_samples);

    int j=0;
    max_index=0;
    
    for (int i=0; i < fun.num_samples; i++) {
        
        sample_max_index=0;
        
        readLine(fp);
        //std::cout<<line<<std::endl;
        fun.X[i]=&X_space[j];
        
        label = strtok(line," \t\n");
        
        fun.y[i] = strtod(label,&endptr);
        
        while (true) {
            
            index = strtok(NULL,":");
            value = strtok(NULL," \t");
            
            if (index == NULL || value == NULL) {
                break;
            }
            
            X_space[j].index = (int) strtol(index,&endptr,10);
            X_space[j].value = strtod(value,&endptr);
            if(X_space[j].index > sample_max_index){
                sample_max_index = X_space[j].index;
            }
            
            j++;
        }
        
        if(sample_max_index > max_index){
            max_index = sample_max_index;
        }
        
        X_space[j].value = fun.bias;
        X_space[j].index = -1;
        j++;
    }
    
    fun.num_features = max_index;
    
    fclose(fp);
}

#endif
