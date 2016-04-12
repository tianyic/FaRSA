//
//  Utils.h
//  fPlusl1cpp
//
//  Created by chentianyi on 1/15/16.
//  Copyright (c) 2016 ChenTianyi. All rights reserved.
//

#ifndef fPlusl1cpp_Utils_h
#define fPlusl1cpp_Utils_h


#include <math.h>

double oneNorm(double *x, int n){
    double norm1 = 0.0;
    for (int i = 0; i < n; i++) {
        norm1 += x[i]>0?x[i]:-x[i];
    }
    return(norm1);
}

double Norm(double *x, int n){
    double norm2 = 0.0;
    for (int i = 0; i < n; i++) {
        norm2 += x[i]*x[i];
    }
    return(sqrt(norm2));
}

double dotProduct(double *v1,double *v2, int n){
    double value = 0.0;
    for (int i = 0; i < n; i++) {
        value += v1[i]*v2[i];
    }
    return(value);
}

void axPlusv(const double a, const double *x, double *v, int n){
    int i;
    for (i = 0; i < n; i++) {
        v[i] += a*x[i];
    }
}

void axPlusby(const double a, const double *x, const double b, double *v, int n){
    int i;
    for (i = 0; i < n; i++) {
        v[i] = a*x[i] + b*v[i];
    }
}

void Copy(double *v, double *v1, double a, int n){
    int i;
    for (i = 0; i < n; i++) {
        v1[i] = a*v[i];
    }
}

void Copy_int(int *v, int *v1, int a, int n){
    int i;
    for (i = 0; i < n; i++) {
        v1[i] = a*v[i];
    }
}
//void QuikSort(int arr[], int low, int high){
//    if(low>=high){
//        return;
//    }
//
//    int pivot = arr[low];
//    int i=low;
//    for(int j=low+1;j<=high;j++){
//        if(arr[j]<=pivot){        //a[j] is smaller than pivot
//            i++;    //a[i] is bigger than pivot
//            if(i!=j){
//                Swap(arr[i],arr[j]);
//            }
//        }
//    }
//    Swap(arr[low],arr[i]);    //Swap pivot to middle position
//    QuikSort(arr,low,i-1);        //a[i] is the pivot now
//    QuikSort(arr,i+1,high);
//}
#endif
