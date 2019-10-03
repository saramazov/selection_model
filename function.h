#ifndef __function_h__
#define __function_h__

#include <cmath>
#include <iostream>
#include <vector>
#include <ctime>
#include <fstream>
#include <random>

using namespace std;


double rand2();

int uniform(int a, int b);

int bernoulli();

double* probability(double** x,int j,int n,double* sel,int n_sel,double* sigma_sel, double* prob, int choice2, double a, double s, int conf, int bins, double b, int power);

double select(double* prob,int n);

double blending(double* prob,int n,int models,double** x,double* w,int j);

double attract(double** x,int i,int j,double* att,int n_att,double* sigma_att,int choice,int det,double* threshold,double* charge,double* elastic,double limi,double* range,int* strength,double* coefficient);

double gauss(double mean,double sigma);

#endif
