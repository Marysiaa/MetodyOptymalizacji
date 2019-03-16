#ifndef OPT_ALG_H
#define OPT_ALG_H

#include"solution.h"
#include<random>
#include<chrono>

double *expansion(double x0, double d, double alpha, int Nmax);

solution fib(double a, double b, double epsilon);

solution lag(double a, double b, double epsilon, int Nmax, double gamma);

solution HJ(matrix x0, double s, double alpha, double epsilon, int Nmax);

solution HJ_trial(solution XB, double s);

solution Rosenbrock(matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax);

solution sym_NM(matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix A, bool wewnetrzna = false);

solution penalty(matrix x0, double c, double a, double epsilon, int Nmax, bool wewnetrzna = false);

solution golden(double a, double b, double epsilon, int Nmax, matrix P = 0.0);

solution SD(matrix x0, double epsilon, int Nmax, double h0 = 0);

solution Newton(matrix x0, double epsilon, int Nmax, double h0 = 0);

solution CG(matrix x0, double epsilon, int Nmax, double h0 = 0);

double compute_b(matrix x, matrix d);

solution Powell(matrix x0, double epsilon, int Nmax, matrix w = 0.0);

double * compute_ab(matrix x, matrix d);

solution EA(int N, matrix lb, matrix ub, double epsilon, int Nmax, double sigma = 1);

#endif