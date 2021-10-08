//Do not edit the code below (unless you know what you are doing)

#ifndef OPT_ALG_H
#define OPT_ALG_H

#include"solution.h"
#include<chrono>

#if LAB_NO>1
double *expansion(double x0, double d, double alpha, int Nmax, matrix *ud = nullptr, matrix *ad = nullptr);
solution fib(double a, double b, double epsilon, matrix *ud = nullptr, matrix *ad = nullptr);
solution lag(double a, double b, double epsilon, double gamma, int Nmax, matrix *ud = nullptr, matrix *ad = nullptr);
#endif
#if LAB_NO>2
solution HJ(matrix x0, double s, double alpha, double epsilon, int Nmax, matrix *ud = nullptr, matrix *ad = nullptr);
solution HJ_trial(solution XB, double s, matrix *ud = nullptr, matrix *ad = nullptr);
solution Rosen(matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix *ud = nullptr, matrix *ad = nullptr);
#endif
#if LAB_NO>3
solution pen(matrix x0, double c, double dc, double epsilon, int Nmax, matrix *ud = nullptr, matrix *ad = nullptr);
solution sym_NM(matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix *ud = nullptr, matrix *ad = nullptr);
#endif
#if LAB_NO>4
solution SD(matrix x0, double h0, double epsilon, int Nmax, matrix *ud = nullptr, matrix *ad = nullptr);
solution CG(matrix x0, double h0, double epsilon, int Nmax, matrix *ud = nullptr, matrix *ad = nullptr);
solution Newton(matrix x0, double h0, double epsilon, int Nmax, matrix *ud = nullptr, matrix *ad = nullptr);
solution golden(double a, double b, double epsilon, int Nmax, matrix *ud = nullptr, matrix *ad = nullptr);
#endif
#if LAB_NO>5
solution Powell(matrix x0, double epsilon, int Nmax, matrix *ud = nullptr, matrix *ad = nullptr);
#endif
#if LAB_NO>6
solution EA(int N, matrix limits, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix *ud = nullptr, matrix *ad = nullptr);
#endif

#endif