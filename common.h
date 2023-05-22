#pragma once
#include <vector>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <thread>

#define phi 1.61803398874989484820
typedef std::vector<double>   vec_n;
typedef std::vector<vec_n>    mat_mn;
typedef double(*func_1d)(const double&);
double alpha = 1.0;
typedef double(*func_nd)(const vec_n&);
#define N_DIM_ACCURACY        1e-6
#define N_DIM_DERIVATIVE_STEP 1e-6
#define N_DIM_ITERS_MAX       1000
#define ACCURACY              1e-6
#define ITERS_MAX             50
#define MAX_DENOMINATOR       1000
#define MAX(A, B) (A > B ? A : B)
#define MIN(A, B) (A < B ? A : B)
#define SIMPLEX_MAX  0
#define SIMPLEX_MIN  1
#define DISPLAY_PROGRES _DEBUG? 1 : 0
const char* LINE_UP = "\033[1A";
const char* LINE_CLEAR = "\x1b[2K\r";
#define GOLDEN_SECTION  ((1.0 + sqrt(5.0)) * 0.5)


#pragma once
