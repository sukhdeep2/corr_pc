//calculating correlations

/*
How it all works:

function corels1 calculate DD or RR, avoids double counting. Corels2 count DR. Finally give sum(x) in each bin and sum(x^2) in each error bin which is sorted later on.

corel called from corels1 or 2 just does calculations for a pair of galaxies. measures their angular diameter distance(dA_r) and average line of sight distance(dP). Then calculates e+ eX along line joining the two galaxies and e+^2,e++ etc. from that and bins all the quantities in dA_r and dP.


*/

#include <omp.h>
#include <iostream>
#include <stdio.h>
#include <math.h>
#define _USE_MATH_DEFINES
#include <time.h>
#include "data_def.h"
#include "calcs.h"
#include "initialization.h"
#include "corels.h"
#include "outp.h"
#include "read_dat.h"
using namespace std;

long corels1(gal g1[], bin_jk &bjk,data_info &data_inf);
//correlations for 1 data set

long corels2(gal g1[],gal g2[], bin_jk &bjk,data_info &data_inf);
// cross correlatio

long wl_corels(gal l[],gal R[], bin_jk &bjk,bin_jk &bjkR, data_info &data_inf);
