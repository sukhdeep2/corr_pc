#include <omp.h>
#include <iostream>
#include <stdio.h>
#include <math.h>
#define _USE_MATH_DEFINES
#include <time.h>
#include "data_def.h"
#include "calcs.h"
#include "initialization.h"
#include "outp.h"

using namespace std;


int density_density_corel(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct); // correlations b/w 2 galaxies
int density_convergence_corel(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct);//wl type correlations
int shape_convergence_corel(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct);//wl type correlations

int shape_shape_corel(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct);
int shape_density_corel(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct);
int ED_corel(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct);