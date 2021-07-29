#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include "data_def.h"
#include "calcs.h"
using namespace std;


int read_data(gal g[],data_info &data_inf);// read data sets from files

int read_patch(int patch_id, int &patch_gal_size_max, gal* &patch_gal, data_info data_inf); // read patch data sets from Rachel's catalog
