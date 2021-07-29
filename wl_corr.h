#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <string>
#include <sstream>
#include <omp.h>
#include <pthread.h>
#include "data_def.h"
#include "calcs.h"
#include "initialization.h"
#include "corr.h"
#include "bins_calcs.h"
#include "outp.h"
#include <mpi.h>
using namespace std;

int wl_corr(data_info &data_inf);
