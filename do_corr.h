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
#include "calcs.h"
#include "initialization.h"
#include "read_dat.h"
#include "corr.h"
#include "corels.h"
#include "bins_calcs.h"
#include "outp.h"
#include <mpi.h>
using namespace std;

int do_corr(data_info &data_inf);

int do_bining(bins &B,data_info &data_inf);
