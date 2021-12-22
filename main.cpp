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
// #include "sky_calcs.h"
// #include "PB_calcs.h"
#include "initialization.h"
#include "corels.h"
#include "do_corr.h"
#include "wl_corr.h"
#include "bins_calcs.h"
#include "outp.h"
#include <mpi.h>

using namespace std;

//double pi=3.14159265359;

int main (int argc,char *argv[])
{

  cout<<"Beginning corr_pc. Last compiled on: "<<__DATE__<<"  at: "<<__TIME__<<endl;
  MPI::Init ( argc, argv );

  string args[argc];

  for (int i=0;i<argc;i++)
    {
      std::stringstream ss (argv[i]); // to conver from *char to double
      ss>>args[i];
    }

  data_info data_inf;
  int bbye=0;

  if(argc>1)
    bbye=data_ini(data_inf,args[1]); // read input

  else
    {
      cout<<"no input file given!!"<<endl;
      bbye=1;
    }

  if(bbye!=0)
    {
      cout<<"Bad initialization!! Bye"<<endl;
      return 1;
    }

  else if(MPI::COMM_WORLD.Get_rank ( )==0)
    {
      string FileName=data_inf.out_file+"inp.info";
      bbye=outp_data_info(data_inf,FileName);
    }

  if(bbye!=0)
    {
      cout<<"OOPS!! Something went wrong... Bye"<<endl;
      return 1;
    }

  int corr=0;
  if (data_inf.which_corr<6)
    corr=do_corr(data_inf);
  else if (data_inf.which_corr>=7)
    corr=wl_corr(data_inf);
  else
    cout<<"Wrong choice of corr... Bye"<<endl;
  MPI::Finalize();
  string FileName=data_inf.out_file+"inp2.info";
  bbye=outp_data_info(data_inf,FileName);
  return 0;
}
