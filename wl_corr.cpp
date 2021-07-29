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

int wl_corr(data_info &data_inf)
{
  int mpi_num = MPI::COMM_WORLD.Get_size ( );
  int mpi_id = MPI::COMM_WORLD.Get_rank ( );

  string out_file="a";

  gal *gl=new gal[data_inf.n_shape];    //array of galaxies for Data set
  gal *Rl=new gal[data_inf.n_shape_rand];    //array of galaxies for randoms

  bins B;
  bins_ini(B,data_inf);

  bins B_cl;
  bins_ini(B_cl,data_inf);

  int bbye=0;

  data_inf.do_SS=1;

  if (mpi_id==0)
    {
      cout<<"working with data and randoms"<<endl;
    }
  
  bbye=read_data(gl,data_inf);
  data_inf.do_SS=0;
  data_inf.do_RsRs=1;
  bbye=read_data(Rl,data_inf);  
  data_inf.do_RsRs=0;
  if(bbye!=0)
    {
      cout<<"wl_corr:: OOPS!! something went wrong. Bye"<<endl;
      return 1;
    }     
  cout<<"wl_corr data intialised. No of lens gals="<<data_inf.n_shape<<endl;
  cout<<"wl_corr randoms intialised. No of lens gals="<<data_inf.n_shape_rand<<endl;
  
  double start2=omp_get_wtime(); // to get hang of run time
  long sd_pairs=0;
  sd_pairs+=wl_corels(gl,Rl,B.SD,B.DRs,data_inf);
  delete []gl;
  delete []Rl;  
  // double n_ph=data_inf.n_density;
  //  data_inf.n_density=0;//to prevent double counting when doing rands

  /*
  data_inf.do_SS=0;
data_inf.do_RsRs=1;
    start2=omp_get_wtime(); // to get hang of run time
  long rd_pairs=0;
  rd_pairs+=wl_corels(Rl,B.DRs,data_inf);
  delete []Rl;
  data_inf.do_RsRs=0;
  */
  if(data_inf.coordinates==2||data_inf.coordinates==4)
    int_calc(B.SD.b,B.DRs.b,data_inf);

  bool sub=false;
  add_bins(B_cl.DRs.b,B.DRs.b,B_cl.DRs.b,data_inf,sub);  
  add_bins(B_cl.SD.b,B.SD.b,B_cl.SD.b,data_inf,sub);

  if (data_inf.do_jk)
    {
      for(int i=0;i<data_inf.n_jk;i++)
	{
	  if(data_inf.coordinates==2||data_inf.coordinates==4)
	    int_calc(B.SD.jk[i].b,B.DRs.jk[i].b,data_inf);

	  add_bins(B_cl.SD.jk[i].b,B.SD.jk[i].b,B_cl.SD.jk[i].b,data_inf,sub);
	  add_bins(B_cl.DRs.jk[i].b,B.DRs.jk[i].b,B_cl.DRs.jk[i].b,data_inf,sub);
	}
    }

  final_calc_mean_projected(B.SD.b,B.DRs.b,B.final_cross.b,data_inf); // calculations for p bins.
  final_calc(B_cl.SD.b,B_cl.DRs.b,B_cl.final_cross.b,data_inf); // calculations for p bins.

  out_file=data_inf.out_file+"bins_wl"+".dat";
  outp_bins(B.final_cross.b,data_inf,out_file);

  out_file=data_inf.out_file+"bins_cl"+".dat";
  outp_bins(B_cl.final_cross.b,data_inf,out_file);

  out_file=data_inf.out_file+"bins_ls"+".dat";
  outp_bins(B.SD.b,data_inf,out_file);

  out_file=data_inf.out_file+"bins_Rs"+".dat";
  outp_bins(B.DRs.b,data_inf,out_file);
  if(data_inf.coordinates==2||data_inf.coordinates==4)
    {
      out_file=data_inf.out_file+"bins2D_wl"+".dat";
      outp_2Dbins(B.final_cross.b,data_inf,out_file);

      out_file=data_inf.out_file+"bins2D_cl"+".dat";
      outp_2Dbins(B_cl.final_cross.b,data_inf,out_file);

      out_file=data_inf.out_file+"bins2D_ls"+".dat";
      outp_2Dbins(B.SD.b,data_inf,out_file);

      out_file=data_inf.out_file+"bins2D_Rs"+".dat";
      outp_2Dbins(B.DRs.b,data_inf,out_file);
    }

  if(data_inf.do_jk)
    {
      for(int i=0;i<data_inf.n_jk;i++)
	{
	  final_calc_mean_projected(B.SD.jk[i].b,B.DRs.jk[i].b,B.final_cross.jk[i].b,data_inf); // calculations for p bins.
	  final_calc(B_cl.SD.jk[i].b,B_cl.DRs.jk[i].b,B_cl.final_cross.jk[i].b,data_inf);

	  std::ostringstream s;
	  s<<i;
	  out_file=data_inf.out_file+"bins_wl_jk"+s.str()+".dat";
	  outp_bins(B.final_cross.jk[i].b,data_inf,out_file);

	  out_file=data_inf.out_file+"bins_cl_jk"+s.str()+".dat";
          outp_bins(B_cl.final_cross.jk[i].b,data_inf,out_file);

	  out_file=data_inf.out_file+"bins_ls_jk"+s.str()+".dat";
          outp_bins(B.SD.jk[i].b,data_inf,out_file);
	  out_file=data_inf.out_file+"bins_Rs_jk"+s.str()+".dat";
          outp_bins(B.DRs.jk[i].b,data_inf,out_file);
	  if(data_inf.coordinates==2||data_inf.coordinates==4)
	    {
	      out_file=data_inf.out_file+"bins2D_wl_jk"+s.str()+".dat";
	      outp_2Dbins(B.final_cross.jk[i].b,data_inf,out_file);

	      out_file=data_inf.out_file+"bins2D_cl_jk"+s.str()+".dat";
	      outp_2Dbins(B_cl.final_cross.jk[i].b,data_inf,out_file);

	      out_file=data_inf.out_file+"bins2D_ls_jk"+s.str()+".dat";
	      outp_2Dbins(B.SD.jk[i].b,data_inf,out_file);
	      out_file=data_inf.out_file+"bins2D_Rs_jk"+s.str()+".dat";
	      outp_2Dbins(B.DRs.jk[i].b,data_inf,out_file);
	    }
	}

      bin bins_final_jk[data_inf.n_bins];
      bin cl_bins_final_jk[data_inf.n_bins];
      if(data_inf.lin_bin) // linear bins intialisation
	{
	  bin_lin(bins_final_jk,data_inf);
	  bin_lin(cl_bins_final_jk,data_inf);
	}
      else{
	bin_log(bins_final_jk,data_inf);
	bin_log(cl_bins_final_jk,data_inf);
      }
      jk_final(bins_final_jk,B.final_cross.b,B.final_cross.jk,data_inf);
      jk_final(cl_bins_final_jk,B_cl.final_cross.b,B_cl.final_cross.jk,data_inf);
      out_file=data_inf.out_file+"bins_wl_jk_final.dat";
      outp_bins(bins_final_jk,data_inf,out_file);
      out_file=data_inf.out_file+"bins_cl_jk_final.dat";
      outp_bins(cl_bins_final_jk,data_inf,out_file);
    }

  /*  out_file=data_inf.out_file+"gal_lens.dat";
  outp_gal(gl,data_inf.n_lens,data_inf,out_file);
  */
  cout<<"Time taken== "<<omp_get_wtime()-start2<<endl;
  cout<<"Calcs Done!!"<<mpi_id<<endl;

  return 0;
}
