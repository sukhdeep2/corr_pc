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
#include "read_dat.h"
#include "corr.h"
#include "corels.h"
#include "bins_calcs.h"
#include "do_corr.h"
#include "jk.h"
#include "outp.h"
#include <mpi.h>
using namespace std;

int do_corr(data_info &data_inf)
{
  int mpi_num = MPI::COMM_WORLD.Get_size ( );
  int mpi_id = MPI::COMM_WORLD.Get_rank ( );

  string out_file="a";

  gal *gshape;//=new gal[data_inf.n_shape];    //array of galaxies for Data set
  gal *grand_density;//=new gal[data_inf.n_density_rand];   //randoms
  gal *grand_shape;//=new gal[data_inf.n_shape_rand];   //randoms
  gal *gdensity;//=new gal[data_inf.n_density];
  int bbye=0;
  bool RsRd_same=data_inf.RsRd_same;//data_inf.shape_rand_pos_file.compare(data_inf.density_rand_pos_file);
  bins B;
  bins_ini(B,data_inf);

  data_inf.do_SS=1;
  gshape=new gal[data_inf.n_shape];
  bbye=read_data(gshape,data_inf);
  data_inf.do_SS=0;

  if(!data_inf.periodic_box||(data_inf.periodic_box && data_inf.n_jk>0)){
    data_inf.do_RsRs=1;
    grand_shape=new gal[data_inf.n_shape_rand];
    bbye+=read_data(grand_shape,data_inf);
    data_inf.do_RsRs=0;
  }
  if(bbye!=0)
    {
      cout<<"auto_corr:: OOPS!! something went wrong. Bye"<<endl;
      return 1;
    }

  long ss_pairs=0,sr_pairs=0,rr_pairs=0,sd_pairs=0,dr_pairs=0;

  if(data_inf.do_auto_corr)
    {
      data_inf.do_SS=1;
      cout<<"Doing auto SS"<<endl;
      ss_pairs+=corels1(gshape,B.SS,data_inf);
      data_inf.do_SS=0;
    }
  if((data_inf.do_auto_corr||RsRd_same)&& (!data_inf.periodic_box)) //||(data_inf.periodic_box && data_inf.n_jk>0)))
    {
      data_inf.do_SRs=1;
      cout<<"Doing auto SR"<<endl;
      sr_pairs+=corels2(grand_shape,gshape,B.SRs,data_inf);
      data_inf.do_SRs=0;

      data_inf.do_RsRs=1;
      cout<<"Doing auto RR"<<endl;
      if(!data_inf.RsRd_same)
	     rr_pairs+=corels1(grand_shape,B.RsRs,data_inf); //NOT VALID FOR CROSS CORRELATION... ONE SIDED P IS PROBLEM
      else
	     rr_pairs+=corels2(grand_shape,grand_shape,B.RsRd,data_inf);
      data_inf.do_RsRs=0;
    }

  if(data_inf.do_cross_corr)
    {
      data_inf.do_DD=1;
      gdensity=new gal[data_inf.n_density];
      bbye+=read_data(gdensity,data_inf);
      data_inf.do_DD=0;

      cout<<!RsRd_same<<!data_inf.periodic_box<<(!RsRd_same && (!data_inf.periodic_box))<<endl;

      if(!RsRd_same)// && (!data_inf.periodic_box))//||(data_inf.periodic_box && data_inf.n_jk>0)))
	{
	  data_inf.do_RsRd=1;
	  grand_density=new gal[data_inf.n_density_rand];
	  bbye+=read_data(grand_density,data_inf);
	  data_inf.do_RsRd=0;

	  if(bbye!=0)
	    {
	      cout<<"cross_corr:: OOPS!! something went wrong. Bye"<<endl;
	      return 1;
	    }

	  cout<<"do_corr: here"<<endl;

	  sr_pairs=0;
	  data_inf.do_SRd=1;
	  cout<<"Doing cross SRd"<<endl;
	  //sr_pairs+=corels2(gshape,grand_density,B.SRd,data_inf);
	  sr_pairs+=corels2(grand_density,gshape,B.SRd,data_inf);
	  data_inf.do_SRd=0;

	  data_inf.do_RsRd=1;
	  cout<<"Doing cross RsRd"<<endl;
	  rr_pairs+=corels2(grand_density,grand_shape,B.RsRd,data_inf);
	  //	  out_file=data_inf.out_file+"gal_density_rand.dat";
	  //outp_gal(grand_density,data_inf,out_file);
	  data_inf.do_RsRd=0;
	}
      else
	     grand_density=grand_shape;

      data_inf.do_SD=1;
      cout<<"Doing cross SD"<<endl;
      sd_pairs+=corels2(gdensity,gshape,B.SD,data_inf);
      data_inf.do_SD=0;
      if(!data_inf.periodic_box)//||(data_inf.periodic_box && data_inf.n_jk>0))
	     {
      	  data_inf.do_DRs=1;
      	  cout<<"Doing cross DRs"<<endl;
      	  dr_pairs+=corels2(grand_shape,gdensity,B.DRs,data_inf); //this is symmetric
      	  data_inf.do_DRs=0;
	      }
    }
  if(RsRd_same && (!data_inf.periodic_box))//||(data_inf.periodic_box && data_inf.n_jk>0)))
    {
      bool sub=false;
      add_bins(B.RsRs.b,B.RsRd.b,B.RsRs.b,data_inf,sub);
      add_bins(B.SRd.b,B.SRs.b,B.SRd.b,data_inf,sub);
      data_inf.RsRs_wt.val*=2.0;//data_inf.RsRd_wt.val;
      data_inf.RsRs_wt.err*=2.0;//data_inf.RsRd_wt.err;
      //      data_inf.SRd_wt.val=data_inf.SRs_wt.val;
      // data_inf.SRd_wt.err=data_inf.SRs_wt.err;
      data_inf.Rd_wt.err=data_inf.Rs_wt.err;
      data_inf.Rd_wt.val=data_inf.Rs_wt.val;
      if (data_inf.do_jk)
	{
	  for (int i=0;i<data_inf.n_jk;i++)
	    {
	      add_bins(B.RsRs.jk[i].b,B.RsRd.jk[i].b,B.RsRs.jk[i].b,data_inf,sub);
	      add_bins(B.SRd.jk[i].b,B.SRs.jk[i].b,B.SRd.jk[i].b,data_inf,sub);
	    }
	}
    }


  do_bining(B,data_inf);
  cout <<"Binning done"<<endl;

  string FileName=data_inf.out_file+"inp2.info";
  bbye=outp_data_info(data_inf,FileName);

  if(data_inf.do_auto_corr)
    {
      data_inf.do_SS=1;
      //out_file=data_inf.out_file+"gal_shape.dat";
      //outp_gal(gshape,data_inf,out_file);
      data_inf.do_SS=0;

      delete []gshape;
      if(!data_inf.periodic_box)
        {
          data_inf.do_RsRs=1;
	  //  out_file=data_inf.out_file+"gal_shape_rand.dat";
	  //outp_gal(grand_shape,data_inf,out_file);
          data_inf.do_RsRs=0;
	  delete []grand_shape;
        }
    }

  if (data_inf.do_cross_corr)
    {
      if(!data_inf.do_auto_corr)
	{
	  data_inf.do_SS=1;
	  //out_file=data_inf.out_file+"gal_shape.dat";
	  //outp_gal(gshape,data_inf,out_file);
	  data_inf.do_SS=0;

	  delete []gshape;
	}
      data_inf.do_DD=1;
      //out_file=data_inf.out_file+"gal_density.dat";
      //outp_gal(gdensity,data_inf,out_file);
      data_inf.do_DD=0;
      delete []gdensity;
      if (!data_inf.RsRd_same && !data_inf.periodic_box)
	delete []grand_density;
    }
  return 0;
}

int do_bining(bins &B,data_info &data_inf)
{
  cout<<"Bin calculations"<<endl;

  string of="";
  of=data_inf.out_file+"bins_auto_SR.temp";
  outp_2Dbins(B.SRs.b,data_inf,of);

  of=data_inf.out_file+"bins_cross_SR.temp";
  outp_2Dbins(B.SRd.b,data_inf,of);

  omp_set_dynamic(0); //Open Mp threads
  if(data_inf.n_threads==0||data_inf.n_threads>omp_get_num_procs()) //no input (or 0) from user about number of threads. Use all available cores
    omp_set_num_threads(omp_get_num_procs());
  else
    omp_set_num_threads(data_inf.n_threads); //Number of threads input by user

  if (data_inf.do_auto_corr)
  {
    if (data_inf.do_auto_rebin)
	   {
	      re_bin(B.SS.b,data_inf);
        re_bin(B.SRs.b,data_inf);
	      re_bin(B.RsRs.b,data_inf);

	      of=data_inf.out_file+"bins_auto_SR2.temp";
	      outp_2Dbins(B.SRs.b,data_inf,of);
	   }
    final_calc_auto(B.SS.b,B.SRs.b,B.RsRs.b,B.final_auto.b,data_inf); // calculations for p bins.

    //int_calc(B.SS.b,B.RsRs.b,B.SRs.b,B.final_auto.b,data_inf); // integration over p bins
    string out_file="";
    out_file=data_inf.out_file+"bins_auto_final.dat";
    outp_bins(B.final_auto.b,data_inf,out_file);

    out_file=data_inf.out_file+"bins_auto_SS.dat";
    outp_bins(B.SS.b,data_inf,out_file);

    out_file=data_inf.out_file+"bins_auto_SR.dat";
    outp_bins(B.SRs.b,data_inf,out_file);

    out_file=data_inf.out_file+"bins_auto_RR.dat";
    outp_bins(B.RsRs.b,data_inf,out_file);

    out_file=data_inf.out_file+"bins2D_auto_final.dat";
    outp_2Dbins(B.final_auto.b,data_inf,out_file);

    out_file=data_inf.out_file+"bins2D_auto_SS.dat";
    outp_2Dbins(B.SS.b,data_inf,out_file);

    out_file=data_inf.out_file+"bins2D_auto_SR.dat";
    outp_2Dbins(B.SRs.b,data_inf,out_file);

    out_file=data_inf.out_file+"bins2D_auto_RR.dat";
    outp_2Dbins(B.RsRs.b,data_inf,out_file);

    if(data_inf.do_jk)
      {
        #pragma omp parallel
        {
        #pragma omp for schedule(dynamic,10)
        for(int i=0;i<data_inf.n_jk;i++)
          {
            if(data_inf.do_auto_rebin)
              {
                re_bin(B.SS.jk[i].b,data_inf);
                re_bin(B.RsRs.jk[i].b,data_inf);
                re_bin(B.SRs.jk[i].b,data_inf);
              }

              final_calc_auto(B.SS.jk[i].b,B.SRs.jk[i].b,B.RsRs.jk[i].b,B.final_auto.jk[i].b,
                                                      data_inf); // calculations for p bins.
              //int_calc(B.SS.jk[i].b,B.RsRs.jk[i].b,B.SRs.jk[i].b,B.final_auto.jk[i].b,data_inf); // integration over p bins

	            std::ostringstream s;
              if (data_inf.do_jk==1)
                s<<i;
              else{
                int ij[2]={0,0};
                indx_jk2(data_inf,ij,i);
                s<<ij[0]<<"_"<<ij[1]<<"_"<<i;
              }

              out_file=data_inf.out_file+"bins_auto_jk"+s.str()+"_SS.dat";
              outp_bins(B.SS.jk[i].b,data_inf,out_file);
              out_file=data_inf.out_file+"bins_auto_jk"+s.str()+"_RR.dat";
              outp_bins(B.RsRs.jk[i].b,data_inf,out_file);
              out_file=data_inf.out_file+"bins_auto_jk"+s.str()+"_SR.dat";
              outp_bins(B.SRs.jk[i].b,data_inf,out_file);
              out_file=data_inf.out_file+"bins_auto_jk"+s.str()+"_final.dat";
              outp_bins(B.final_auto.jk[i].b,data_inf,out_file);

	            out_file=data_inf.out_file+"bins2D_auto_jk"+s.str()+"_final.dat";
              outp_2Dbins(B.final_auto.jk[i].b,data_inf,out_file);
              out_file=data_inf.out_file+"bins2D_auto_jk"+s.str()+"_SS.dat";
              outp_2Dbins(B.SS.jk[i].b,data_inf,out_file);
              out_file=data_inf.out_file+"bins2D_auto_jk"+s.str()+"_RR.dat";
              outp_2Dbins(B.RsRs.jk[i].b,data_inf,out_file);
              out_file=data_inf.out_file+"bins2D_auto_jk"+s.str()+"_SR.dat";
              outp_2Dbins(B.SRs.jk[i].b,data_inf,out_file);

            }}

	      bin bins_final_jk[data_inf.n_bins];

	      if(data_inf.lin_bin) // linear bins intialisation
	       bin_lin(bins_final_jk,data_inf);
	      else
	       bin_log(bins_final_jk,data_inf);

        if (data_inf.do_jk==1){
	         jk_final(bins_final_jk,B.final_auto.b,B.final_auto.jk,data_inf);
           jk_final2D(bins_final_jk,B.final_auto.b,B.final_auto.jk,data_inf);

          out_file=data_inf.out_file+"bins_auto_jk_final.dat";
          outp_bins(bins_final_jk,data_inf,out_file);

          out_file=data_inf.out_file+"bins2D_auto_jk_final.dat";
          outp_2Dbins(bins_final_jk,data_inf,out_file);
        }
      }
  }

  if (data_inf.do_cross_corr)
    {
      if (data_inf.estimator==4)
	     final_calc_mean_projected(B.SD.b,B.SRd.b,B.final_cross.b,data_inf);
      else
	     final_calc(B.SD.b,B.SRd.b,B.DRs.b,B.RsRd.b,B.final_cross.b,data_inf); // calculations for p bins.

      //int_calc(B.SD.b,B.SRd.b,B.DRs.b,B.RsRd.b,B.final_cross.b,data_inf); // integration over p bins

      string out_file="";
      out_file=data_inf.out_file+"bins_cross_final.dat";
      outp_bins(B.final_cross.b,data_inf,out_file);

      out_file=data_inf.out_file+"bins2D_cross_final.dat";
      outp_2Dbins(B.final_cross.b,data_inf,out_file);

      out_file=data_inf.out_file+"bins_cross_SD.dat";
      outp_bins(B.SD.b,data_inf,out_file);

      out_file=data_inf.out_file+"bins2D_cross_SD.dat";
      outp_2Dbins(B.SD.b,data_inf,out_file);

      out_file=data_inf.out_file+"bins_cross_SR.dat";
      outp_bins(B.SRd.b,data_inf,out_file);

      out_file=data_inf.out_file+"bins2D_cross_SR.dat";
      outp_2Dbins(B.SRd.b,data_inf,out_file);

      out_file=data_inf.out_file+"bins_cross_DR.dat";
      outp_bins(B.DRs.b,data_inf,out_file);

      out_file=data_inf.out_file+"bins2D_cross_DR.dat";
      outp_2Dbins(B.DRs.b,data_inf,out_file);

      out_file=data_inf.out_file+"bins_cross_RR.dat";
      outp_bins(B.RsRd.b,data_inf,out_file);

      out_file=data_inf.out_file+"bins2D_cross_RR.dat";
      outp_2Dbins(B.RsRd.b,data_inf,out_file);

      if(data_inf.do_jk)
	     {
         //#pragma omp parallel
         {
         #pragma omp for schedule(dynamic,10)
          for(int i=0;i<data_inf.n_jk;i++)
            {
	             if (data_inf.estimator==4)
		             final_calc_mean_projected(B.SD.jk[i].b,B.SRd.jk[i].b,B.final_cross.jk[i].b,
                                            data_inf);
	             else
		             final_calc(B.SD.jk[i].b,B.SRd.jk[i].b,B.DRs.jk[i].b,B.RsRd.jk[i].b,
                            B.final_cross.jk[i].b,data_inf); // calculations for p bins.

	            std::ostringstream s;
              if (data_inf.do_jk==1)
                s<<i;
              else{
                int ij[2]={0,0};
                indx_jk2(data_inf,ij,i);
                s<<ij[0]<<"_"<<ij[1]<<"_"<<i;
              }
              cout<<"Bin calculations:jk "<<i<<s.str()<<endl;

              out_file=data_inf.out_file+"bins_cross_jk"+s.str()+"_SD.dat";
              outp_bins(B.SD.jk[i].b,data_inf,out_file);
              out_file=data_inf.out_file+"bins_cross_jk"+s.str()+"_RR.dat";
              outp_bins(B.RsRd.jk[i].b,data_inf,out_file);
              out_file=data_inf.out_file+"bins_cross_jk"+s.str()+"_DR.dat";
              outp_bins(B.DRs.jk[i].b,data_inf,out_file);
	            out_file=data_inf.out_file+"bins_cross_jk"+s.str()+"_SR.dat";
              outp_bins(B.SRd.jk[i].b,data_inf,out_file);
	            out_file=data_inf.out_file+"bins_cross_jk"+s.str()+"_final.dat";
              outp_bins(B.final_cross.jk[i].b,data_inf,out_file);

              out_file=data_inf.out_file+"bins2D_cross_jk"+s.str()+"_final.dat";
              outp_2Dbins(B.final_cross.jk[i].b,data_inf,out_file);
              out_file=data_inf.out_file+"bins2D_cross_jk"+s.str()+"_SD.dat";
              outp_2Dbins(B.SD.jk[i].b,data_inf,out_file);
              out_file=data_inf.out_file+"bins2D_cross_jk"+s.str()+"_RR.dat";
              outp_2Dbins(B.RsRd.jk[i].b,data_inf,out_file);
              out_file=data_inf.out_file+"bins2D_cross_jk"+s.str()+"_DR.dat";
              outp_2Dbins(B.DRs.jk[i].b,data_inf,out_file);
              out_file=data_inf.out_file+"bins2D_cross_jk"+s.str()+"_SR.dat";
              outp_2Dbins(B.SRd.jk[i].b,data_inf,out_file);
            }}
	        bin bins_final_jk[data_inf.n_bins];

          if(data_inf.lin_bin) // linear bins intialisation
            bin_lin(bins_final_jk,data_inf);
          else
            bin_log(bins_final_jk,data_inf);

          if (data_inf.do_jk==1){
            jk_final(bins_final_jk,B.final_cross.b,B.final_cross.jk,data_inf);

            jk_final2D(bins_final_jk,B.final_cross.b,B.final_cross.jk,data_inf);

            out_file=data_inf.out_file+"bins_cross_jk_final.dat";
            outp_bins(bins_final_jk,data_inf,out_file);

  	        out_file=data_inf.out_file+"bins2D_cross_jk_final.dat";
            outp_2Dbins(bins_final_jk,data_inf,out_file);
          }
	      }
    }
  return 0;
}
