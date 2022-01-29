// intialize data with random values and bins with values 0 ;

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
#include "sky_calcs.h"
#include "PB_calcs.h"
#include "corels.h"
#include "jk.h"
using namespace std;


//==================================================================================================================

int da_Z(data_info &data_inf) // read and store angular diameter distance as array
{
  double z_temp=0,z=0;
  ifstream data;
  data.open(data_inf.distance_file.c_str());
  string junk; // to discard unnecessary columns

  if (data.is_open())
    {
      int i=0;

      while (!data.eof()&&(i<(ceil((data_inf.z_max+1)/data_inf.dz)+1)))
	{

	  if(!data_inf.use_comoving)
	    {
	      data>>z>>data_inf.H_z[i]>>data_inf.dA_z[i]>>data_inf.dC_z[i]>>z_temp;  //values from files.. angular diameter distance
	    }
	  else
	    {
	      data>>z>>data_inf.H_z[i]>>z_temp>>data_inf.dC_z[i]>>data_inf.dA_z[i];  //values from files.. transverse comoving distance
	    }
	  i++;
	}

      if(i<ceil(data_inf.z_max/data_inf.dz)+1)
	{
	  cout<<"da_Z:: distances file (dA vs Z) ended unexpectedly"<<endl;
	  return 1;
	}
    }
  else
    {
      cout<<"da_Z:: distances file (dA vs Z) not opened: "<<data_inf.distance_file<<endl;return 1;
    }
  data.close();
  return 0;
}

int data_ini(data_info &data_inf,string FileName)//read input file and intialise main data_info struct
{
  ifstream data;
  data.open(FileName.c_str()); // input file from user
  string junk;// discard first column of inp file

  if (data.is_open())
    { // read data
      data>>junk>>data_inf.which_corr;
      data>>junk>>data_inf.coordinates;  // new line in ver3
      data>>junk>>data_inf.estimator;
      data>>junk>>data_inf.data_sorted;
      data>>junk>>data_inf.use_comoving;
      data>>junk>>data_inf.do_jk;
      data>>junk>>data_inf.sig_crit; //ver3

      data>>junk>>data_inf.shape_pos_file;
      data>>junk>>data_inf.shape_z_file;
      data>>junk>>data_inf.shape_e_file;
      data>>junk>>data_inf.shape_wt_file;
      data>>junk>>data_inf.shape_jk_file;
      data>>junk>>data_inf.shape_patch_file; // ver3

      data>>junk>>data_inf.density_pos_file;
      data>>junk>>data_inf.density_z_file;
      data>>junk>>data_inf.density_wt_file;
      data>>junk>>data_inf.density_jk_file;
      data>>junk>>data_inf.density_e_file;

      data>>junk>>data_inf.shape_rand_pos_file;
      data>>junk>>data_inf.shape_rand_z_file;
      data>>junk>>data_inf.shape_rand_wt_file;
      data>>junk>>data_inf.shape_rand_jk_file;
      data>>junk>>data_inf.shape_rand_patch_file; //ver3

      data>>junk>>data_inf.density_rand_pos_file;
      data>>junk>>data_inf.density_rand_z_file;
      data>>junk>>data_inf.density_rand_wt_file;
      data>>junk>>data_inf.density_rand_jk_file;

      data>>junk>>data_inf.distance_file;
      data>>junk>>data_inf.patches_file;

      data>>junk>>data_inf.out_file;

      data>>junk>>data_inf.n_threads;

      data>>junk>>data_inf.n_shape;
      data>>junk>>data_inf.n_density;
      data>>junk>>data_inf.n_shape_rand;
      data>>junk>>data_inf.n_density_rand;
      data>>junk>>data_inf.rand_subsample;//ver3
      data>>junk>>data_inf.n_jk_regions;
      data>>junk>>data_inf.n_patch; //ver3

      data>>junk>>data_inf.bin_r_min;
      data>>junk>>data_inf.bin_r_max;
      data>>junk>>data_inf.n_bins;
      data>>junk>>data_inf.lin_bin;
      data>>junk>>data_inf.n_p_bin;
      data>>junk>>data_inf.p_min;
      data>>junk>>data_inf.p_max;
      data>>junk>>data_inf.z_min;
      data>>junk>>data_inf.z_max;
      data>>junk>>data_inf.dz;
      data>>junk>>data_inf.z_sep_min; //ver3
      data>>junk>>data_inf.z_sep_max; //ver3
      data>>junk>>data_inf.periodic_box;
      data>>junk>>data_inf.periodic_box_size;
    }

  else
    {
      cout<<"data_ini:: initialisation File not open"<<FileName<<endl;return 1;
    }
  data.clear();
  data.close();

  data_inf.CMB_DIST=9820.428;
  data_inf.n_bin_var=7;

  if (data_inf.periodic_box!=0)
    {
      seed(); //seeding random number generator.. function defined in corels.cpp
      //data_inf.xyz_coord=1;
      if (data_inf.coordinates!=6&&data_inf.coordinates!=7)
	     {
	        cout<<" wrong coordinates choice for periodic box:"<< data_inf.periodic_box_size<<" "<<data_inf.coordinates<<endl;
	         return 1;
	     }
      data_inf.dz=0;
      cout<<"doing xyz box: size="<<data_inf.periodic_box_size<<"   periodic="<<data_inf.periodic_box<<" njk: "<<data_inf.n_jk_regions<<endl;
      if (data_inf.bin_r_max>data_inf.periodic_box_size||data_inf.p_max>data_inf.periodic_box_size)
	{
	  cout<<"bin limits > the periodic box size"<<endl;
	}
      if (data_inf.data_sorted==2){
	cout<<"data_sorted==2 doesnot work with periodic box.. use 1"<<endl;
	return 1;}
    }

  switch (data_inf.estimator)
    {
    case 0:
      data_inf.do_cross_corr=1;
      data_inf.do_auto_corr=0;
      break;
    case 2:
      data_inf.do_auto_rebin=1;
    case 1:
      data_inf.do_auto_corr=1;
      data_inf.do_cross_corr=0;
      break;
    case 3:
      data_inf.do_auto_corr=1;
      data_inf.do_cross_corr=1;
      data_inf.do_auto_rebin=1;
      break;
    default:
      data_inf.do_cross_corr=1;
      data_inf.do_auto_corr=0;
      data_inf.do_auto_rebin=0;
      break;
    }

  data_inf.RsRd_same=0;

  if(data_inf.n_density_rand==0&&data_inf.which_corr<=2)
    {
      data_inf.RsRd_same=1;
      data_inf.n_density_rand=data_inf.n_shape_rand;
    }
  data_inf.include_prob_RR=1;
  if (data_inf.rand_subsample)
    data_inf.include_prob_RR=sqrt(sqrt(double(data_inf.n_density)/data_inf.n_density_rand*double(data_inf.n_shape)/data_inf.n_shape_rand))*1.0;

  cout<<"RsRd same  "<<data_inf.RsRd_same<<"  RR-prob:"<<data_inf.include_prob_RR<<endl;
  //if(data_inf.do_jk!=0){
    jk_initilize(data_inf);
    //}

  if (data_inf.which_corr>=7)
    {
      if (data_inf.z_sep_max==0)
	     data_inf.z_sep_max=10;
      data_inf.n_density=0;
      data_inf.patches=new int [data_inf.n_patch];
      ifstream data_patches;
      data_patches.open(data_inf.patches_file.c_str());
      for (int i=0;i<data_inf.n_patch;i++)
	     {
	        data_patches>>data_inf.patches[i];
	         if (!data_patches.good())
      	    {
      	      cout<<"patches file problem: "<<data_inf.patches_file<<endl;
      	      return 1;
      	    }
	     }
      data_patches.close();
    }

  int z_idx=ceil(data_inf.z_max/data_inf.dz);

  /*  if (data_inf.p_max==0 && data_inf.p_min==0)
    {
      cout<<"Will work in r-mu plane"<<endl;
      data_inf.do_r_mu=1;
      data_inf.p_max=1;
      data_inf.p_min=-1;
    }
  */

  if (data_inf.z_sep_max==0)
    {
      data_inf.z_sep_max=9;
    }
  //data_inf.z_sep_min=0;
  data_inf.PB_sorted_starts_S=0;
  data_inf.PB_sorted_starts_D=0;

  if (data_inf.data_sorted==0)
    {
      data_inf.PB_sorted_starts_S=0;
      data_inf.PB_sorted_starts_D=0;
    }

  //  if (!data_inf.xyz_coord){
  if (data_inf.dz==0 && data_inf.periodic_box==0)
    {
      cout<<"dz=0"<<endl;
      return 1;
    }
    int k_n=1;
    if (data_inf.periodic_box==0){
        k_n=ceil((data_inf.z_max+1)/data_inf.dz)+20;
    }
    const unsigned int k=k_n;
//20 is some buffer since we use higher values sometime for interpolation

  data_inf.dA_z=new double[k];// make array for angular diameter distance
  data_inf.dC_z=new double[k];
  data_inf.H_z=new double[k];// make array for hubble parameter

  int ch=0;
  if (data_inf.periodic_box==0)
      ch=da_Z(data_inf);	//reading da_z values from file.. look in initialization


  if(ch==1)return 1;

  if (data_inf.coordinates==0)
    {
      data_inf.z_sep_max=data_inf.p_max*data_inf.H_z[z_idx]/300000.0+0.001;
      data_inf.z_sep_min=data_inf.p_min*data_inf.H_z[z_idx]/300000.0-0.001;
    }
  else if (data_inf.coordinates==1)
    {
      data_inf.z_sep_max=data_inf.bin_r_max*data_inf.H_z[z_idx]/300000.0+0.001;
      data_inf.z_sep_min=-1*data_inf.z_sep_max;
    }
  else if (data_inf.coordinates==6) //sims, xyz, rp-pi
  {
    data_inf.z_sep_max=data_inf.p_max+1;
    data_inf.z_sep_min=data_inf.p_min-1;
  }
  else if (data_inf.coordinates==7) //sims, xyz, rp-pi
  {
    data_inf.z_sep_max=data_inf.bin_r_max+0.1;
    data_inf.z_sep_min=data_inf.z_sep_max*-1;
  }
  /*  else
    {
      data_inf.z_sep_max=10;
      data_inf.z_sep_min=-10;
      }*/
    //}

  if(data_inf.bin_r_max==0||data_inf.bin_r_min>data_inf.bin_r_max)
    {
      cout<<"Bad Bin limits"<<endl;
      return 1;
    }

  data_inf.p_size=(data_inf.p_max-data_inf.p_min)/(data_inf.n_p_bin);

  if (data_inf.coordinates==3||data_inf.coordinates==5)
    data_inf.p_size=1;

  if(data_inf.lin_bin)data_inf.r_size=(data_inf.bin_r_max-data_inf.bin_r_min)/data_inf.n_bins;
  else data_inf.r_size=(log10(data_inf.bin_r_max)-log10(data_inf.bin_r_min))/data_inf.n_bins;

  data_inf.g_with_shape=0;
  data_inf.do_DD=0;
  data_inf.do_DRs=0;
  data_inf.do_RsRs=0;

  data_inf.do_SS=0;
  data_inf.do_SRs=0;
  data_inf.do_SD=0;

  data_inf.do_SRd=0;
  data_inf.do_RsRd=0;

  data_inf.SS_wt.val=0;//double(data_inf.n_shape)*(double(data_inf.n_shape)-1.0)/2.0;
  data_inf.RsRs_wt.val=0;//double(data_inf.n_shape_rand)*(double(data_inf.n_shape_rand)-1.0)/2.0;
  data_inf.RsRd_wt.val=0;//double(data_inf.n_shape_rand)*double(data_inf.n_density_rand);
  data_inf.DD_wt.val=0;//double(data_inf.n_density)*(double(data_inf.n_density)-1.0)/2.0;

  data_inf.SRs_wt.val=0;//double(data_inf.n_shape)*double(data_inf.n_shape_rand);
  data_inf.SRd_wt.val=0;//double(data_inf.n_shape)*double(data_inf.n_density_rand);
  data_inf.SD_wt.val=0;//double(data_inf.n_shape)*double(data_inf.n_density);
  data_inf.DRs_wt.val=0;//double(data_inf.n_density)*double(data_inf.n_shape_rand);

  data_inf.DD_wt.err=data_inf.DD_wt.val;
  data_inf.RsRs_wt.err=data_inf.RsRs_wt.val;
  data_inf.SRd_wt.err=data_inf.SRd_wt.val;
  data_inf.SS_wt.err=data_inf.SS_wt.val;
  data_inf.DRs_wt.err=data_inf.DRs_wt.val;
  data_inf.RsRd_wt.err=data_inf.RsRd_wt.val;
  data_inf.SRs_wt.err=data_inf.SRs_wt.val;
  data_inf.SD_wt.err=data_inf.SD_wt.val;
  data_inf.S_wt.val=0;
  data_inf.D_wt.val=0;
  data_inf.Rs_wt.val=0;
  data_inf.Rd_wt.val=0;
  data_inf.S_wt.err=0;
  data_inf.D_wt.err=0;
  data_inf.Rs_wt.err=0;
  data_inf.Rd_wt.err=0;
  return 0;
}


//======================bins======================================================

void P_bin_ini(P_bin p[],data_info &data_inf) // intialize bins
{
  float p_size=data_inf.p_size; // bin size
  float temp=data_inf.p_min;
  for(int i=0;i<data_inf.n_p_bin;i++)
    {
      p[i].p_min=temp;// lower limit of the bin
      temp+=p_size;
      p[i].p_max=temp;// upper limit of bin

      p[i].ls_err=0.0; // other values assinged 0
      //p[i].num_den=0;
      p[i].num=0;
      p[i].wt_num=0;
      //p[i].wt_num.err=0;

      for(int j=0;j<7;j++)
	{
	  p[i].data[j].val=0;
	  p[i].data[j].err=0;
	}
    }
}

void bin_lin(bin bi[],data_info &data_inf)	//making linear bins
{
  float b_size=(data_inf.bin_r_max-data_inf.bin_r_min)/(float)data_inf.n_bins; // bin size
  float temp=data_inf.bin_r_min;

  for(int i=0;i<data_inf.n_bins;i++)
    {
      bi[i].b_min=temp; // limits of the bin
      temp+=b_size;
      bi[i].b_max=temp;

      bi[i].wt_num=0;
    //  bi[i].wt_num.err=0;
      bi[i].ls_err_int=0; // other values 0
      // bi[i].num_den=0;
      bi[i].num=0;
      // bi[i].b_err=0;
      // bi[i].num_den=0;
      bi[i].sig_crit=0;
      const unsigned int kn=data_inf.n_p_bin;

      bi[i].p_bin=new P_bin[kn]; // p bins
      P_bin_ini(bi[i].p_bin,data_inf);
      for(int j=0;j<7;j++)
	     {
	        bi[i].data[j].val=0;
	        bi[i].data[j].err=0;
	      }
    }
}

void bin_log(bin bi[],data_info &data_inf) // log bins
{
  float d_min_log=log10(data_inf.bin_r_min);
  float d_max_log=log10(data_inf.bin_r_max);

  float b_size=(d_max_log-d_min_log)/(float)data_inf.n_bins; // bin size

  for(int i=0;i<data_inf.n_bins;i++)
    {
      bi[i].ls_err_int=0;
      bi[i].b_min=pow(10.0,d_min_log);
      d_min_log+=b_size;
      bi[i].b_max=pow(10.0,d_min_log);
      bi[i].num=0;
      bi[i].wt_num=0;
      //bi[i].wt_num.err=0;
      //bi[i].num_den=0;
      //bi[i].b_err=0;
      bi[i].sig_crit=0;

      const unsigned int kn=data_inf.n_p_bin;

      bi[i].p_bin=new P_bin[kn];

      P_bin_ini(bi[i].p_bin,data_inf);

      for(int j=0;j<7;j++)
	     {
	        bi[i].data[j].val=0;
	        bi[i].data[j].err=0;
	      }
    }
}


void jk_ini(jackknife jk[], data_info &data_inf)
{
  for (int i=0;i<data_inf.n_jk;i++)
    {
      //jk[i].region=data_inf.jk_regions[i];
      jk[i].b=new bin[data_inf.n_bins];

      if (data_inf.lin_bin)
	     bin_lin(jk[i].b,data_inf);
      else bin_log(jk[i].b,data_inf);
    }
}

void bin_jk_ini(bin_jk &bjk,data_info &data_inf){
  bjk.b=new bin [data_inf.n_bins];
  bjk.jk=new jackknife[data_inf.n_jk];
  //  void (*bin_ini)(bin b[],data_info &data_inf);
  if(data_inf.lin_bin) // linear bins intialisation
    bin_lin(bjk.b,data_inf);
  else
    bin_log(bjk.b,data_inf);
  //bin_ini(bjk.b,data_inf);
  jk_ini(bjk.jk,data_inf);
}

void bins_ini(bins &B,data_info &data_inf)
{
  if (data_inf.do_auto_corr){
    bin_jk_ini(B.SS,data_inf);
    bin_jk_ini(B.final_auto_jk,data_inf);
    bin_jk_ini(B.final_auto,data_inf);
  }

  bin_jk_ini(B.SD,data_inf);
  bin_jk_ini(B.SRs,data_inf);
  bin_jk_ini(B.SRd,data_inf);
  bin_jk_ini(B.DRs,data_inf);
  bin_jk_ini(B.RsRd,data_inf);
  bin_jk_ini(B.RsRs,data_inf);
  bin_jk_ini(B.final_cross_jk,data_inf);
  bin_jk_ini(B.final_cross,data_inf);
}

void choose_corel_func(calc_temp &ct,data_info &data_inf)
{
  if ((data_inf.do_RsRs==1||data_inf.do_RsRd==1))
    //&&data_inf.which_corr<7&&(data_inf.coordinates!=2&&data_inf.coordinates!=4))
    {
      ct.corel=&density_density_corel;
      return;
    }

  switch(data_inf.which_corr)
    {
    case 0:
      ct.corel=&density_density_corel;
      break;
    case 1:
      ct.corel=&shape_density_corel;
      if(data_inf.do_SS==1)
	      ct.corel=&shape_shape_corel;
      break;
    case 2:
      ct.corel=&shape_shape_corel;
        if(data_inf.do_SD!=1 && data_inf.do_SS!=1)
	        ct.corel=&shape_density_corel;
        if(data_inf.coordinates==2||data_inf.coordinates==4)
	      {
        if(data_inf.do_DRs==1||data_inf.do_SD==1)
          ct.corel=&shape_shape_corel;
        else
          ct.corel=&density_density_corel;
        }
      break;
    case 3:
      ct.corel=&density_convergence_corel;
      break;
    case 4:
      ct.corel=&shape_convergence_corel;
      break;
    case 5:
      ct.corel=&ED_corel;
      break;
    case 7:
      ct.corel=&shape_density_corel;
      break;
    case 8:
      ct.corel=&shape_shape_corel;
      break;
    case 9:
      ct.corel=&shape_convergence_corel;
      break;
    case 10:
      ct.corel=&density_density_corel;
      break;
    }
}

void choose_rp_pi_func(calc_temp &ct,data_info &data_inf)
{
  ct.rp_calc=&rp_calc_sky;
  ct.pi_calc=&pi_calc;

  if (data_inf.p_max>600||data_inf.coordinates==3||data_inf.coordinates==5)
    {
      ct.pi_calc=&pi_calc_projected;
      ct.rp_calc=&rp_calc_sky_projected; //not for coordinates=5. Changed below
    }

  if (data_inf.periodic_box || data_inf.coordinates==6||data_inf.coordinates==7)
    {
      ct.rp_calc=&rp_calc_PB;
      ct.pi_calc=&pi_calc_PB;
    }

  if (data_inf.coordinates==2||data_inf.coordinates==4)
    ct.pi_calc=&pi_calc_phi;
  if (data_inf.coordinates==4||data_inf.coordinates==5)
    ct.rp_calc=&rp_calc_sky_angular;

  ct.r_calc=ct.rp_calc;

  if(data_inf.coordinates==1||data_inf.coordinates==7)
    {
      ct.r_calc=&r_mu_calc;
    }
}

void calc_temp_ini(calc_temp &ct,data_info &data_inf)
{
  bin_jk_ini(ct.bjk,data_inf);
  ct.n1=-10;ct.n2=-10;
  ct.n_patch_gal=0;
  ct.e_new_i[0]=0;ct.e_new_j[0]=0;
  ct.e_new_i[1]=0;ct.e_new_j[1]=0;
  ct.temp[1]=0;ct.temp[1]=0;
  ct.gjk[1]=0;ct.gjk[1]=0;
  ct.jk_prob[1]=1;ct.jk_prob[1]=1;
  ct.delta_z=0;ct.da=0;ct.z_avg=0;ct.d_P=0;ct.da_R=0;ct.wt_f.val=0;ct.wt_f.err=0;ct.theta=0;ct.z_max=0;ct.z_min=0;
  ct.ep.val=0;ct.ep.err=0;
  ct.ex.val=0;ct.ex.err=0;
  ct.epp.val=0;ct.epp.err=0;
  ct.exx.val=0;ct.exx.err=0;
  ct.epx.val=0;ct.epx.err=0;
  ct.etheta.val=0;ct.etheta.err=0;
  ct.dec_max=0;ct.dec_min=0;
  ct.ra_max=0;ct.ra_min=0;ct.ra_lim=0;
  ct.PBx=0;ct.PBy=0;ct.PBz=0;ct.fail=0;
  ct.da_max=0;
  ct.sig_crit=1;ct.sig_crit_inv=0;ct.dl=0;ct.ds=0;ct.dls=0; //sig_crit is used by default but not updated.. it should be 1 here
  ct.include_prob_max=1;
  if (data_inf.do_RsRs or data_inf.do_RsRd)
    ct.include_prob_max=data_inf.include_prob_RR;
  choose_corel_func(ct,data_inf);
  choose_rp_pi_func(ct,data_inf);

  ct.jk_check=&jk_check_sky;
  ct.angle_orient_axis=&angle_orient_axis_sky;
  /*  if (data_inf.which_corr>=7)
    ct.jk_check=&jk_check_sky_wl;
  */

  ct.ret_z=1;
  ct.ret_dec=1;
  if (data_inf.data_sorted==1)
    ct.ret_z=1.e9;
  else if (data_inf.data_sorted==2)
    ct.ret_dec=1.e9;

  if (data_inf.periodic_box || data_inf.coordinates==6)
    {
      ct.jk_check=&jk_check_PB;
      ct.angle_orient_axis=&angle_orient_axis_PB;
    }
  if (data_inf.do_jk==0||data_inf.n_jk==0)
    {
      ct.jk_check=&no_jk;
      ct.jk_prob[1]=0;ct.jk_prob[1]=0;
    }
  ct.sig_crit_func=sig_crit_calc;
  if (data_inf.sig_crit>1)
    ct.sig_crit_func=sig_crit_g1;
}
