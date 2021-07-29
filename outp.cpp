#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <math.h>
#include <sstream>
#include "data_def.h"
#include <iomanip>
#include <limits>

int outp_gal(gal g[],data_info data_inf,string FileName)	//output galaxy data set to file
{
  ofstream writ1;
  writ1.open(FileName.c_str()); // open file

  int n_gal=0;

  if (data_inf.do_SS)n_gal=data_inf.n_shape;
  else if (data_inf.do_RsRs)n_gal=data_inf.n_shape_rand;
  else if (data_inf.do_RsRd)n_gal=data_inf.n_density_rand;
  else if (data_inf.do_DD)n_gal=data_inf.n_density;

  if (!writ1.is_open())
    {
      cout<<FileName<<" not open"<<endl;
      return 1;
    }

  writ1<<"# RA        DEC         Z        e1       e2     wt    jkr    phi  lens_wt"<<endl;

  for (int i=0;i<n_gal;i++)// write galaxies out
    {
      writ1<<g[i].ra.val_deg<<"  "<<g[i].dec.val_deg<<"		"<<g[i].redshift;

      writ1<<"  	"<<g[i].e[0]<<"		"<<g[i].e[1]<<"        "<<g[i].wt<<"      "<<g[i].jk<<"    "<<g[i].phi<<"   "<<g[i].lens_wt<<endl;

    }

  writ1.close();
  return 0;
}


int outp_bins(bin b[],data_info data_inf, string FileName)	//output bins to file
{
	ofstream writ2;
	writ2.open(FileName.c_str());
	//	writ2<<std::setprecision(std::numeric_limits<long double>::digits10);
	if (!writ2.is_open())
	  {
	    cout<<FileName<<" not open"<<endl;
	    return 1;
	  }

//column headings
 writ2<<"# rp	npairs	wgg	wgg_err	wgp	wgp_err	wgx wgx_err	wpp 	wpp_err		wxx	wxx_err		wpx  wpx_err 	theta   theta_err   LS_err     wt_npair    rp_low   rp_high sig_crit"<<endl;

//output data
 for (int i=0;i<data_inf.n_bins;i++)
   {
     writ2<<(b[i].b_max+b[i].b_min)/2<<"	"<<b[i].num;//<<"	"<<b[i].b_err;

     for(int j=0;j<7;j++)
       {
	        writ2<<"	"<<b[i].data[j].val<<"	"<<b[i].data[j].err;
       }

     writ2<<"	"<<b[i].ls_err_int<<"   "<<b[i].wt_num<<"      "<<b[i].b_min<<"    "<<b[i].b_max<<"    "<<b[i].sig_crit<<endl;
   }

 writ2.close();
 return 0;
}


int outp_data_info(data_info data_inf,string FileName)
{
  ofstream data;
  data.open(FileName.c_str());

  if (!data.is_open())
    {
      cout<<FileName<<" not open"<<endl;
      return 1;
    }

  else{
    data<<"do_auto_corr         "<<data_inf.do_auto_corr<<endl;
    data<<"do_cross_corr        "<<data_inf.do_cross_corr<<endl;
    data<<"data_sorted          "<<data_inf.data_sorted<<endl;
    data<<"use_comoving         "<<data_inf.use_comoving<<endl;
    data<<"do_auto_rebin        "<<data_inf.do_auto_rebin<<endl;
    data<<"do_jk                "<<data_inf.do_jk<<endl;

    data<<"shape_pos            "<<data_inf.shape_pos_file<<endl;
    data<<"shape_z              "<<data_inf.shape_z_file<<endl;
    data<<"shape_e              "<<data_inf.shape_e_file<<endl;
    data<<"shape_wt             "<<data_inf.shape_wt_file<<endl;
    data<<"shape_jk             "<<data_inf.shape_jk_file<<endl;

    data<<"density_pos          "<<data_inf.density_pos_file<<endl;
    data<<"density_z            "<<data_inf.density_z_file<<endl;
    data<<"density_wt           "<<data_inf.density_wt_file<<endl;
    data<<"density_jk           "<<data_inf.density_jk_file<<endl;
    data<<"density_e            "<<data_inf.density_e_file<<endl;

    data<<"shape_rand_pos       "<<data_inf.shape_rand_pos_file<<endl;
    data<<"shape_rand_z         "<<data_inf.shape_rand_z_file<<endl;
    data<<"shape_rand_wt        "<<data_inf.shape_rand_wt_file<<endl;
    data<<"shape_rand_jk        "<<data_inf.shape_rand_jk_file<<endl;

    data<<"density_rand_pos     "<<data_inf.density_rand_pos_file<<endl;
    data<<"density_rand_z       "<<data_inf.density_rand_z_file<<endl;
    data<<"density_rand_wt      "<<data_inf.density_rand_wt_file<<endl;
    data<<"density_rand_jk      "<<data_inf.density_rand_jk_file<<endl;

    data<<"distances            "<<data_inf.distance_file<<endl;

    data<<"output               "<<data_inf.out_file<<endl;

    data<<"n_threads            "<<data_inf.n_threads<<endl;

    data<<"n_shape              "<<data_inf.n_shape<<endl;
    data<<"n_density            "<<data_inf.n_density<<endl;

    data<<"n_shape_rand         "<<data_inf.n_shape_rand<<endl;
    data<<"n_density_rand       "<<data_inf.n_density_rand<<endl;
    data<<"n_jk                 "<<data_inf.n_jk<<endl;
    /*    data<<"ra_max               "<<data_inf.ra_max<<endl;
    data<<"ra_min               "<<data_inf.ra_min<<endl;
    data<<"dec_max              "<<data_inf.dec_max<<endl;
    data<<"dec_min              "<<data_inf.dec_min<<endl;*/
    data<<"bin_r_min            "<<data_inf.bin_r_min<<endl;
    data<<"bin_r_max            "<<data_inf.bin_r_max<<endl;
    data<<"n_bins               "<<data_inf.n_bins<<endl;
    data<<"lin_bin              "<<data_inf.lin_bin<<endl;
    data<<"n_p_bin              "<<data_inf.n_p_bin<<endl;
    data<<"p_min                "<<data_inf.p_min<<endl;
    data<<"p_max                "<<data_inf.p_max<<endl;
    data<<"z_min                "<<data_inf.z_min<<endl;
    data<<"z_max                "<<data_inf.z_max<<endl;
    data<<"dz                   "<<data_inf.dz<<endl;
    //    data<<"da_max               "<<data_inf.da_max.val_deg<<"    "<<data_inf.da_max.val_rad<<endl;
    data<<"p_bin_siz            "<<data_inf.p_size<<endl;
    data<<"z_sep_max            "<<data_inf.z_sep_max<<endl;
    data<<"z_sep_min            "<<data_inf.z_sep_min<<endl;
    data<<"DD_Wt                "<<data_inf.DD_wt.val<<endl;
    data<<"DRs_wt               "<<data_inf.DRs_wt.val<<endl;
    data<<"SS_wt                "<<data_inf.SS_wt.val<<endl;
    data<<"SRs_wt               "<<data_inf.SRs_wt.val<<endl;
    data<<"SD_wt                "<<data_inf.SD_wt.val<<endl;
    data<<"RsRs_Wt              "<<data_inf.RsRs_wt.val<<endl;
    data<<"RsRd_wt              "<<data_inf.RsRd_wt.val<<endl;
    data<<"S_wt                 "<<data_inf.S_wt.val<<endl;
    data<<"D_wt                 "<<data_inf.D_wt.val<<endl;
    data<<"Rs_wt                "<<data_inf.Rs_wt.val<<endl;
    data<<"Rd_wt                "<<data_inf.Rd_wt.val<<endl;
    data<<"RsRd_same            "<<data_inf.RsRd_same<<endl;
  }
  return 0;
}

int outp_temp_bins(bin b[],data_info data_inf,string FileName)
{

  ofstream writ2;
  writ2<<std::setprecision(std::numeric_limits<long double>::digits10);
  writ2.open(FileName.c_str());

  if (!writ2.is_open())
    {
      cout<<FileName<<" not open"<<endl;
      return 1;
    }

  for (int i=0;i<data_inf.n_bins;i++)
    {
      writ2<<b[i].num<<"    "<<b[i].wt_num<<"       "<<b[i].ls_err_int;

      for(int j=0;j<6;j++)
	{
	  writ2<<"      "<<b[i].data[j].val<<"        "<<b[i].data[j].err;
	}

      writ2<<endl;

      for(int k=0;k<data_inf.n_p_bin;k++)
	{
	  writ2<<b[i].p_bin[k].num<<"      "<<b[i].p_bin[k].wt_num<<"      "<<b[i].p_bin[k].ls_err;

	  for(int j=0;j<7;j++)
	    {
	      writ2<<"      "<<b[i].p_bin[k].data[j].val<<"        "<<b[i].p_bin[k].data[j].err;
	    }
	  writ2<<endl;
	}
    }
  writ2.close();
  return 0;
}

int inp_bins(bin b[],data_info data_inf, string FileName)     //read bins from file
{
  ifstream writ2;
  writ2.open(FileName.c_str());

  if (!writ2.is_open())
    {
      cout<<FileName<<" not open"<<endl;
      return 1;
    }

  else
    {
      for (int i=0;i<data_inf.n_bins;i++)
	{
	  writ2>>b[i].num>>b[i].wt_num>>b[i].ls_err_int;

	  for(int j=0;j<7;j++)
	    {
	      writ2>>b[i].data[j].val>>b[i].data[j].err;
	    }


	  for(int k=0;k<data_inf.n_p_bin;k++)
	    {
	      writ2>>b[i].p_bin[k].num>>b[i].p_bin[k].wt_num>>b[i].p_bin[k].ls_err;

	      for(int j=0;j<7;j++)
		{
		  writ2>>b[i].p_bin[k].data[j].val>>b[i].p_bin[k].data[j].err;
		}
	    }
	}
    }

  writ2.close();
  return 0;
}

int add_bins(bin b1[],bin b2[],bin final[], data_info data_inf, bool subtract)     //add two bins
{
  int b2_fact=1;
  if (subtract)
    b2_fact=-1;

  for (int i=0;i<data_inf.n_bins;i++)
    {
      final[i].num=b1[i].num+b2[i].num*b2_fact;

      final[i].wt_num=b1[i].wt_num+b2[i].wt_num*b2_fact;
      //final[i].wt_num.err=b1[i].wt_num.err+b2[i].wt_num.err*b2_fact;

      // final[i].num_den=b1[i].num_den+b2[i].num_den*b2_fact;

      // final[i].b_err=b1[i].b_err+b2[i].b_err*b2_fact;

      final[i].ls_err_int=b1[i].ls_err_int+b2[i].ls_err_int*b2_fact;

      final[i].sig_crit=b1[i].sig_crit+b2[i].sig_crit*b2_fact;

      for(int k=0;k<data_inf.n_p_bin;k++)
	{
	  final[i].p_bin[k].num=b1[i].p_bin[k].num+b2[i].p_bin[k].num*b2_fact;

	  final[i].p_bin[k].wt_num=b1[i].p_bin[k].wt_num+b2[i].p_bin[k].wt_num*b2_fact;
	  //final[i].p_bin[k].wt_num.err=b1[i].p_bin[k].wt_num.err+b2[i].p_bin[k].wt_num.err*b2_fact;

	  // final[i].p_bin[k].num_den=b1[i].p_bin[k].num_den+b2[i].p_bin[k].num_den*b2_fact;

	  final[i].p_bin[k].ls_err=b1[i].p_bin[k].ls_err+b2[i].p_bin[k].ls_err*b2_fact;

          for(int j=0;j<7;j++)
            {
	      final[i].p_bin[k].data[j].val=b1[i].p_bin[k].data[j].val+b2[i].p_bin[k].data[j].val*b2_fact;
	      final[i].p_bin[k].data[j].err=b1[i].p_bin[k].data[j].err+b2[i].p_bin[k].data[j].err*b2_fact;
	    }
	}

      for(int j=0;j<7;j++)
	{
	  final[i].data[j].val=b1[i].data[j].val+b2[i].data[j].val*b2_fact;
	  final[i].data[j].err=b1[i].data[j].err+b2[i].data[j].err*b2_fact;
	}
    }
  return 0;
}


int outp_2Dbins(bin b[],data_info data_inf, string FileName)      //output bins to file
{
  ofstream writ2;
  writ2.open(FileName.c_str());
  //  writ2<<std::setprecision(std::numeric_limits<long double>::digits10);
  if (!writ2.is_open())
    {
      cout<<FileName<<" not open"<<endl;
      return 1;
    }

  //column headings
  writ2<<"# rp     pi   npairs     wgg      wgg_err  wgp     wgp_err  wgx  wgx_err     wpp  wpp_err    wxx   wxx_err  wpx   wpx_err    theta    theta_err  wt_npairs "<<endl;

  //output data
  for (int i=0;i<data_inf.n_bins;i++)
    {
      for (int k=0;k<data_inf.n_p_bin;k++)
	     {
	        writ2<<(b[i].b_max+b[i].b_min)/2<<"        "<<(b[i].p_bin[k].p_min+b[i].p_bin[k].p_max)/2.0<<"         "<<b[i].p_bin[k].num;

	        for(int j=0;j<7;j++)
	         {
	            writ2<<"       "<<b[i].p_bin[k].data[j].val<<"      "<<b[i].p_bin[k].data[j].err;
	         }
	        writ2<<"   "<<b[i].p_bin[k].wt_num<<endl;
	     }
    }
  writ2.close();
  return 0;
}
