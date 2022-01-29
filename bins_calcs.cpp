#include <omp.h>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "data_def.h"
#include "calcs.h"
#include "sky_calcs.h"
#include "PB_calcs.h"
#include "bins_calcs.h"

void final_calc(bin sd[], bin sr[],bin final[], data_info &data_inf)// SD/RD-1
{
  data_inf.SD_wt.val=data_inf.S_wt.val;// we only need normalization by relative weights of lens and randoms
  data_inf.SD_wt.err=data_inf.SD_wt.val*sqrt(data_inf.S_wt.err/pow(data_inf.S_wt.val,2));

  data_inf.DRs_wt.val=data_inf.Rs_wt.val;
  data_inf.DRs_wt.err=data_inf.DRs_wt.val*sqrt(data_inf.Rs_wt.err/pow(data_inf.Rs_wt.val,2));

  data_inf.SRd_wt.val=data_inf.DRs_wt.val;
  data_inf.SRd_wt.err=data_inf.DRs_wt.err;

  data_inf.RsRd_wt.val=data_inf.DRs_wt.val;
  data_inf.RsRd_wt.err=data_inf.DRs_wt.err;

  final_calc_bins(sd,sr,sr,sr,final,data_inf);
  int_calc(sd,sr,sr,sr,final,data_inf);
}

void final_calc(bin sd[], bin sr[],bin dr[],bin rr[], bin final[], data_info &data_inf)// calculations to get final answers.
{
  data_inf.SD_wt.val=data_inf.S_wt.val*data_inf.D_wt.val;
  data_inf.SD_wt.err=data_inf.SD_wt.val*sqrt(data_inf.S_wt.err/pow(data_inf.S_wt.val,2)+data_inf.D_wt.err/pow(data_inf.D_wt.val,2));

  data_inf.DRs_wt.val=data_inf.Rs_wt.val*data_inf.D_wt.val;
  data_inf.DRs_wt.err=data_inf.DRs_wt.val*sqrt(data_inf.Rs_wt.err/pow(data_inf.Rs_wt.val,2)+data_inf.D_wt.err/pow(data_inf.D_wt.val,2));

  data_inf.SRd_wt.val=data_inf.S_wt.val*data_inf.Rd_wt.val;
  data_inf.SRd_wt.err=data_inf.SRd_wt.val*sqrt(data_inf.S_wt.err/pow(data_inf.S_wt.val,2)+data_inf.Rd_wt.err/pow(data_inf.Rd_wt.val,2));

  data_inf.RsRd_wt.val=data_inf.Rs_wt.val*data_inf.Rd_wt.val*data_inf.include_prob_RR*data_inf.include_prob_RR;
  data_inf.RsRd_wt.err=data_inf.RsRd_wt.val*sqrt(data_inf.Rs_wt.err/pow(data_inf.Rs_wt.val,2)+data_inf.Rd_wt.err/pow(data_inf.Rd_wt.val,2));

  if (data_inf.periodic_box)
    {
      data_inf.RsRd_wt.val=1;
      data_inf.RsRd_wt.err=0;
      RR_PB_calc(rr,data_inf);
    }
  final_calc_bins(sd,sr,dr,rr,final,data_inf);
  int_calc(sd,sr,dr,rr,final,data_inf);
}

void final_calc_auto(bin ss[], bin sr[], bin rr[], bin final[], data_info &data_inf)// calculations to get final answers.
{
  double SS_wt=data_inf.S_wt.val*data_inf.S_wt.val/2.0;
  double SS_err=SS_wt*sqrt(data_inf.S_wt.err/pow(data_inf.S_wt.val,2)+data_inf.S_wt.err/pow(data_inf.S_wt.val,2));

  double SR_wt=data_inf.S_wt.val*data_inf.Rs_wt.val;
  double SR_err=SR_wt*sqrt(data_inf.S_wt.err/pow(data_inf.S_wt.val,2)+data_inf.Rs_wt.err/pow(data_inf.Rs_wt.val,2));

  double RR_wt=data_inf.Rs_wt.val*data_inf.Rs_wt.val*data_inf.include_prob_RR*data_inf.include_prob_RR;
  double RR_err=RR_wt*sqrt(data_inf.Rs_wt.err/pow(data_inf.Rs_wt.val,2)+data_inf.Rs_wt.err/pow(data_inf.Rs_wt.val,2));

  if (!data_inf.RsRd_same)
    {
      RR_wt/=2.0;
      RR_err/=2.0;
    }
  if (data_inf.periodic_box)
    {
      RR_wt=1;
      //      SS_wt*=2.0;
      RR_err=0;
      RR_PB_calc(rr,data_inf);
      if (data_inf.do_auto_rebin) //auto_rebin is useless for periodic box. All pairs are counted.
        re_bin(rr,data_inf);
    }
  data_inf.SD_wt.val=SS_wt;
  data_inf.SRd_wt.val=SR_wt;
  data_inf.RsRd_wt.val=RR_wt;
  data_inf.DRs_wt.val=SR_wt;

  data_inf.SD_wt.err=SS_err;
  data_inf.SRd_wt.err=SR_err;
  data_inf.RsRd_wt.err=RR_err;
  data_inf.DRs_wt.err=SR_err;

  final_calc_bins(ss,sr,sr,rr,final,data_inf); //dr used only in 'gg' calculations, for which using DR=SR same as 2*SR
  int_calc(ss,rr,sr,final,data_inf);
}


void final_calc_bins(bin sd[], bin sr[],bin dr[],bin rr[], bin final[], data_info &data_inf)// calculations to get final answers.
{
  double SD_wt=data_inf.SD_wt.val;
  double SR_wt=data_inf.SRd_wt.val;
  double DR_wt=data_inf.DRs_wt.val;
  double RR_wt=data_inf.RsRd_wt.val;

  double SD_err=data_inf.SD_wt.err;
  double SR_err=data_inf.SRd_wt.err;
  double DR_err=data_inf.DRs_wt.err;
  double RR_err=data_inf.RsRd_wt.err;

  for(int i=0;i<data_inf.n_bins;i++)
    {
      for(int k=0;k<data_inf.n_p_bin;k++)
	{
	  final[i].num+=sd[i].p_bin[k].num;//sum of all p bins. To be compared with value calculated separately as cross check
	  //final[i].b_err+=(1.0/(double)sd[i].p_bin[k].num);
    //      final[i].b_err+=(1.0/(double)sr[i].p_bin[k].num);

	  double sd_num=(double)sd[i].p_bin[k].num;
 	  double sr_num=(double)sr[i].p_bin[k].num;
          double dr_num=(double)dr[i].p_bin[k].num;
          double rr_num=(double)rr[i].p_bin[k].num;

          double sd_wt_num=sd[i].p_bin[k].wt_num;//weighted pair count
          double sr_wt_num=sr[i].p_bin[k].wt_num;
          double dr_wt_num=dr[i].p_bin[k].wt_num;
          double rr_wt_num=rr[i].p_bin[k].wt_num;

          // double sd_wt_err=sd[i].p_bin[k].wt_num.err;
          // double sr_wt_err=sr[i].p_bin[k].wt_num.err;
          // double dr_wt_err=dr[i].p_bin[k].wt_num.err;
          // double rr_wt_err=rr[i].p_bin[k].wt_num.err;

	  if (rr_num!=0)//avoid division by zero
	    {	//wgg
	      double sd_num_rel=(sd_wt_num/rr_wt_num)*(RR_wt/SD_wt); //SD/RR
              double sr_num_rel=(sr_wt_num/rr_wt_num)*(RR_wt/SR_wt); //SR/RR
              double dr_num_rel=(dr_wt_num/rr_wt_num)*(RR_wt/DR_wt); //DR/RR

	      // double sd_rel_err=sd_wt_err*pow(1.0/rr_wt_num*(RR_wt/SD_wt),2.0)+pow(sd_num_rel,2.0)*(rr_wt_err/pow(rr_wt_num,2)+SD_err/pow(SD_wt,2)+RR_err/pow(RR_wt,2));
	      // double sr_rel_err=sr_wt_err*pow(1./rr_wt_num*(RR_wt/SR_wt),2)+pow(sr_num_rel,2.0)*(rr_wt_err/pow(rr_wt_num,2)+SR_err/pow(SR_wt,2)+RR_err/pow(RR_wt,2));
	      // double dr_rel_err=dr_wt_err*pow(1/rr_wt_num*(RR_wt/DR_wt),2)+pow(dr_num_rel,2.0)*(rr_wt_err/pow(rr_wt_num,2)+DR_err/pow(DR_wt,2)+RR_err/pow(RR_wt,2));
	      // double rr_rel_err=0;

	      double w_dp=sd_num_rel-sr_num_rel-dr_num_rel+1.0; //wgg

	      // double num_var=sd_rel_err+sr_rel_err+dr_rel_err+rr_rel_err;//wgg variance

	      if(data_inf.periodic_box)
		{
		  w_dp=sd_num_rel-1.0;
		  //num_var=sd_wt_err*pow(1.0/rr_wt_num*(RR_wt/SD_wt),2.0);//sd_wt_err/rr_wt_num*(RR_wt/SD_wt);
		}

	      final[i].p_bin[k].data[0].val=w_dp;

	      // double w_err=sqrt(num_var);
	      // final[i].p_bin[k].data[0].err=w_err;// W_delta (r,p)
	      //	      final[i].p_bin[k].ls_err=sqrt(pow(1+w_dp,3)/sr_wt_num);// LS estimator delta (r,p)

	      sd[i].p_bin[k].data[0].val=sd_num_rel;
	      sr[i].p_bin[k].data[0].val=sr_num_rel;
	      dr[i].p_bin[k].data[0].val=dr_num_rel;
	      rr[i].p_bin[k].data[0].val=1.0;

	      // sd[i].p_bin[k].data[0].err=sqrt(sd_rel_err);
	      // sr[i].p_bin[k].data[0].err=sqrt(sr_rel_err);
	      // dr[i].p_bin[k].data[0].err=sqrt(dr_rel_err);
	      // rr[i].p_bin[k].data[0].err=sqrt(rr_rel_err);

	      for(int j=1;j<7;j++)//'e' calculations
		{
		  if(data_inf.which_corr==0)
		    {
		      continue;
		    }

		  double sd_val=sd[i].p_bin[k].data[j].val;

		  double sr_val=sr[i].p_bin[k].data[j].val;

		  double sd_sum_sq=sd[i].p_bin[k].data[j].err;
		  double sr_sum_sq=sr[i].p_bin[k].data[j].err;

		  double sd_val_sq=pow(sd_val,2);
		  double sr_val_sq=pow(sr_val,2);

		  double sd_err2= sd_sum_sq - (sd_val_sq/sd_num); //variance
		  double sr_err2= sr_sum_sq - (sr_val_sq/sr_num);
		  double sd_val_norm=0;
		  double sr_val_norm=0;
		  if (data_inf.which_corr==5 && j==6)
                    {
		      if(sd_num>0){
			sd_val_norm=(sd_val/sd_wt_num);
		      }
		      sd_val_norm-=1./3;
                    }

		  else{
		    sd_val_norm=(sd_val/rr_wt_num)*(RR_wt/SD_wt); //SD/RR
		    sr_val_norm=(sr_val/rr_wt_num)*(RR_wt/SR_wt); //SR/RR
		  }

		  double sd_val_norm2=(sd_val); //don't take mean here. Will give higher weight to noise.
		  //to deal with diviion with zero
                  if(sd_num==0)
                    {
                      sd_err2=0;
                    }
                  if(sr_num==0)
		    {
                      sr_err2=0;
                    }

		  // double sd_rel_err=sd_err2*pow((1/rr_wt_num)*(RR_wt/SD_wt),2)+pow(sd_val_norm,2)*(SD_err/pow(SD_wt,2)+RR_err/(pow(RR_wt,2))+rr_wt_err/pow(rr_wt_num,2));
      //             double sr_rel_err=sr_err2*pow(1/rr_wt_num*RR_wt/SR_wt,2.0)+pow(sr_val_norm,2)*(SR_err/pow(SR_wt,2)+RR_err/(pow(RR_wt,2))+rr_wt_err/pow(rr_wt_num,2));
      //
		  // double num_err_t=sd_rel_err+sr_rel_err;//variance
		  double val_f=sd_val_norm-sr_val_norm;//w_aa..
		  // double err_f=sqrt(num_err_t); //final error

		  if(sd_num==0||sr_num==0||rr_num==0)//bins with no data
		    {
		      val_f=0;
		      // err_f=0;
		      if (data_inf.periodic_box)
			{
			  val_f=sd_val_norm;
			  // err_f=sqrt(sd_rel_err);
			}
		    }

		  final[i].p_bin[k].data[j].val=val_f;
		  // final[i].p_bin[k].data[j].err=err_f;
		  sd[i].p_bin[k].data[j].val=sd_val_norm2;
		  // sd[i].p_bin[k].data[j].err=sqrt(sd_rel_err);
		  sr[i].p_bin[k].data[j].val=sr_val_norm;
		  // sr[i].p_bin[k].data[j].err=sqrt(sr_rel_err);

		  if(j==6)
		    {
		      sd[i].p_bin[k].data[j].val=sd_val;// don't take mean here. That will give higher weight to noise.
                      sr[i].p_bin[k].data[j].val=sr_val;//

		      if (sd_num==0)
			{
			  sd[i].p_bin[k].data[j].val=0;
			}

                      if (sr_num==0)
                        {
                          sr[i].p_bin[k].data[j].val=0;
                        }
		    }
		}
	    }
	}
      if(final[i].num!=sd[i].num)cout<<"number counting problem at bin: "<<i<<endl;
    }
}

void final_calc_mean_projected(bin sd[],bin sr[],bin final[],data_info data_inf)// get binned data in order
{
  int_calc(sd,data_inf);
  int_calc(sr,data_inf);

  double sd_sum_sq=0,sd_sq=0,sd_var=0,sd_num=0,final_var=0,wt_num_var=0,wt_num_sq=0;
  double sr_sum_sq=0,sr_sq=0,sr_var=0,sr_num=0,boost_fact=1;
  bool do_boost=0;
  //  double nS=(double) data_inf.n_shape;
  //double nR=(double) data_inf.n_shape_rand;
  double nS=(double) data_inf.n_density;
  double nR=(double) data_inf.n_density_rand;

  if (data_inf.which_corr>=7)
    {
      do_boost=1;
    }
  for(int i=0;i<data_inf.n_bins;i++)
    {
      if(sd[i].num>0 || sr[i].num>0)
	     {
      	  sd_sum_sq=sd[i].data[0].err;
      	  sd_sq=pow(sd[i].data[0].val,2.0);
      	  sd_num=(double)sd[i].num;
      	  sd_var=sd_sum_sq-(sd_sq/sd_num);//var of numerator

          sr_sum_sq=sr[i].data[0].err;
          sr_sq=pow(sr[i].data[0].val,2.0);
          sr_num=(double)sr[i].num;
          sr_var=sr_sum_sq-(sr_sq/sr_num);//var of numerator

      	  wt_num_sq=pow(sd[i].wt_num,2.0);
      	  //wt_num_var=sd[i].wt_num.err-wt_num_sq/sd_num; //variance of denominator

      	  // final_var=sd_var/sd_sq+wt_num_var/wt_num_sq;//total relative variance

          if (sd_num>0)
            final[i].sig_crit=sd[i].sig_crit/sd[i].wt_num;

      	  for(int j=0;j<data_inf.n_bin_var;j++)
      	   {
      	      final[i].data[j].val=0;
      	      if (sd_num>0)
      		      final[i].data[j].val=sd[i].data[j].val/sd[i].wt_num;// normalized value
      	      final[i].data[j].err=final[i].data[j].val*sqrt(final_var); //final error (poisson)... This is wrong
      	      if (sr_num>0)
      		      {
      		        if (do_boost)
      		          boost_fact=sd[i].wt_num/sr[i].wt_num*nR/nS;
      		        final[i].data[j].val*=boost_fact;//default boost_fact 1
      		        final[i].data[j].val-=sr[i].data[j].val/sr[i].wt_num;
      		      }
      	    }
      	}
      else //in case of empty bin
	     {
    	  for(int j=0;j<data_inf.n_bin_var;j++)
        {
    	      final[i].data[j].val=0;
    	      final[i].data[j].err=0;
    	    }
	     }

      final[i].num=sd[i].num; //number of pairs
      final[i].wt_num=sd[i].wt_num; //weighted number of pairs

      if(final[i].data[0].val!= final[i].data[0].val)//check of nan values
	     {
	        cout<<"NaN Problem in bin "<<i<<"  "<<endl;
	      }
    }
}

void int_calc(bin dat[],data_info data_inf)// integrate over p bins
{
  double dC=data_inf.p_size;
  if (data_inf.coordinates!=0&&data_inf.coordinates!=6)
    dC=1;
  for(int i=0;i<data_inf.n_bins;i++)
    {
      for(int k=0;k<data_inf.n_p_bin;k++)//integrating over all p_bins
	{
	  // if (data_inf.periodic_box && data_inf.do_RR)
	  dat[i].wt_num+=dat[i].p_bin[k].wt_num;
	  for(int j=0;j<7;j++)
	    {
	      dat[i].data[j].val+=dat[i].p_bin[k].data[j].val*dC;
	      dat[i].data[j].err+=pow(dat[i].p_bin[k].data[j].err*dC,2);
	    }
	  dat[i].ls_err_int+=pow(dat[i].p_bin[k].ls_err*dC,2);
	}
      for(int j=0;j<7;j++) 	//sqrt to get final errors
	{
	  //	  if (data_inf.do_density_only)continue;
	  dat[i].data[j].err=sqrt(dat[i].data[j].err);
	}
      dat[i].ls_err_int=sqrt(dat[i].ls_err_int);
      //dat[i].b_err=sqrt(dat[i].b_err)*dC;
    }
}


void int_calc(bin sd[], bin sr[],bin dr[],bin rr[], bin final[], data_info &data_inf)// integrate over p bins
{
  data_inf.do_SD=1;
  int_calc(sd,data_inf);
  data_inf.do_SD=0;
  int_calc(sr,data_inf);
  int_calc(dr,data_inf);
  //data_inf.do_RR=1;
  int_calc(rr,data_inf);
  //data_inf.do_RR=0;
  int_calc(final,data_inf);
}

void int_calc(bin ss[], bin sr[],bin rr[], bin final[], data_info &data_inf)// integrate over p bins
{
  data_inf.do_SS=1;
  int_calc(ss,data_inf);
  data_inf.do_SS=0;
  int_calc(sr,data_inf);
  //data_inf.do_RR=1;
  int_calc(rr,data_inf);
  //data_inf.do_RR=0;
  int_calc(final,data_inf);
}

void int_calc(bin ss[], bin sr[], data_info &data_inf)// integrate over p bins
{
  data_inf.do_SS=1;
  int_calc(ss,data_inf);
  data_inf.do_SS=0;
  int_calc(sr,data_inf);
}

void jk_final(bin bins_final_jk[] ,bin bins_final[],jackknife jk_final[],data_info data_inf)
{
  for(int i=0;i<data_inf.n_bins;i++)
    {
      for(int j=0;j<data_inf.n_bin_var;j++)
	{
	  double bias=0;
    double var=0;
	  double t=0;
	  double mean=0;

	  for(int k=0;k<data_inf.n_jk;k++)
	    {
	      t=(jk_final[k].b[i].data[j].val-bins_final[i].data[j].val);
	      bias+=t;
	      mean+=jk_final[k].b[i].data[j].val;
	    }

	  mean/=data_inf.n_jk;

	  for(int k=0;k<data_inf.n_jk;k++)
            {
              t=(jk_final[k].b[i].data[j].val-mean);
              var+=pow(t,2.0);
            }

	  bias*=((data_inf.n_jk-1.0)/(data_inf.n_jk));
	  var*=((data_inf.n_jk-1.0)/(data_inf.n_jk));
	  bins_final_jk[i].data[j].val=bins_final[i].data[j].val; //-bias;
	  bins_final_jk[i].data[j].err=sqrt(var);
	}
    }
}

void jk_final2D(bin bins_final_jk[] ,bin bins_final[],jackknife jk_final[],data_info data_inf)
{
  for(int i=0;i<data_inf.n_bins;i++)
    {
      for (int l=0;l<data_inf.n_p_bin;l++)
	     {
    	  for(int j=0;j<7;j++)
    	    {
    	      double bias=0;
    	      double var=0;
    	      double t=0;
    	      double mean=0;

    	      for(int k=0;k<data_inf.n_jk;k++)
    		      {
    		          t=(jk_final[k].b[i].p_bin[l].data[j].val-bins_final[i].p_bin[l].data[j].val);
    		          bias+=t;
    		          mean+=jk_final[k].b[i].p_bin[l].data[j].val;
    		      }

	          mean/=data_inf.n_jk;

	          for(int k=0;k<data_inf.n_jk;k++)
		          {
		              t=(jk_final[k].b[i].p_bin[l].data[j].val-mean);
		              var+=pow(t,2.0);
		          }

	          bias*=((data_inf.n_jk-1.0)/(data_inf.n_jk));
	          var*=((data_inf.n_jk-1.0)/(data_inf.n_jk));
	          bins_final_jk[i].p_bin[l].data[j].val=bins_final[i].p_bin[l].data[j].val; //-bias;
	          bins_final_jk[i].p_bin[l].data[j].err=sqrt(var);
	        }
	      }
    }
}


void re_bin(bin ss[], data_info &data_inf)
{
  if(data_inf.n_p_bin<=1)return ;
  for(int i=0;i<data_inf.n_bins;i++)
    {
      for(int k=0;k<data_inf.n_p_bin;k++)
	{

	  if (ss[i].p_bin[k].p_min>=0)
	    {
	      break;//only shift -ve side bins to +ve side
	    }

	  double p_val=0.0-(ss[i].p_bin[k].p_min+ss[i].p_bin[k].p_max)/2.0; //bin center
	  int p=floor((p_val-data_inf.p_min)/data_inf.p_size); //+ve side bin
	binning:

	  if(p<0||p>=data_inf.n_p_bin)
	    {
	      break;
	    }

	  else if (p_val<ss[i].p_bin[p].p_min)
	    {
	      p--;
	      goto binning;
	    }

	  else if(p_val>ss[i].p_bin[p].p_max)
	    {
	      p++;
	      goto binning;
	    }
	  if ( fabs(ss[i].p_bin[k].p_min+ss[i].p_bin[p].p_max)>(data_inf.p_size/100.)  || fabs(ss[i].p_bin[k].p_max+ss[i].p_bin[p].p_min)>(data_inf.p_size/100.)||ss[i].p_bin[p].p_max<0)
	    {
	      cout<<"rebinning problem:"<<i<<"   "<<p<<"   "<<k<<"   "<<p_val<<"   "<<ss[i].p_bin[p].p_min<<endl;
	      //return 1;
	      continue;
	    }

	  ss[i].p_bin[p].num+=ss[i].p_bin[k].num;
	  ss[i].p_bin[k].num=0;
	  ss[i].p_bin[p].wt_num+=ss[i].p_bin[k].wt_num;
          ss[i].p_bin[k].wt_num=0;

          // ss[i].p_bin[p].wt_num.err+=ss[i].p_bin[k].wt_num.err;
          // ss[i].p_bin[k].wt_num.err=0;

	  // ss[i].p_bin[p].num_den+=ss[i].p_bin[k].num_den;
	  ss[i].p_bin[k].num=0;
	  ss[i].p_bin[p].ls_err+=ss[i].p_bin[k].ls_err;
          ss[i].p_bin[k].ls_err=0;

	  for(int j=0;j<7;j++)
	    {
	      ss[i].p_bin[p].data[j].val=ss[i].p_bin[p].data[j].val+ss[i].p_bin[k].data[j].val;
	      ss[i].p_bin[k].data[j].val=0;
              ss[i].p_bin[p].data[j].err=ss[i].p_bin[p].data[j].err+ss[i].p_bin[k].data[j].err;
              ss[i].p_bin[k].data[j].err=0;
	    }
	}
    }
}
