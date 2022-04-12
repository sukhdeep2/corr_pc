#include <omp.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define _USE_MATH_DEFINES
#include <time.h>
#include "data_def.h"
#include "calcs.h"
#include "sky_calcs.h"
#include "PB_calcs.h"
#include "initialization.h"
#include "outp.h"
#include <algorithm>    // std::min

using namespace std;

const double pi=M_PI;
const double r2d=180.0/pi;

int corel_bin(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct) // get distances and find bins
{
  ct.pi_calc(g1,g2,ct,data_inf);
  ct.r_calc(g1,g2,data_inf,ct); //ideally should be placed below for speed. Putting here because d_P changes for some coordinates

  ct.n2=floor((ct.d_P-data_inf.p_min)/data_inf.p_size);//bin number
  if(ct.n2<0||ct.n2>=data_inf.n_p_bin)
    return 1;

  //  ct.rp_calc(g1,g2,data_inf,ct);
  if(data_inf.lin_bin)
    ct.n1=floor((ct.da_R-data_inf.bin_r_min)/data_inf.r_size);// get the bin number
  else
    ct.n1=floor(log10(ct.da_R/data_inf.bin_r_min)/data_inf.r_size); //Same, Log bins
  //cout<<"corel_bin:"<<ct.n1<<"  "<<ct.da_R<<"  "<<endl;

  if(ct.n1<0||ct.n1>=data_inf.n_bins)
    return 1;

  ct.gjk[1]=g2.jk;
  ct.jk_prob[1]=g2.jk_prob;
  ct.jk_prob[0]=g1.jk_prob;
  ct.wt_f.val=g1.wt*g2.wt;
  //  ct.sig_crit=1; //set by default.. not necessary here
  return 0;
}

int density_density_corel(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct)
{
  ct.fail=corel_bin(g1,g2,data_inf,ct);

  if (data_inf.sig_crit>0)
    {
      //if(data_inf.g_with_shape==2)
      ct.sig_crit_func(g1,g2,data_inf,ct);//in wl g1 is always lens
      // if(data_inf.g_with_shape==1)
      //ct.sig_crit_func(g2,g1,data_inf,ct);
      if (ct.sig_crit>1.e5)
        cout<<"large sig_crit"<<ct.sig_crit<<endl;
      ct.bjk.b[ct.n1].sig_crit+=ct.wt_f.val*ct.sig_crit;
    }
  if (ct.fail!=0)
    return ct.fail;

  ct.wt_f.err=pow(ct.wt_f.val,2);

  ct.bjk.b[ct.n1].num++;
  //  ct.bjk.b[ct.n1].wt_num+=ct.wt_f.val; //adding pair to dA_R bin
  // ct.bjk.b[ct.n1].wt_num.err+=ct.wt_f.err;

  ct.bjk.b[ct.n1].p_bin[ct.n2].num++;

  ct.bjk.b[ct.n1].p_bin[ct.n2].wt_num+=ct.wt_f.val;//adding pair to dA_R-d_P bin...
  //ct.bjk.b[ct.n1].p_bin[ct.n2].wt_num.err+=ct.wt_f.err;

  return 0;
}

int density_convergence_corel(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct)
{
  ct.fail=density_density_corel(g1,g2,data_inf,ct);
  if (ct.fail!=0)
      return ct.fail;
  if(data_inf.g_with_shape==1)
    ct.ep.val=ct.wt_f.val*g1.e[0]*ct.sig_crit;
  else if(data_inf.g_with_shape==2)
    ct.ep.val=ct.wt_f.val*g2.e[0]*ct.sig_crit;
  ct.bjk.b[ct.n1].p_bin[ct.n2].data[1].val+=ct.ep.val;

  return 0;
}

int shape_convergence_corel(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct) // correlations b/w 2 galaxies
{
  ct.fail=density_density_corel(g1,g2,data_inf,ct);
  if (ct.fail!=0)
    return ct.fail;

  ct.theta=ct.angle_orient_axis(g1,g2);  //angle of line joining centre of two galaxies with x (RA) axis
  if(data_inf.g_with_shape==2)
    {
      elip_rotate(g2.e,ct.e_new_i,ct.theta); //lensing... g2 is the source
      ct.ep.val=ct.wt_f.val*g1.e[0]*ct.sig_crit;//convergence
    }
  else if(data_inf.g_with_shape==1)
    {
      elip_rotate(g1.e,ct.e_new_i,ct.theta);   //rotation when e1 and e2 are known
      ct.ep.val=ct.wt_f.val*g2.e[0]*ct.sig_crit;
    }
  ct.epp.val=ct.ep.val*ct.e_new_i[0];
  ct.epx.val=ct.ep.val*ct.e_new_i[1];

  ct.bjk.b[ct.n1].p_bin[ct.n2].data[1].val+=ct.ep.val;

  ct.bjk.b[ct.n1].p_bin[ct.n2].data[3].val+=ct.epp.val;

  ct.bjk.b[ct.n1].p_bin[ct.n2].data[5].val+=ct.epx.val;

  return 0;
}

int shape_shape_corel(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct)
{
  ct.fail=density_density_corel(g1,g2,data_inf,ct);

  if (ct.fail!=0)
    {
      return ct.fail;
    }

  ct.theta=ct.angle_orient_axis(g1,g2);  //angle of line joining centre of two galaxies with x (RA) axis

  elip_rotate(g1.e,ct.e_new_i,ct.theta);   //rotation when e1 and e2 are known
  //	    elip_rotate(g1.e0,e_new_i,theta-g1.phi);	//rotation when e0 and orientation angle known

  elip_rotate(g2.e,ct.e_new_j,ct.theta);
  //                elip_rotate(g2.e0,e_new_j,theta-g2.phi);

  if (data_inf.estimator==0||data_inf.coordinates==2||data_inf.coordinates==4) //cross correlation, wgp is only with shape sample.
  //  if(data_inf.coordinates==2||data_inf.coordinates==4)
   {//use the galaxy other than one used in phi bin calc
      if(data_inf.g_with_shape==2)
	{
	  ct.ep.val=ct.sig_crit*ct.wt_f.val*(ct.e_new_j[0]); //e+
	  ct.ep.err=pow(ct.wt_f.val*ct.e_new_j[0],2);
	  ct.etheta.val=ephi_calc(g2,ct);
	  ct.etheta.err=pow(ct.etheta.val,2.0);
	}
      else if(data_inf.g_with_shape==1)
	{
	  ct.ep.val=ct.sig_crit*ct.wt_f.val*(ct.e_new_i[0]); //e+
          ct.ep.err=pow(ct.wt_f.val*ct.e_new_i[0],2);
	  ct.etheta.val=ephi_calc(g1,ct);
          ct.etheta.err=pow(ct.etheta.val,2.0);
	}
   }
  else
    {
      ct.ep.val=ct.sig_crit*ct.wt_f.val*(ct.e_new_i[0]+ct.e_new_j[0])/2.0; //e+
      ct.ep.err=pow(ct.wt_f.val*ct.e_new_i[0]/2.0,2)+pow(ct.wt_f.val*ct.e_new_j[0]/2.0,2); //e+

      ct.ex.val=ct.sig_crit*ct.wt_f.val*(ct.e_new_i[1]+ct.e_new_j[1])/2.0; //eX
      ct.ex.err=pow(ct.wt_f.val*ct.e_new_i[1]/2.0,2)+pow(ct.wt_f.val*ct.e_new_j[1]/2.0,2); //eX

      ct.etheta.val=ephi_calc(g1,ct);
      double e2=ephi_calc(g2,ct);

      ct.etheta.err=pow(ct.etheta.val,2.0)+pow(e2,2.0);
      ct.etheta.val+=e2;
    }
      
  ct.epp.val=ct.sig_crit*ct.wt_f.val*(ct.e_new_i[0]*ct.e_new_j[0]); //e++
  ct.epp.err=pow(ct.wt_f.val*(ct.e_new_i[0]*ct.e_new_j[0]),2); //e++
  ct.exx.val=ct.sig_crit*ct.wt_f.val*(ct.e_new_i[1]*ct.e_new_j[1]); // eXX
  ct.exx.err=pow(ct.wt_f.val*(ct.e_new_i[1]*ct.e_new_j[1]),2); // eXX
  ct.epx.val=ct.sig_crit*ct.wt_f.val*(ct.e_new_i[1]*ct.e_new_j[0]+ct.e_new_i[0]*ct.e_new_j[1])/2.0; //e+X
  ct.epx.err=pow(ct.wt_f.val*(ct.e_new_i[1]*ct.e_new_j[0])/2.0,2)+pow(ct.wt_f.val*ct.e_new_i[0]*ct.e_new_j[1]/2.0,2); //e+X

  ct.bjk.b[ct.n1].p_bin[ct.n2].data[1].val+=ct.ep.val; //e+

  ct.bjk.b[ct.n1].p_bin[ct.n2].data[1].err+=ct.ep.err; //e+


  ct.bjk.b[ct.n1].p_bin[ct.n2].data[2].val+=ct.ex.val; //ex

  ct.bjk.b[ct.n1].p_bin[ct.n2].data[2].err+=ct.ex.err; //ex

  ct.bjk.b[ct.n1].p_bin[ct.n2].data[3].val+=ct.epp.val;

  ct.bjk.b[ct.n1].p_bin[ct.n2].data[3].err+=ct.epp.err;

  ct.bjk.b[ct.n1].p_bin[ct.n2].data[4].val+=ct.exx.val;

  ct.bjk.b[ct.n1].p_bin[ct.n2].data[4].err+=ct.exx.err;

  ct.bjk.b[ct.n1].p_bin[ct.n2].data[5].val+=ct.epx.val;

  ct.bjk.b[ct.n1].p_bin[ct.n2].data[5].err+=ct.epx.err;

  ct.bjk.b[ct.n1].p_bin[ct.n2].data[6].val+=ct.etheta.val;
  ct.bjk.b[ct.n1].p_bin[ct.n2].data[6].err+=ct.etheta.err;

  return 0;
}

int shape_density_corel(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct)
{
  ct.fail=density_density_corel(g1,g2,data_inf,ct);

  if (ct.fail!=0)
    {
      return ct.fail;
    }

  ct.theta=ct.angle_orient_axis(g1,g2);  //angle of line joining centre of two galaxies with x (RA) axis
  if (data_inf.g_with_shape==2)
    elip_rotate(g2.e,ct.e_new_i,ct.theta); //lensing... g2 is the source
  else
    elip_rotate(g1.e,ct.e_new_i,ct.theta);   //rotation when e1 and e2 are known

  ct.ep.val=ct.sig_crit*ct.wt_f.val*(ct.e_new_i[0]); //e+

  ct.ep.err=pow(ct.wt_f.val*ct.e_new_i[0],2); //e+

  ct.ex.val=ct.sig_crit*ct.wt_f.val*(ct.e_new_i[1]); //eX

  ct.ex.err=pow(ct.wt_f.val*ct.e_new_i[1],2); //eX

  ct.etheta.val=ephi_calc(g1,ct);
  ct.etheta.err=pow(ct.etheta.val,2.0);

  ct.bjk.b[ct.n1].p_bin[ct.n2].data[1].val+=ct.ep.val; //e+

  ct.bjk.b[ct.n1].p_bin[ct.n2].data[1].err+=ct.ep.err; //e+

  ct.bjk.b[ct.n1].p_bin[ct.n2].data[6].val+=ct.etheta.val;
  ct.bjk.b[ct.n1].p_bin[ct.n2].data[6].err+=ct.etheta.err;

  ct.bjk.b[ct.n1].p_bin[ct.n2].data[2].val+=ct.ex.val; //ex

  ct.bjk.b[ct.n1].p_bin[ct.n2].data[2].err+=ct.ex.err; //ex

  return 0;
}

int ED_corel(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct)
{
  ct.fail=density_density_corel(g1,g2,data_inf,ct);

  if (ct.fail!=0)
    {
      return ct.fail;
    }
  
  ct.etheta.val=ED_calc_PB(g1,g2,data_inf,ct);
  ct.etheta.err=pow(ct.etheta.val,2.0);

  ct.bjk.b[ct.n1].p_bin[ct.n2].data[6].val+=ct.etheta.val;
  ct.bjk.b[ct.n1].p_bin[ct.n2].data[6].err+=ct.etheta.err;

  return 0;
}
