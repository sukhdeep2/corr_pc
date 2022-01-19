#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "data_def.h"
#include "calcs.h"

const double r2d=180.0/M_PI;

void lim_calc_xyz(gal &g1,calc_temp &ct,data_info &data_inf)
{
  ct.PBx=0;  ct.PBy=0;   ct.PBz=0;
  if (data_inf.coordinates<3||data_inf.coordinates==6) //FIXME
    {
      ct.z_max=g1.redshift+data_inf.p_max;
      ct.z_min=g1.redshift+data_inf.p_min;
    }
  else{
    ct.z_max=g1.redshift+data_inf.bin_r_max;
    ct.z_min=g1.redshift-data_inf.bin_r_max;
  }
  ct.ra_max=g1.ra.val_deg+data_inf.bin_r_max;
  ct.ra_min=g1.ra.val_deg-data_inf.bin_r_max;
  ct.dec_max=g1.dec.val_deg+data_inf.bin_r_max;
  ct.dec_min=g1.dec.val_deg-data_inf.bin_r_max;
}

void lim_calc_PB(gal &g1,calc_temp &ct,data_info &data_inf)
{
  ct.gjk[0]=g1.jk;
  ct.jk_prob[0]=g1.jk_prob;
  lim_calc_xyz(g1,ct,data_inf);

  if (ct.z_min<0)
    {
      ct.PBz=-1;
    }
  if (ct.z_max>data_inf.periodic_box_size)
    {
      ct.PBz=abs(ct.PBz)+1;
    }

  if (ct.ra_min<0)
    {
      ct.PBx=-1;
    }
  if (ct.ra_max>data_inf.periodic_box_size)
    {
      ct.PBx=abs(ct.PBx)+1;
    }
  if (ct.dec_min<0)
    {
      ct.PBy=-1;
    }
  if (ct.dec_max>data_inf.periodic_box_size)
    {
      ct.PBy=abs(ct.PBy)+1;
    }
}

void RR_PB_calc(bin rr[], data_info &data_inf) //periodic box RR pair count expectation  //doesnot work for r-mu
{
  double area=0,volume=0;
  double total_volume=pow(data_inf.periodic_box_size,3);
  for(int i=0;i<data_inf.n_bins;i++)
    {
      if (data_inf.coordinates==6)
        area=M_PI*(rr[i].b_max*rr[i].b_max-rr[i].b_min*rr[i].b_min);
      else if(data_inf.coordinates==7)
        area=4./3.*M_PI*(rr[i].b_max*rr[i].b_max*rr[i].b_max-rr[i].b_min*rr[i].b_min*rr[i].b_min);
      rr[i].num=1;
      for(int k=0;k<data_inf.n_p_bin;k++)
        {
          if (data_inf.coordinates==6)
            volume=area*(rr[i].p_bin[k].p_max-rr[i].p_bin[k].p_min);
          else if(data_inf.coordinates==7)
            volume=area; //FIXME
          rr[i].p_bin[k].wt_num=volume/total_volume; //no need to use pair count. SD_wt (or SS_wt) will do the job
	         rr[i].p_bin[k].num=1;
	  //	  cout<<"RR_PB_calc  "<<volume<<"   "<<total_volume<<"    "<<rr[i].p_bin[k].wt_num<<endl;
        }
    }
}

void pi_calc_PB(gal &g1,gal &g2,calc_temp &ct, data_info &data_inf)//periodic boundary
{
  ct.delta_z=((g2.redshift+data_inf.periodic_box_size*ct.PBz)-g1.redshift);
  ct.d_P=ct.delta_z;
}

void rp_calc_PB(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct)
{
  ct.temp[0]=g1.ra.val_deg-(g2.ra.val_deg+data_inf.periodic_box_size*ct.PBx);
  ct.temp[1]=g1.dec.val_deg-(g2.dec.val_deg+data_inf.periodic_box_size*ct.PBy);
  ct.da_R=sqrt( ct.temp[0]*ct.temp[0] + ct.temp[1]*ct.temp[1]);
}


double ED_calc_PB(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct)
{
  double xyz[3];
  xyz[0]=g1.ra.val_deg-(g2.ra.val_deg+data_inf.periodic_box_size*ct.PBx);
  xyz[1]=g1.dec.val_deg-(g2.dec.val_deg+data_inf.periodic_box_size*ct.PBy);
  xyz[2]=g1.redshift- (g2.redshift+data_inf.periodic_box_size*ct.PBz);
  

  //xyz[0]=g1.ra.val_deg-g2.ra.val_deg;
  //xyz[1]=g1.dec.val_deg-g2.dec.val_deg;
  //xyz[2]=g2.redshift-g1.redshift;
  
  double xyz_norm=sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
  double cos_ang=-10;
  if (data_inf.g_with_shape==2){
    cos_ang=g2.e[0]*xyz[0]+g2.e[1]*xyz[1]+g2.e[2]*xyz[2];
  }
  else
    cos_ang=g1.e[0]*xyz[0]+g1.e[1]*xyz[1]+g1.e[2]*xyz[2];
  
  cos_ang/=xyz_norm; //make sure g.e is normalized as well

  /*
  if (abs(cos_ang)>1){
    cout<<"doing ED calc "<<cos_ang<<"  "<<xyz_norm<<"  "<<endl;
    }

    cout<<"doing ED calc "<<cos_ang<<endl;
     cout<<"doing ED calc dis: "<<xyz_norm<<endl;
     cout<<"doing ED calc dis: "<<xyz[0]<<"  "<<xyz[1]<<"  "<<xyz[2]<<endl;
     cout<<"doing ED calc dis: "<<g2.e[0]<<"  "<<g2.e[1]<<"  "<<g2.e[2]<<endl;
*/
  return cos_ang*cos_ang;
}
