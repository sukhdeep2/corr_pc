//some calculations: like angular distace b/w two galaxies

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "data_def.h"

extern const double r2d=180.0/M_PI;

float angle_rad(angl &a) // angel in radians from struct angl
{
  return 0;   //struct angle doesn't support this now
  //  return ((float)a.deg+((float)a.min/60)+((float)a.sec/3600))*(3.142/180);
}

double angle_deg_to_rad(angl &a) // from struct angl.
{
  return (a.val_deg/r2d);
}


float angle_deg_to_rad(float &a) // from angle given in degrees (not struct angl)
{
  return (a/r2d);
}

double angle_deg_to_rad(double &a) // from angle given in degrees (not str t angl)
{
  return (a/r2d);
}

double arctan(double &x,double &y)//return atan in range [-pi,pi] instead of pi/2
{
  double t=atan(y/x);
  if (x<0)
    {
      if (y<0)
	        t-=M_PI;
      else
	        t+=M_PI;
   }
  return t;
}

double absolute(double n)	//to return |n|
{
  if (n<0) return -1*n;
  else  return n;
}

double angle_orient_axis_PB(gal &g1,gal &g2)   //angle of line joining two galaxies with x-axis
{
  return atan((g2.dec.val_deg-g1.dec.val_deg)/(g2.ra.val_deg-g1.ra.val_deg));
}


void elip_rotate(double e_initial[], double e_final[], double &theta) //rotate ellipse in e_initial (2 ellipticites) by theta and store in e_final
{
  e_final[0]=cos(2.0*theta)*e_initial[0]-sin(2.0*theta)*e_initial[1]; //ellipticites along new axis
  e_final[1]=sin(2.0*theta)*e_initial[0]+cos(2.0*theta)*e_initial[1];
}


void elip_rotate(double e_initial, double e_final[], double &theta) //rotate ellipse in e_initial (e0) by theta and store in e_final
{
  e_final[0]=cos(2.0*theta)*e_initial; //ellipticites along new axis
  e_final[1]=-1.0*sin(2.0*theta)*e_initial;
}

double elip_angle(gal &g1)
{
  double e1=-1*g1.e[1];
  g1.phi=0.5*arctan(g1.e[0],e1);//arctan defined above.. returns in range [-pi,pi]
  if (g1.phi!=g1.phi)
    {
      //cout<<"calcs: gal phi infinity"<<endl;
      g1.phi=0;
    }
  return g1.phi;
}

void sig_crit_inv(calc_temp& ct,double &zl, data_info &data_inf) //sig_crit inverse
{
  if (data_inf.use_comoving)
    ct.sig_crit_inv=ct.dl*ct.dls*(1+zl)/ct.ds;
                                        /* 1+zl factor enters here only once when using comoving
                                            as D_A=D_C/(1+z). (1+z) cancels out one of the 1+zl
                                            factors. ds and dls cancel out 1+zs factors. The
                                            additional 1+zl come from constants.*/
  else
    ct.sig_crit_inv=ct.dl*ct.dls/ct.ds;//physical coordinates
}

double ephi_calc(gal &g1,calc_temp&ct)//return smallest angle between elipticity and line, in range [0,pi/2]
{
  ct.temp[0]=absolute(g1.phi-ct.theta);//both phi and theta are in range [-pi/2,pi/2]
  if (ct.temp[0]>M_PI/2.0)
    {
      ct.temp[0]=M_PI-ct.temp[0];
    }
  return ct.temp[0]; //cos(ct.temp[0]);
}


void pi_calc_projected(gal &g1,gal &g2,calc_temp &ct, data_info &data_inf)//same as pi_calc except for z_avg
{
  ct.delta_z=(g2.redshift-g1.redshift);
  ct.z_avg=g1.redshift;
  //ct.d_P=g2.DC-g1.DC;
  ct.d_P=0;//comoving_p(ct.delta_z,ct.z_avg,data_inf);  // comoving line of sight distance
}

void pi_calc_phi(gal &g1,gal &g2,calc_temp &ct, data_info &data_inf)
{
  pi_calc_projected(g1,g2,ct,data_inf);
  ct.theta=ct.angle_orient_axis(g1,g2);
  if(data_inf.g_with_shape==2)
    ct.d_P=ephi_calc(g2,ct);
  else
    ct.d_P=ephi_calc(g1,ct);
}

void pi_calc(gal &g1,gal &g2,calc_temp &ct, data_info &data_inf)
{
  ct.delta_z=(g2.redshift-g1.redshift);
  ct.z_avg=(g1.redshift+g2.redshift)/2.0;
  //  ct.d_P=comoving_dist(g1.redshift,g2.redshift,data_inf);
  ct.d_P=g2.DC-g1.DC;//comoving distance
}


void r_mu_calc(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct)
{
  ct.rp_calc(g1,g2,data_inf,ct);
  ct.da_R=sqrt(ct.da_R*ct.da_R+ct.d_P*ct.d_P);
  ct.d_P/=ct.da_R;
}

void sig_crit_calc(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct)
{/* Comoving distances only*/
  ct.temp[0]=0;
  ct.dl=g1.DA;//comoving_dist(ct.temp[0],g1.redshift,data_inf);//lens
  ct.ds=g2.DA;//comoving_dist(ct.temp[0],g2.redshift,data_inf); //source
  ct.dls=ct.ds-ct.dl;//comoving_dist(g1.redshift,g2.redshift,data_inf);

  if (!data_inf.use_comoving)
    ct.dls=1./(1+g2.redshift)*(g2.DC-g1.DC);//valid only in flat LCDM

  sig_crit_inv(ct,g1.redshift,data_inf); //sig_crit_inv

  ct.wt_f.val*=pow(ct.sig_crit_inv,2);//*g1.wt*g2.wt;
  ct.sig_crit=1./ct.sig_crit_inv;
  g1.lens_wt+=ct.wt_f.val;
  // cout<<g1.redshift<<"  "<<g2.redshift<<endl;
  // cout<<ct.sig_crit<<"  "<<ct.dl<<"  "<<ct.ds<<"   "<<ct.dls<<endl;
}

void sig_crit_g1(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct)
{/* Comoving distances only*/
  ct.temp[0]=0;
  ct.dl=g1.DA;//comoving_dist(ct.temp[0],g1.redshift,data_inf);//lens
  ct.ds=data_inf.CMB_DIST;
  ct.dls=ct.ds-ct.dl;//comoving_dist(g1.redshift,g2.redshift,data_inf);
  sig_crit_inv(ct,g1.redshift,data_inf); //sig_crit_inv
  ct.wt_f.val*=pow(ct.sig_crit_inv,2);//*g1.wt*g2.wt;
  ct.sig_crit=1./ct.sig_crit_inv;
  g1.lens_wt+=ct.wt_f.val;
  if (isinf(ct.sig_crit))
    {
      ct.sig_crit=0;//because 0*inf=nan
      ct.fail=1;//useless pair
    }
}
