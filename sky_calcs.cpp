#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "data_def.h"
#include "calcs.h"

const double r2d=180.0/M_PI;


double angle_orient_axis_sky(gal &g1,gal &g2)	//angle of line joining two galaxies with x-axis (RA)
{
  /*  double nume=(g2.dec.val_deg-g1.dec.val_deg);
      double denom=((g2.ra.val_deg-g1.ra.val_deg)*cos((g1.dec.val_rad + g2.dec.val_rad) /2.0 ));*/

  return atan((g2.dec.val_deg-g1.dec.val_deg)/((g2.ra.val_deg-g1.ra.val_deg)
					       *cos((g1.dec.val_rad + g2.dec.val_rad) /2.0 )));

    //  return atan(nume/denom);
}

double dA_calc(gal &g1, gal &g2)	//angular distance in degrees. Wikipedia Great_circle formula
{
  double sum=sin(g1.dec.val_rad) * sin(g2.dec.val_rad) +  cos(g1.dec.val_rad)*cos(g2.dec.val_rad) * cos((g1.ra.val_rad-g2.ra.val_rad));

  if(sum>0.9999999999)return acos(1.0); //floating point error

  return acos(sum);
}

double dAR_calc(double &da, double &z, data_info &data_inf)
//calculate angular diameter distance at redshift z for angular separation da
{
  /*  if (z>data_inf.z_max-0.001)
    {
      cout<<z<<" greater than zmax"<<endl;
      exit(1);
      }*/
  int i_t=floor(z/data_inf.dz);

  return (da)*(data_inf.dA_z[i_t]+((data_inf.dA_z[i_t+1]-data_inf.dA_z[i_t])*(z-data_inf.dz*i_t)/data_inf.dz));

  // return da*data_inf.dA_z[i_t];
}


double comoving_dist(double &z_l,double &z_h,data_info &data_inf)            // comoving distance b/w two objects
{
  int i_t=floor(z_l/data_inf.dz);//i_t is index number for array of H_z, hubble scale at redshift z (H as func of z, with z implied from index number).

  double dl=(data_inf.dC_z[i_t]+((data_inf.dC_z[i_t+1]-data_inf.dC_z[i_t])*(z_l-data_inf.dz*i_t)/data_inf.dz));

  i_t=floor(z_h/data_inf.dz);
  double dh=(data_inf.dC_z[i_t]+((data_inf.dC_z[i_t+1]-data_inf.dC_z[i_t])*(z_h-data_inf.dz*i_t)/data_inf.dz));

  return dh-dl;
}

double comoving_p(double &d_z,double &z,data_info &data_inf)		// comoving distance b/w two objects wth separation d_z, at mean redshift z
{
  int i_t=floor(z/data_inf.dz);//i_t is index number for array of H_z, hubble scale at redshift z (H as func of z, with z implied from index number).

  return 300000.0*d_z/(data_inf.H_z[i_t]+((data_inf.H_z[i_t+1]-data_inf.H_z[i_t])*(z-data_inf.dz*i_t)/data_inf.dz));

  //  return 300000.0*d_z/data_inf.H_z[i_t];
}


double dAR_to_dA(double &dar, double &z, data_info &data_inf)//angular diameter distance to angular distance
{
  int i_t=floor(z/data_inf.dz);
  return dar/(data_inf.dA_z[i_t]+((data_inf.dA_z[i_t+1]-data_inf.dA_z[i_t])*(z-data_inf.dz*i_t)/data_inf.dz));
}

double p_to_dz(double &p, double &z, data_info &data_inf)//comoving distance to dz
{
  int i_t=floor(z/data_inf.dz);
  return p*(data_inf.H_z[i_t]+((data_inf.H_z[i_t+1]-data_inf.H_z[i_t])*(z-data_inf.dz*i_t)/data_inf.dz))/300000.0;
}


void rp_calc_sky_projected(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct)
{
  ct.da=dA_calc(g1,g2);//angular distance
  //  ct.da_R=dAR_calc(ct.da,ct.z_avg,data_inf);//angular diameter distance at mean redshift of galaxies.
  ct.da_R=ct.da*g1.DA;//accurate for comoving transverse distance in flat universe
}

void rp_calc_sky(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct)
{
  ct.da=dA_calc(g1,g2);//angular distance
  //  ct.da_R=dAR_calc(ct.da,ct.z_avg,data_inf);//angular diameter distance at mean redshift of galaxies.
  ct.da_R=ct.da*(g1.DA+g2.DA)/2.;//accurate for comoving transverse distance in flat universe
}

void rp_calc_sky_angular(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct)
{
  ct.da=dA_calc(g1,g2);//angular distance
  ct.da_R=ct.da*r2d;
}


void lim_calc_ra_dec(gal &g1,calc_temp &ct,data_info &data_inf)
{
  if (ct.z_min<0.001) //to avoid zeros...
    ct.z_min=0.001;
  if (ct.z_min>data_inf.z_max)
    {
      cout<<"ct z-min > data_inf z-max "<<endl;
      exit(1);
    }
  switch (data_inf.coordinates)
    {
    case 4:
    case 5:
      ct.da_max=data_inf.bin_r_max/r2d;
      break;
    default:
      ct.da_max=dAR_to_dA(data_inf.bin_r_max,ct.z_min,data_inf);
//da_max is greater at lower redshift. using z_min
      break;
    }

  ct.ra_lim=acos((pow(cos(ct.da_max),1)-pow(sin(g1.dec.val_rad),2))/pow(cos(g1.dec.val_rad),2))*1.02*r2d;
// This come from great distance formula after taking derivative wrt g2.dec...
// It is wrong to assume that for ra_max separation g1.dec=g2.dec

  ct.ra_max=g1.ra.val_deg+ct.ra_lim;
  ct.ra_min=g1.ra.val_deg-ct.ra_lim;


  if (ct.ra_max>360){
    ct.ra_min=0;
    ct.ra_lim=360;
  }
  if (ct.ra_min<0){
    ct.ra_max=360;
    ct.ra_lim=360;
  }

  ct.dec_max=g1.dec.val_deg+ct.da_max*r2d*1.05;
  ct.dec_min=g1.dec.val_deg-ct.da_max*r2d*1.05;
}




void lim_calc_sky(gal &g1,calc_temp &ct,data_info &data_inf)
{
  ct.gjk[0]=g1.jk;
  ct.jk_prob[0]=g1.jk_prob;
  ct.PBx=0;ct.PBz=0;ct.PBy=0;
  ct.z_max=g1.redshift+data_inf.z_sep_max;
  switch (data_inf.coordinates)
    {
    case 2:
    case 3:
    case 4:
    case 5:
      ct.z_min=g1.redshift; //calculate angular separation with respect to g1
      lim_calc_ra_dec(g1,ct,data_inf);

      if (data_inf.z_sep_min==0)
	     ct.z_min=data_inf.z_min;
      else
	     ct.z_min=g1.redshift+data_inf.z_sep_min;

      if (data_inf.z_sep_max==0)
	     ct.z_max=data_inf.z_max;
      else
	     ct.z_max=g1.redshift+data_inf.z_sep_max;

      break;

    case 1: //r-mu
      ct.z_max=g1.redshift+p_to_dz(data_inf.bin_r_max,ct.z_max,data_inf);
      ct.z_min=g1.redshift*2-ct.z_max;
      lim_calc_ra_dec(g1,ct,data_inf);
      break;

    default:
      ct.z_max=g1.redshift+p_to_dz(data_inf.p_max,ct.z_max,data_inf);// increases with z. z_max provides some sort of insurance
      ct.z_min=g1.redshift+p_to_dz(data_inf.p_min,g1.redshift,data_inf);
      lim_calc_ra_dec(g1,ct,data_inf);
      break;
    }

  /*  if (ct.z_min<data_inf.z_min)
    {
      ct.z_min=data_inf.z_min;
      }*/

  //lim_calc_ra_dec(g1,ct,data_inf); //this needs to be in switch due to case 2,3,4
}

void lim_calc_sky_wl(gal &g1,calc_temp &ct,data_info &data_inf)
{
  ct.gjk[0]=g1.jk;
  ct.jk_prob[0]=g1.jk_prob;
  ct.PBx=0;ct.PBz=0;ct.PBy=0;

  ct.z_min=g1.redshift;
  lim_calc_ra_dec(g1,ct,data_inf);
  ct.z_min+=data_inf.dz;//to avoid nan problem in sigma crit calcs
  ct.z_min+=data_inf.z_sep_min;
  ct.z_max=g1.redshift+data_inf.z_sep_max;
}
