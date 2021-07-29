#include <math.h>
#define _USE_MATH_DEFINES
#include "data_def.h"

double angle_orient_axis_sky(gal &g1,gal &g2);	//angle of line joining two galaxies with x-axis (RA)

void lim_calc_ra_dec(gal &g1,calc_temp &ct,data_info &data_inf);

void lim_calc_sky(gal &g1,calc_temp &ct,data_info &data_inf);

void lim_calc_sky_wl(gal &g1,calc_temp &ct,data_info &data_inf);

double dA_calc(gal &g1, gal &g2);	//angular distance in degrees. Wikipedia Great_circle formula

double dAR_calc(double &da, double &z, data_info &data_inf);	//calculate angular diameter distance at redshift z for angular separation da

double comoving_p(double &d_z,double &z,data_info &data_inf); // comoving distance b/w two objects wth separation d_z, at mean redshift z

double comoving_dist(double &z_l,double &z_h,data_info &data_inf);

double dAR_to_dA(double &dar, double &z, data_info &data_inf);

double p_to_dz(double &p, double &z, data_info &data_inf);

void rp_calc_sky(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct);
void rp_calc_sky_projected(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct);
void rp_calc_sky_angular(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct);
