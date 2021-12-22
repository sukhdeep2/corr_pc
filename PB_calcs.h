#include "data_def.h"

void lim_calc_xyz(gal &g1,calc_temp &ct,data_info &data_inf);
void lim_calc_PB(gal &g1,calc_temp &ct,data_info &data_inf);
void RR_PB_calc(bin rr[], data_info &data_inf);
void pi_calc_PB(gal &g1,gal &g2,calc_temp &ct, data_info &data_inf);//periodic boundary
void rp_calc_PB(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct);
double ED_calc_PB(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct);