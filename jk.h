
void jk_initilize(data_info &data_inf);

int corel_jk(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct);

double jk_prob_PB();
void seed();

void jk_check_PB(calc_temp &ct);
void jk_check_sky(calc_temp &ct);
//void jk_check_sky_wl(calc_temp &ct);
void no_jk(calc_temp &ct);
int jk2_indx(data_info &data_inf,int ij[2]);
void indx_jk2(data_info &data_inf,int ij[2],int k);
void final_jk_bins(bin_jk &bjk, data_info &data_inf); //subtract jk region from full to get jk sample
