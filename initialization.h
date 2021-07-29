
// intialize data with random values and bins with values 0 ;
int da_Z(data_info &data_inf);// read and store angular diameter distance as array

int data_ini(data_info &data_inf,string FileName); // read input file and intialise main data_info struct

int read_data(gal g[],data_info &data_inf);// read data set and randoms set from file

void P_bin_ini(P_bin p[],data_info &data_inf); // intialize bins

void bin_lin(bin bi[],data_info &data_inf);	//making linear bins

void bin_log(bin bi[],data_info &data_inf); // log bins

void jk_ini(jackknife jk[], data_info &data_inf);

void bin_jk_ini(bin_jk &bjk,data_info &data_inf);

void bins_ini(bins &B,data_info &data_inf);

void calc_temp_ini(calc_temp &ct,data_info &data_inf);
