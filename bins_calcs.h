#include <omp.h>
#include "data_def.h"

using namespace std;

void W_calc_cross(bin sd[], bin s_r[], data_info data_inf); //W_delta calculation. This measurement is independent of p bins

void bin_calc_cross(bin sd[],bin sr[],bin dr[],bin rr[],data_info data_inf);// get binned data in order

void final_calc(bin sd[], bin sr[],bin final[], data_info &data_inf);// SD/RD-1

void final_calc(bin sd[], bin sr[],bin dr[],bin rr[], bin final[], data_info &data_inf);// calculations to get final answers

void final_calc_auto(bin ss[], bin sr[],bin rr[], bin final[], data_info &data_inf);// calculations to get final answers

void final_calc_bins(bin sd[], bin sr[],bin dr[],bin rr[], bin final[], data_info &data_inf);

void int_calc(bin dat[],data_info data_inf);// integrate over p bins

void int_calc(bin sd[], bin sr[],bin dr[],bin rr[], bin final[], data_info &data_inf);// for multiple data sets. Calls other int_calc for each data set.

void int_calc(bin ss[], bin sr[],bin rr[], bin final[], data_info &data_inf);

void int_calc(bin ss[], bin sr[], data_info &data_inf);

void jk_final(bin bins_final_jk[] ,bin bins_final[],jackknife jk_final[],data_info data_inf);

void jk_final2D(bin bins_final_jk[] ,bin bins_final[],jackknife jk_final[],data_info data_inf);

void re_bin(bin ss[], data_info &data_inf);

void final_calc_mean_projected(bin sd[],bin sr[],bin final[],data_info data_inf);
