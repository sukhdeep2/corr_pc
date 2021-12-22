//define the data types used

#ifndef data_def
#define data_def

#include <string>
using namespace std;

struct var_st{
  double val,err;
};

struct angl{    //angle definition
  /*  int deg;
  int min;
  int sec;*/

  double val_rad;       //value calc in degrees
  double val_deg;
};


struct gal{     // each galaxy

  int jk,patch;

  angl ra; //position in sky
  angl dec;

  //  angl orient;      //orientation angle of major axis

  double *e,phi;     //elipticity
  //  double e0;        //ellipticity defined along major axis

  //  var_st redshift; //**********************change to double

  double redshift,DC,DA; //comoving distances

  //  double zl_68,zh_68;

  double wt,jk_prob,lens_wt;//weight assigned to each galaxy

  double include_prob;//probability that galaxy should be included in pair counting.. used to subsample randoms for RR term
};

struct data_info{ // to store information about data set, file names, etc.
  bool do_auto_corr,do_cross_corr,do_auto_rebin,use_comoving,RsRd_same,periodic_box,rand_subsample;
  int do_jk;
  int data_sorted;//0: no sorting, 1:sorted by z, 2:sorted by dec... for wl sorting by patch, always
  int which_corr; /*0: density-density, 1: ia: shape-density, 2 ia:shape-shape,
		    3:kappa-density 4: ia: shape-kappa 5: ia: ED
		    7: wl:source shape-lens, 8: wl:source shape-lens shape,
		    9: wl: shape-convergence, 10: wl: density-density*/

  int coordinates;// 0: rp-pi, 1:r-mu, 2: rp-phi 3: rp, 4: theta-phi, 5:theta  6: xyz, rp-pi  7: xyz, r-mu

  int sig_crit; //0: no sig crit; 1: calc sig crit for each pair, 2: use g1 wt as sig crit

  int estimator; //0:landy-szalay (LS): cross, 1: LS:auto, 2:LS:auto+rebin, 3: LS:auto+cross+auto_rebin, 4:mean

  string shape_pos_file,shape_z_file,shape_e_file,shape_wt_file,shape_jk_file,shape_patch_file;
  string density_pos_file,density_z_file,density_wt_file,density_jk_file,density_e_file;
  string shape_rand_pos_file,shape_rand_z_file,out_file,shape_rand_wt_file,shape_rand_jk_file,shape_rand_patch_file;//,shape_rand_e_file;
  string density_rand_pos_file,density_rand_z_file,density_rand_wt_file,density_rand_jk_file;
  //  int (*corel_func)(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct);
  string jk_file,patches_file;

  string distance_file;

  int n_bins,n_threads,n_jk,n_jk_regions,n_bin_var,n_patch;	//number of bins, omp threads
  int *jk_regions,*patches;
  double bin_r_min,bin_r_max,r_size;

  bool lin_bin;		//if true, use linear bins. otherwise log bins

  //  double n_rand_ratio; //ration of n_rand/n_gal

  double dz,z_min,z_max,z_sep_max,z_sep_min; //z_sep_max=max redshift separation (min in -ve). z_max=max redshift. dz=increment in z for redshift distance (DA(z) and H(z)) calculations
  double *dA_z,*H_z,*dC_z;

  int n_p_bin,PB_sorted_starts_S,PB_sorted_starts_D,PB_sorted_starts_Rs;   //number of p_bin
  double p_max,p_min,p_size;
  double periodic_box_size;
  angl da_max,da_min;
  double CMB_DIST;
  int n_shape,n_density,n_density_rand,n_shape_rand;	//number of galaxies
  int g_with_shape;//when doing shape X density, which galaxy has shape.. for lensing always g2.. for IA like, 1 if shape sample has more galaxies otherwise 2

  //double ra_max,ra_min,dec_max,dec_min;	//size of the grid
  //double rang_taken;	//area taken within the grid

  bool do_DD,do_SS,do_DRs,do_SRs,do_SD,do_RsRd,do_DRd,do_SRd,do_RsRs; //so program knows which corr it is doing

  var_st DD_wt,RsRs_wt,DRs_wt,DRd,SS_wt,SD_wt,SRs_wt,SRd_wt,RsRd_wt,D_wt,S_wt,Rs_wt,Rd_wt;
  double include_prob_RR;
};

struct P_bin{
  double p_min;
  double p_max;
  long num; //number of pairs in the bin
  double wt_num;
  // var_st wt_num;
  // double num_den;		//pairs /binsize
  var_st data[7];	//0=W,1=e+,2=eX,3=e++,4=eXX,5=e+X,6=theta
  double ls_err;
};

struct bin{                    //bin info
  double b_min;
  double b_max;
  long num;//number of pairs in the bin
  double wt_num;
  // var_st wt_num;
  // double num_den;		//pairs /binsize(area)
  // double b_err; //poisson error in pair count
  double sig_crit; //weighted sigma_crit when doing lensing
  P_bin *p_bin;
  var_st data[7];	//0=W,1=e+,2=eX,3=e++,4=eXX,5=e+X,6=angle
  double ls_err_int;
};

struct jackknife{
  //int region;// region excluded
  bin *b;
};


struct bin_jk{
  bin *b;
  jackknife *jk;
};

struct bins{
  bin_jk SS,RsRs,SRs; // bins
  bin_jk SD,DRs,RsRd,SRd;
  bin_jk final_auto,final_cross;   // bin for final output
  bin_jk final_auto_jk,final_cross_jk;   // bin for final output
};

struct calc_temp{
  bin_jk bjk;
  var_st ep,ex,epp,exx,epx,etheta,wt_f;

  int n1,n2,gjk[2],n_patch_gal,ret_z,ret_dec;
  double delta_z,da,da_max,z_avg,d_P,da_R,theta,sig_crit,sig_crit_inv,dl,ds,dls;
  double z_max,z_min,ra_max,ra_min,dec_max,dec_min,ra_lim,include_prob_max;

  double e_new_i[2];double e_new_j[2];double temp[2];double jk_prob[2];
  int PBx,PBy,PBz,fail;
  void (*pi_calc)(gal &g1,gal &g2,calc_temp &ct, data_info &data_inf);
  void (*rp_calc)(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct);
  void (*r_calc)(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct);
  void (*jk_check)(calc_temp &ct);
  double (*angle_orient_axis)(gal &g1,gal &g2);
  int (*corel)(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct);
  void (*sig_crit_func)(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct);
};

#endif
