//some calculations: like angular distace b/w two galaxies.. same for sky and PB

double angle_rad(angl &a); // angel in radians from struct angl

double angle_deg_to_rad(angl &a); //Degree to radians conversion, for struct angle

float angle_deg_to_rad(float &a); // degree to radians, for angle in degree (not the struct)

double angle_deg_to_rad(double &a);

double absolute(double n);	//to return |n|

double round_dec(double &num, int dec_places);


double angle_orient_axis_PB(gal &g1,gal &g2);  //angle of line joining two galaxies with x-axis (RA)


void elip_rotate(double e_initial[], double e_final[], double &theta); //rotate ellipse in e_initial (2 ellipticites) by theta and store in e_final

void elip_rotate(double e_initial, double e_final[], double &theta); //rotate ellipse in e_initial (e0) by theta and store in e_final

double elip_angle(gal &g1);

void sig_crit_inv(calc_temp& ct,double &zl, data_info &data_inf);

double ephi_calc(gal &g1,calc_temp&ct);//return smallest angle between elipticity and line, in range [0,pi/2]
void pi_calc(gal &g1,gal &g2,calc_temp &ct, data_info &data_inf);
void pi_calc_phi(gal &g1,gal &g2,calc_temp &ct, data_info &data_inf);
void pi_calc_projected(gal &g1,gal &g2,calc_temp &ct, data_info &data_inf);
void r_mu_calc(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct);
void sig_crit_calc(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct);
void sig_crit_g1(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct);
