int outp_gal(gal g[],data_info data_inf,string FileName);	//output galaxy data set to file

int outp_bins(bin b[],data_info data_inf, string FileName);	//output to file

int outp_temp_bins(bin b[],data_info data_inf, string FileName);    //output to file                                                                 

int outp_data_info(data_info data_inf,string FileName);
      
int inp_bins(bin b[],data_info data_inf, string FileName);    //output to file                                                                       

int outp_2Dbins(bin b[],data_info data_inf, string FileName);

int add_bins(bin b1[],bin b2[],bin final[], data_info data_inf,bool subtract);     //add two bins


