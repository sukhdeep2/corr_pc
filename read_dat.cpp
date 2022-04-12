#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include "data_def.h"
#include "calcs.h"
using namespace std;


//==================================================================================================================

int read_data(gal g[],data_info &data_inf) // read data sets from files
{
  string p_f,z_f,e_f,w_f,jk_f,patch_f;
  int n_gal=0,it=0;
  double *wt,*wt_err;
  int *PB_starts;
  bool cross_shape_shape=0,read_phi=0,read_kappa=0,read_e3=0,read_e=1,read_jk=1,read_z=1,rand_phi=0; //read_e3 for ED
  double jk_prob=0,temp_en=0;

  if (data_inf.coordinates==4||data_inf.coordinates==5)
    read_z=0;

  switch(data_inf.which_corr)
    {
    case 0:
      read_e=0;
      break;
    case 3:
      cross_shape_shape=0;
      read_e=0;
      read_kappa=1;
      //if (data_inf.do_SS) //shape has kappa
      // read_e=0;
      break;
    case 4:
      read_e=0;
      if (data_inf.do_SS) //shape sample has e, desnity has kappa
	      read_e=1;
      else
	      read_kappa=1;
      cross_shape_shape=1;
      break;
    case 5:
      read_e=0;
      if (data_inf.do_SS)
        read_e3=1;
      cross_shape_shape=0;
    case 8:
      read_kappa=0;
      read_e=0;
      cross_shape_shape=1;
      break;
    case 9:
      read_kappa=1;
      read_e=0;
      cross_shape_shape=1;
      break;
    case 2:
      cross_shape_shape=1;
      break;
    case 7:
      read_e=0;
      break;
    default:
      cross_shape_shape=0;
    }
  //  cout<<"which_corr:"<<data_inf.which_corr<<" read_e:"<<read_e<<endl;
  if(data_inf.do_SS)
    {
      p_f=data_inf.shape_pos_file;
      patch_f=data_inf.shape_patch_file;
      z_f=data_inf.shape_z_file;
      e_f=data_inf.shape_e_file;
      w_f=data_inf.shape_wt_file;
      jk_f=data_inf.shape_jk_file;
      n_gal=data_inf.n_shape;
      wt=&data_inf.S_wt.val;
      wt_err=&data_inf.S_wt.err;
      PB_starts=&data_inf.PB_sorted_starts_S;
      if (data_inf.which_corr==4)
    	{
	      read_e=1;read_kappa=0;
	    }
      if (data_inf.coordinates==2||data_inf.coordinates==4)
	      read_e=1;
    }

  else if(data_inf.do_DD)
    {
      p_f=data_inf.density_pos_file;
      z_f=data_inf.density_z_file;
      w_f=data_inf.density_wt_file;
      jk_f=data_inf.density_jk_file;
      n_gal=data_inf.n_density;
      wt=&data_inf.D_wt.val;
      wt_err=&data_inf.D_wt.err;
      PB_starts=&data_inf.PB_sorted_starts_D;
      if(cross_shape_shape)
	    {
	      e_f=data_inf.density_e_file;
	    }
        else
	    {
	      read_e=0;read_kappa=0;read_phi=0;read_e3=0;
	    }
    }

  else if(data_inf.do_RsRs)
    {
      p_f=data_inf.shape_rand_pos_file;
      patch_f=data_inf.shape_rand_patch_file;
      z_f=data_inf.shape_rand_z_file;
      w_f=data_inf.shape_rand_wt_file;
      jk_f=data_inf.shape_rand_jk_file;
      wt=&data_inf.Rs_wt.val;
      wt_err=&data_inf.Rs_wt.err;
      PB_starts=&data_inf.PB_sorted_starts_Rs;
      n_gal=data_inf.n_shape_rand;
      read_e=0;read_kappa=0;read_phi=0;
      if (data_inf.coordinates==2||data_inf.coordinates==4)
      {
        //read_e=1;
        //e_f=data_inf.shape_rand_e_file;
        rand_phi=1;
      }
    }

  else if(data_inf.do_RsRd)
    {
      p_f=data_inf.density_rand_pos_file;
      z_f=data_inf.density_rand_z_file;
      w_f=data_inf.density_rand_wt_file;
      jk_f=data_inf.density_rand_jk_file;
      wt=&data_inf.Rd_wt.val;
      wt_err=&data_inf.Rd_wt.err;
      PB_starts=&data_inf.PB_sorted_starts_Rs;
      n_gal=data_inf.n_density_rand;
      read_e=0;read_kappa=0;read_phi=0;
      if ((data_inf.coordinates==2||data_inf.coordinates==4)&&data_inf.which_corr==2)
        {
	        rand_phi=1;
	      }
    }

  if (n_gal==0)
    {
      cout<<"n gal zero. Nothing to read: "<<p_f<<endl;
      return 0;
    }

  ifstream data_p,data_z,data_e,data_wt,data_jk,data_patch;
  data_p.open(p_f.c_str());// RA-Dec. sky positions
  if (read_z)data_z.open(z_f.c_str()); // redshift
  if(jk_f!="0")data_jk.open(jk_f.c_str()); //jackknife
  else {read_jk=0; cout<<"no jk file given"<<endl;}

  if(!data_inf.do_jk){
    read_jk=0;}
  
  if(w_f!="0")
    {
      cout<<"doing weighted calcs"<<endl;
      data_wt.open(w_f.c_str()); // weights
      if(!data_wt.is_open())
	{
	  cout<<"wt file not open: "<<w_f<<endl;
	  return 1;
	}
    }

  if(!data_p.is_open()){cout<<"initialisation::pos File not open: "<<p_f<<endl;return 1;}
  if(read_z&&!data_z.is_open()){cout<<"initialisation::z File not open:"<<z_f<<endl;return 1;}

  if(data_inf.do_jk && (!data_jk.is_open()) && read_jk){cout<<"initialisation::jk File not open: "<<w_f<<data_inf.do_SS<<data_inf.do_DD<<data_inf.do_jk<<endl;return 1;}

  if (data_inf.do_jk && !read_jk)
    jk_prob=-1;


  if (read_e||read_kappa||read_phi||read_e3)
    {
      data_e.open(e_f.c_str()); // ellipticities
      if(!data_e.is_open()){cout<<"initialisation:: e File not open:"<<e_f<<endl;return 1;}
      if (e_f.find("phi") != string::npos||read_phi)
	{
	  cout<<"Found phi in file name.. Will read e file as 1 column and store in phi.."<<e_f<<endl;
	  read_phi=1;
	  read_e=0;
	  read_kappa=0;
	}
      if (e_f.find("kappa") != string::npos||read_kappa) //change this to kappa.dat
        {
          cout<<"Found kappa in file name.. Will read e file as 1 column and store in e[0].."<<e_f<<endl;
          read_kappa=1;
	  read_e=0;
	  read_phi=0;
        }
    }

  if(data_inf.which_corr>=7) //need patches for lensing
    {
      data_patch.open(patch_f.c_str()); // patches
      if(!data_patch.is_open()){cout<<"initialisation:: patch File not open:"<<patch_f<<endl; return 1;}
    }

  string junk;

  int i=0;
  while(data_p.good() && (data_z.good()||!read_z) && (data_jk.good()|| !read_jk) && (data_wt.good()||w_f=="0") &&(data_patch.good()|| data_inf.which_corr<7) && (data_e.good()||!read_e||!read_kappa||!read_phi) && i<n_gal)
    {
      //			g[i].e=new double[2];
      data_p>>g[i].ra.val_deg>>g[i].dec.val_deg;

      g[i].ra.val_rad=angle_deg_to_rad(g[i].ra.val_deg);
      g[i].dec.val_rad=angle_deg_to_rad(g[i].dec.val_deg);

      g[i].redshift=0;

      if (read_z)
	{
	  data_z>>g[i].redshift;//>>g[i].zl_68>>g[i].zh_68;  //values from files
	  if (g[i].redshift<0)
	    {
	      g[i].redshift=0;
	      cout<<"gal redshift <0. changed to z=0"<<endl;
	    }
	  if(g[i].redshift>data_inf.z_max){cout<<"initialisation:: wrong z at data "<<i<<"     z=="<<g[i].redshift<<"  vs  "<<data_inf.z_max<<endl;return 1;}
	}

      g[i].wt=1;
      g[i].lens_wt=0;
      g[i].jk_prob=1;//jk_prob

      if (data_inf.periodic_box)
      {
        g[i].jk_prob=rand() / (RAND_MAX + 1.);
        g[i].DC=g[i].redshift;
        g[i].DA=0;
      }
      else
      {
        it=floor(g[i].redshift/data_inf.dz);
        g[i].DC=data_inf.dC_z[it];
        g[i].DA=data_inf.dA_z[it];
      }

      if (data_inf.do_jk==0||data_inf.n_jk==0)
	      g[i].jk_prob=0;

      if (w_f!="0")
      	data_wt>>g[i].wt;

      if(read_jk)
      {
        data_jk>>g[i].jk;
        //g[i].jk_prob=1;
      }
      else g[i].jk=-1;

      // g[i].e[0]=0;
      // g[i].e[1]=0;
      g[i].phi=-10;

      g[i].patch=-1;
      if (data_inf.which_corr>=7)
	      data_patch>>g[i].patch;

      if(read_e && !read_kappa && !read_phi)
      {
        g[i].e=new double [2];
        data_e>>g[i].e[0]>>g[i].e[1];
        elip_angle(g[i]);
      }

      if(read_e3)
      {
        g[i].e=new double [3];
        data_e>>g[i].e[0]>>g[i].e[1]>>g[i].e[2];
        temp_en=sqrt(g[i].e[0]*g[i].e[0]+g[i].e[1]*g[i].e[1]+g[i].e[2]*g[i].e[2]);
        if (temp_en>0)
        {
          g[i].e[0]/=temp_en;
          g[i].e[1]/=temp_en;
          g[i].e[2]/=temp_en;
        }
      }

      if(!read_e && read_phi &&!read_kappa)
        {
          data_e>>g[i].phi;
        }

      if(read_kappa && !read_e && !read_phi)
        {
	  g[i].e=new double [1];
          data_e>>g[i].e[0];
        }
      if (rand_phi)
        {
          g[i].phi=rand() / (RAND_MAX + 1.)*M_PI/2-M_PI/2;//g.phi is computed in range [-pi/2,pi/2]
        }

      if(data_inf.periodic_box && data_inf.data_sorted)
        {
          if (*PB_starts<0)
            {
              if(g[i].redshift-data_inf.periodic_box_size > data_inf.p_min)
          *PB_starts=i;
            }
        }
      g[i].include_prob=0;//galaxies with prob greater than some limit are dropped.. only for random auto correlations
      if (data_inf.do_RsRs||data_inf.do_RsRd)
      	g[i].include_prob=rand() / (RAND_MAX + 1.);

      *wt+=g[i].wt;
      *wt_err+=g[i].wt*g[i].wt;
      i++;
    }

  if(i<n_gal)
    {
      cout<<"gal data reading ended unexpectedly: "<<i<<"   "<<n_gal<<endl;
      if(!data_p.good()){cout<<"initialisation::pos File not open: "<<p_f<<endl;}
      if(!data_z.good()&&read_z){cout<<"initialisation::z File not open:"<<z_f<<endl;}
      if(data_inf.do_jk && !data_jk.good()){cout<<"initialisation::jk File not open: "<<jk_f<<endl;}
      if(!data_wt.good()){cout<<"initialisation::wt File not open:"<<w_f<<endl;}
      if(!data_e.good() && (read_e||read_phi||read_kappa)){cout<<"initialisation::e File not open:"<<w_f<<endl;}
      //      cout<<"do_DD= "<<data_inf.do_DD<<"   do_RR= "<<data_inf.do_RR<<"   do_SS= "<<data_inf.do_SS<<endl;
      return 1;
    }

  data_p.close();
  if (read_z)
    data_z.close();
  if(w_f!="0")data_wt.close();
  if(read_e||read_kappa||read_phi)
    data_e.close();

  return 0;
}

int read_patch(int patch_id, int &patch_gal_size_max, gal* &patch_gal, data_info data_inf) // read data sets from files
{
  int n=patch_id;
  string filename;
  stringstream ss;
  for (int i=0;i<3;i++)
    {
      ss << n%10;
      n/=10;
    }
  string str = ss.str();
  std::reverse(str.begin(),str.end());
  //filename=""+data_inf.patch_file+str+".tbl";
  filename="/home/rmandelb.proj/data-shared/weaklens-010-goodsigma/largecat/bgCat-010-"+str+".tbl";

  int n_patch_gal=0,ni=0,it=0;
  double ra,dec,z,e[2],e_err,tempelate,z_ll,z_ul,magr,junk2;
  string junk,fline;
  int gal_counter=0;

  double R=0.87;
  double iSN=0.36;
  double iSN_sq=pow(iSN,2.0);

  ifstream data;
  data.open(filename.c_str());//
  if(!data.is_open()){cout<<"patch File not open: "<<filename<<endl;return 1;}

  data>>n_patch_gal;
  getline(data,fline);//need to move pointer to next line

  if (patch_gal_size_max<n_patch_gal)
    {
      patch_gal=new gal[n_patch_gal];
      cout<<"reading patch"<<patch_id<<"  patch size not enough for "<<n_patch_gal<<endl;
    }
  //  cout<<"reading patch"<<patch_id<<"  "<<n_patch_gal<<endl;
  for (int i=0;i<n_patch_gal;i++)
    {
      data>>junk2>>junk2>>ra>>dec>>e[0]>>e[1]>>e_err>>z;
      getline(data,fline);
      /*std::istringstream line(fline);

      line>>junk2>>junk2>>ra>>dec>>e[0]>>e[1]>>e_err>>z;//>>tempelate>>z_ll>>z_ul>>junk>>junk>>junk>>junk>>magr;*/ //this is very slow

      if (z>data_inf.z_max||z<data_inf.z_min)
	{
	  i-=1;
	  n_patch_gal-=1;
	  continue;
	}
      e_err*=2.0;
      e[0]*=1.0;
      //      if(ra<0||ra>360) cout<<"RA problem"<<endl;
      // if(dec<-90||dec>90) cout<<"DEC problem"<<endl;
      it=floor(z/data_inf.dz);
      patch_gal[gal_counter].DC=data_inf.dC_z[it];
      patch_gal[gal_counter].DA=data_inf.dA_z[it];

      patch_gal[gal_counter].redshift=z;
      patch_gal[gal_counter].ra.val_deg=ra;
      patch_gal[gal_counter].dec.val_deg=dec;
      patch_gal[gal_counter].ra.val_rad=angle_deg_to_rad(ra);
      patch_gal[gal_counter].dec.val_rad=angle_deg_to_rad(dec);

      patch_gal[gal_counter].e=new double [2];
      patch_gal[gal_counter].e[0]=e[0]/(2.0*R);
      patch_gal[gal_counter].e[1]=e[1]/(2.0*R);
      patch_gal[gal_counter].wt=1.0/(pow(e_err,2.0)+iSN_sq);
      patch_gal[gal_counter].jk=-1;
      patch_gal[gal_counter].jk_prob=0;
      patch_gal[gal_counter].include_prob=0;//galaxies with prob greater than some limit are dropped.. only for random auto correlation
      patch_gal[gal_counter].lens_wt=0;
      gal_counter++;
    }
  std::ostringstream s;
  s<<patch_id;
  string out_file=data_inf.out_file+"source_gal1_"+s.str()+".dat";
  //  outp_gal(patch_gal,n_patch_gal,data_inf,out_file);

  return n_patch_gal;
}
