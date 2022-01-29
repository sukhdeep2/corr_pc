//calculating correlations

 /*
 How it all works:

 function corels1 calculate DD or RR, avoids double counting. Corels2 count DR. Finally give sum(x) in each bin and sum(x^2) in each error bin which is sorted later on.

 corel called from corels1 or 2 just does calculations for a pair of galaxies. measures their angular diameter distance(dA_r) and average line of sight distance(dP). Then calculates e+ eX along line joining the two galaxies and e+^2,e++ etc. from that and bins all the quantities in dA_r and dP.


 */

#include <omp.h>
#include <iostream>
#include <stdio.h>
#include <math.h>
#define _USE_MATH_DEFINES
#include <time.h>
#include <stdlib.h>
#include "data_def.h"
#include "calcs.h"
#include "corr.h"
#include "initialization.h"
#include "bins_calcs.h"
#include "outp.h"
#include "read_dat.h"
#include "sky_calcs.h"
#include "PB_calcs.h"
#include "jk.h"
using namespace std;

 //double pi=3.14159265359;
//extern const double pi;
const double r2d=180.0/M_PI;

int do_pair_unsorted(gal &g1,gal &g2,calc_temp &ct,data_info &data_inf, long &n_out,long &n_total, int &n_num)
{
  if (data_inf.coordinates!=4&&data_inf.coordinates!=5)
    {
      //leave these cuts in z due to use of z_sep_min/max..
      if((g2.redshift+data_inf.periodic_box_size*ct.PBz)>ct.z_max||
	 (g2.redshift+data_inf.periodic_box_size*ct.PBz)<ct.z_min)
	return 1;
    }

  if((g2.ra.val_deg+data_inf.periodic_box_size*ct.PBx)<ct.ra_min||
     (g2.ra.val_deg+data_inf.periodic_box_size*ct.PBx)>ct.ra_max||
    (g2.dec.val_deg+data_inf.periodic_box_size*ct.PBy)>ct.dec_max||
      (g2.dec.val_deg+data_inf.periodic_box_size*ct.PBy)<ct.dec_min)//addition doesnot add much overhead
   return 1;

  /*if(abs(g2.ra.val_deg+data_inf.periodic_box_size*ct.PBx-g1.ra.val_deg)<ct.ra_lim)
     return 1;//slower*/

  /*  if((g2.dec.val_deg+data_inf.periodic_box_size*ct.PBy)>ct.dec_max||
     (g2.dec.val_deg+data_inf.periodic_box_size*ct.PBy)<ct.dec_min)
    return 1;*/

  n_out+=corel_jk(g1,g2,data_inf,ct);//for shape density first g should be shape
  n_total++;
  return 1;
}

int one=1;
int inf=1.e9;

int do_pair_sorted(gal &g1,gal &g2,calc_temp &ct,data_info &data_inf, long &n_out,long &n_total, int &n_num)
{
  if (data_inf.coordinates!=4&&data_inf.coordinates!=5){
    //leave these cuts in z due to use of z_sep_min/max..
    if((g2.redshift+data_inf.periodic_box_size*ct.PBz)>ct.z_max||
       (g2.redshift+data_inf.periodic_box_size*ct.PBz)<ct.z_min)
      return ct.ret_z;//n_num; //breaks the while loop
  }
  if((g2.dec.val_deg+data_inf.periodic_box_size*ct.PBy)>ct.dec_max||
     (g2.dec.val_deg+data_inf.periodic_box_size*ct.PBy)<ct.dec_min)
    return ct.ret_dec;

  if((g2.ra.val_deg+data_inf.periodic_box_size*ct.PBx)<ct.ra_min||
     (g2.ra.val_deg+data_inf.periodic_box_size*ct.PBx)>ct.ra_max)
    return 1;

  n_out+=corel_jk(g1,g2,data_inf,ct);
  n_total++;
  return 1;
}

void comb_thread_bins(bin_jk &bjk,data_info &data_inf,calc_temp &ct)
{
  bool sub=false;
  add_bins(bjk.b,ct.bjk.b,bjk.b,data_inf,sub);

  if (data_inf.do_jk!=0)
    {
      for(int i=0;i<data_inf.n_jk;i++)
        {
          add_bins(bjk.jk[i].b,ct.bjk.jk[i].b,bjk.jk[i].b,data_inf,sub);
        }
    }
}


long corels1(gal g1[], bin_jk &bjk,data_info &data_inf)
//correlations for 1 data set
{
  long n_out=0,n_total=0,PB=0; //to keep track of pair counts, total and excluded (n_out)

  int n_num=0; // total number of galaxies. value updated below.

  double start=omp_get_wtime();// to get hang of how much time it took

  data_inf.g_with_shape=1;

  if (data_inf.do_SS)
    n_num=data_inf.n_shape; //num of galaxies in data

  else if(data_inf.do_RsRs)
    n_num=data_inf.n_shape_rand; //num of galaxies in shape randoms

  omp_set_dynamic(0); //Open Mp threads
  if(data_inf.n_threads==0||data_inf.n_threads>omp_get_num_procs()) //no input (or 0) from user about number of threads. Use all available cores
    omp_set_num_threads(omp_get_num_procs());
  else
    omp_set_num_threads(data_inf.n_threads); //Number of threads input by user

  int begin=0,end=n_num;

  int (*do_pair)(gal &g1,gal &g2,calc_temp &ct,data_info &data_inf, long &n_out, long &n_total, int &n_num);
  if(data_inf.data_sorted==0)
    do_pair=&do_pair_unsorted;
  else
    do_pair=&do_pair_sorted;

  void (*lim_calc)(gal &g1,calc_temp &ct,data_info &data_inf);
  if (data_inf.periodic_box)
    lim_calc=&lim_calc_PB;
  else if(data_inf.coordinates==6 && !data_inf.periodic_box)
    lim_calc=&lim_calc_xyz;
  /*  else if(data_inf.coordinates==2||data_inf.coordinates==3)
    lim_calc=&lim_calc_sky_projected;
  else if(data_inf.coordinates==4)
  lim_calc=&lim_calc_sky_angular;*/
  else
    lim_calc=&lim_calc_sky;

    cout<<"corels 1 "<<g1[n_num-1].redshift<<endl;
#pragma omp parallel reduction(+:n_out,n_total,PB) //private(temp,delta_z,z_avg,da_R,da,d_P,wt_f,theta,n1,n2,e_new_i,e_new_j) 	//parallel code!!
  {
    calc_temp ct;
    double xt=0,yt=0,zt=0;
    calc_temp_ini(ct,data_inf);
    int j=0;
    long n_out1=0,n_total1=0;
#pragma omp for schedule(dynamic,100)//collapse(2) //schedule(dynamic,10) nowait
    for (int i=begin;i<end;i++)
      {
	if (g1[i].include_prob>ct.include_prob_max)
	  continue;
	lim_calc(g1[i],ct,data_inf);
	xt=ct.PBx;//max(0,ct.PBx);
	yt=ct.PBy;//max(0,ct.PBy);
	zt=ct.PBz;//max(0,ct.PBz);

	ct.PBx=ct.PBx>0;
	for (int xi=0;xi<=abs(xt);xi++)
	  {
	    ct.PBy=yt>0;
	    for(int yi=0;yi<=abs(yt);yi++)
	      {
		ct.PBz=zt>0;
		for(int zi=0;zi<=abs(zt);zi++)
		  {
		    j=i+1; //This probably won't work for periodic box... This actually works... only that half of the pairs are counted... even in periodic box... Think pair distances when periodic boundary comes into play
		    //This doesnot work for large pi correlations

		    if(data_inf.data_sorted!=0)//(data_inf.periodic_box)
		      {
			if (ct.PBz>0)
			  j=0;
			else if (ct.PBz<0)
			  j=i+1;//data_inf.PB_sorted_starts;
		      }
		    PB+=(xi+yi+zi)>0;

		    while(j<n_num) //doing correlations, pair by pair. Note no double counting here
		      {
			if (g1[j].include_prob>ct.include_prob_max)
			  {j+=1;continue;}
			j+=do_pair(g1[i],g1[j],ct,data_inf,n_out1,n_total1,n_num);
		      }
		    ct.PBz--;
		  }
		ct.PBy--;
	      }
	    ct.PBx--;
	  }
      }

    #pragma omp critical
    {
      comb_thread_bins(bjk,data_inf,ct);
      n_out+=n_out1;
      n_total+=n_total1;
    }
  }

  if (data_inf.do_jk==1)
    {
      final_jk_bins(bjk,data_inf);
    }

  long n_in=n_total-n_out;

  cout<<"corel 1 ends: time:" <<omp_get_wtime()-start<<"  total pair= "<<n_total<<"	"<<"pairs out= "<<n_out<<"          pairs in="<<n_in<<"   PB="<<PB<<endl;
  return n_in;
}


long corels2(gal g1[],gal g2[], bin_jk &bjk,data_info &data_inf)
// cross correlations b/w 2 data sets
{
  long n_out=0,n_total=0, PB=0;int n_threads; // to keep track of pair counts, total and excluded (n_out), number of threads

  int n_num1=0,n_num2=0,num_temp=0;// total number of galaxies for each data set. value updated below.
  bool order_changed=false;
  if (data_inf.do_DRs)
    {
      n_num2=data_inf.n_density;
      n_num1=data_inf.n_shape_rand;
      data_inf.g_with_shape=0;
      if (data_inf.coordinates==2||data_inf.coordinates==4)
	data_inf.g_with_shape=1;
    }

  else if (data_inf.do_RsRd)
    {
      n_num2=data_inf.n_shape_rand;
      n_num1=data_inf.n_density_rand;
      data_inf.g_with_shape=0;
      if (data_inf.coordinates==2||data_inf.coordinates==4)
        data_inf.g_with_shape=2;
    }
  else if(data_inf.do_RsRs)
    {
      n_num1=data_inf.n_shape_rand;
      n_num2=n_num1;
      data_inf.g_with_shape=0;
      if (data_inf.coordinates==2||data_inf.coordinates==4)
        data_inf.g_with_shape=1;
    }

  else if (data_inf.do_SRd)
    {
      n_num2=data_inf.n_shape;
      n_num1=data_inf.n_density_rand;
      data_inf.g_with_shape=2;
    }

  else if (data_inf.do_SD)
    {
      n_num2=data_inf.n_shape;
      n_num1=data_inf.n_density;
      data_inf.g_with_shape=2;
    }

  else if (data_inf.do_SRs)
    {
      n_num2=data_inf.n_shape;
      n_num1=data_inf.n_shape_rand;
      data_inf.g_with_shape=2;
    }

  if(n_num1==0 || n_num2==0)
    {
      cout<<"corels2, ngal=0"<<endl;
      return 0;
    }

  gal *g1_orig,*g2_orig;

  if (n_num2>n_num1&&!(data_inf.coordinates==3 && data_inf.sig_crit!=0))
    {//2nd condition excludes the wl case
      g1_orig=g1;
      g2_orig=g2;
      g1=g2;
      g2=g1_orig;
      order_changed=true;
      num_temp=n_num1;
      n_num1=n_num2;
      n_num2=num_temp;
      if(data_inf.g_with_shape==2)
	data_inf.g_with_shape=1;
      else if(data_inf.g_with_shape==1)
	data_inf.g_with_shape=2;
      cout<<"corels2, n2>n1, order changed"<<endl;
    }
  cout<<"corels2, n1, n2: "<<n_num1<<" "<<n_num2<<endl;
  omp_set_dynamic(0);
  if(data_inf.n_threads==0||data_inf.n_threads>omp_get_num_procs()*1.5) //no input from user about number of threads. Use all available cores
    omp_set_num_threads(omp_get_num_procs());
  else
    omp_set_num_threads(data_inf.n_threads);

  double start=omp_get_wtime();// to get hang of how much time it took

  int (*do_pair)(gal &g1,gal &g2,calc_temp &ct,data_info &data_inf, long &n_out,long &n_total, int &n_num);
  if(data_inf.data_sorted==0)
    do_pair=&do_pair_unsorted;
  else
    do_pair=&do_pair_sorted;

  void (*lim_calc)(gal &g1,calc_temp &ct,data_info &data_inf);

  if (data_inf.periodic_box)
    lim_calc=&lim_calc_PB;
  else if(data_inf.coordinates==6 && !data_inf.periodic_box)
    lim_calc=&lim_calc_xyz;
  /*  else if(data_inf.coordinates==2||data_inf.coordinates==3)
    lim_calc=&lim_calc_sky_projected;
  else if(data_inf.coordinates==4)
  lim_calc=&lim_calc_sky_angular;*/
  else if(data_inf.coordinates==3&&data_inf.sig_crit!=0)
    lim_calc=&lim_calc_sky_wl;
  else
    lim_calc=&lim_calc_sky;

  long omp_num=1000;

#pragma omp parallel //reduction(+:n_out,n_total,PB)  //private(temp,delta_z,z_avg,da_R,da,d_P,wt_f,theta,n1,n2,e_new_i,e_new_j)   //parallel code!!
  {
    omp_num=floor(n_num1/(omp_get_num_threads()*40));
    if(omp_num<1)
      omp_num=1;
    //    cout<<"omp:"<<omp_get_num_threads()<<"  "<<omp_num<<endl;
    calc_temp ct;
    calc_temp_ini(ct,data_inf);
    int begin=0,end=n_num1;
    int j=0;
    double xt=0,yt=0,zt=0;
    long n_out1=0,n_total1=0;
    long starts=0,starts_PB=data_inf.PB_sorted_starts_D;
    if (data_inf.do_SRs||data_inf.do_RsRs||data_inf.do_RsRd)
      {
	      starts_PB=data_inf.PB_sorted_starts_Rs;
      }

    double dz=0;//(g1[begin].redshift.val-g2[starts].redshift.val);
    //    cout<<"corels2: starts_pb= "<<starts_PB<<endl;
    #pragma omp for schedule(dynamic,omp_num)//collapse(2) //schedule(dynamic,10) nowait
    for (int i=0;i<n_num1;i++)
          {
	    if (g1[i].include_prob>ct.include_prob_max)
	      {
		      continue;
	      }
	    lim_calc(g1[i],ct,data_inf);
	    xt=ct.PBx;//max(0,ct.PBx);
	    yt=ct.PBy;//max(0,ct.PBy);
	    zt=ct.PBz;//max(0,ct.PBz);

	    if (data_inf.data_sorted==1)
	      {
        dz=(g2[starts].redshift);
        while ((dz>ct.z_max || dz<ct.z_min)&&starts<n_num2)
          {
            starts++;
            dz=(g2[starts].redshift);
          }
	      }
	    if (data_inf.data_sorted==2)
        {
          dz=(g2[starts].dec.val_deg);
          while ((dz>ct.dec_max || dz<ct.dec_min)&&starts<n_num2)
            {
              starts++;
              dz=(g2[starts].dec.val_deg);
            }
        }

	    ct.PBy=yt>0;
	    for(int yi=0;yi<=abs(yt);yi++)
	      {
        ct.PBx=xt>0;
        for(int xi=0;xi<=abs(xt);xi++)
          {
            ct.PBz=zt>0;
            for (int zi=0;zi<=abs(zt);zi++)
              {
              if (ct.PBz==0)
                j=starts;

              else if (ct.PBz<0)
                {
                  //			    starts_PB=data_inf.PB_sorted_starts;
                  if (data_inf.data_sorted==1 && starts_PB<n_num2)
                    {
                    dz=(g2[starts_PB].redshift+ct.PBz*data_inf.periodic_box_size);
                    while ((dz>ct.z_max || dz<ct.z_min)&&starts_PB<(n_num2-1))
                      {
                        starts_PB++;
                        dz=(g2[starts_PB].redshift+ct.PBz*data_inf.periodic_box_size);
                      }
                      if (starts_PB<0)starts_PB=0;
                    }
			            j=starts_PB;
			          }
			        else
			          j=0;

			PB+=(xi+yi+zi)>0; //this is likely to be affected by race conditions.
			//double start2=omp_get_wtime();
			while(j<n_num2)
			  {
			    if (g2[j].include_prob>ct.include_prob_max)
			      {j+=1;continue;}
			    j+=do_pair(g1[i],g2[j],ct,data_inf,n_out1,n_total1,n_num2);
			  }
			ct.PBz--;
		  }
		  ct.PBx--;
		}
		ct.PBy--;
	      }
	  }

   #pragma omp critical
    {
      comb_thread_bins(bjk,data_inf,ct);
      n_out+=n_out1;
      n_total+=n_total1;
      // bool sub=false;
      // add_bins(bjk.b,ct.bjk.b,bjk.b,data_inf,sub);
      // n_out+=n_out1;
      // n_total+=n_total1;
      // if (data_inf.do_jk)
      //   {
      //     for(int i=0;i<data_inf.n_jk;i++)
      //       {
	    //          add_bins(bjk.jk[i].b,ct.bjk.jk[i].b,bjk.jk[i].b,data_inf,sub);
      //       }
      //   }
    }
  }

  if (data_inf.do_jk==1)
    {
      final_jk_bins(bjk,data_inf);
    }
  if (order_changed)
    {
      g1=g1_orig;
      g2=g2_orig;
    }

  long n_in=n_total-n_out;

  cout<<"corel2 ends::time: "<<omp_get_wtime()-start<<"  total pair= "<<n_total<<"	"<<"pairs out= "<<n_out<<"   pairs in="<<n_in<<"    PB="<<PB<<endl;
  return n_in;
}


void do_wl_patch(gal l[],gal patch_gal[],calc_temp &ct, data_info &data_inf,int &n_num,
		int &patch,int &n_patch_gal,
		int &start,int &end,long &n_out,long &n_total)
{
  for (int j=start;j<=end;j++)
    {
      if(l[j].patch!=patch)
	{
	  cout<<"wl_jk: patch problem "<<patch<<endl;
	  break;
	}

      int start2=0;
      lim_calc_sky_wl(l[j],ct,data_inf);

      while (start2<n_patch_gal)
	{
	  if (patch_gal[start2].dec.val_deg>=ct.dec_min)
	    break;
	  start2++;
	}
      //start2=0; gives same pairs for slect 42
      // int source_end=0;
      for (int k=start2;k<n_patch_gal;k++)
	{
	  //patch_gal[k].lens_wt=0;
	  if(patch_gal[k].dec.val_deg>ct.dec_max)
	    {
	      break; //works for select 42
	    }

	  //if(patch_gal[k].dec.val_deg<ct.dec_min)  //not necessary. same for select 42
	  //continue;
	  if(patch_gal[k].redshift<ct.z_min)
	    continue;
	  if(patch_gal[k].ra.val_deg<ct.ra_min)
	    continue;
	  if(patch_gal[k].ra.val_deg>ct.ra_max)
	    continue;

	  n_total++;
	  n_out+=corel_jk(l[j],patch_gal[k],data_inf,ct);
	  //      cout<<n_total<<endl;
	  //source_end=k;
	}
    }
}

long wl_corels(gal l[],gal R[], bin_jk &bjk,bin_jk &bjkR, data_info &data_inf)
{
  long n_out=0,n_total=0,nR_out=0,nR_total=0; //to keep track of pair counts, total and excluded (n_out)
  int n_num=0,nR_num=0; // total number of galaxies. value updated below.
  double start=omp_get_wtime();// to get hang of how much time it took

  //if (data_inf.do_SS)
  n_num=data_inf.n_shape; //num of galaxies in data

  //else if(data_inf.do_RsRs)
  nR_num=data_inf.n_shape_rand; //num of galaxies in shape randoms

  data_inf.g_with_shape=2;
  if(data_inf.which_corr=8)
    {
      if(data_inf.coordinates==2||data_inf.coordinates==4)
	data_inf.g_with_shape=1;
    }
  omp_set_dynamic(0);
  if(data_inf.n_threads==0||data_inf.n_threads>omp_get_num_procs()) //no input from user about number of threads. Use all available core
    {
      omp_set_num_threads(omp_get_num_procs());
      data_inf.n_threads=omp_get_num_procs();
    }
  else
    omp_set_num_threads(data_inf.n_threads);

  int omp_dynamic_num=1;/*data_inf.n_patch/data_inf.n_threads/4;
  if (omp_dynamic_num<1)
    omp_dynamic_num=1;
    cout<<"Openmp num:"<<omp_dynamic_num<<"  threads"<<data_inf.n_threads<<endl;*/

  int patch_gal_size_max=2000000;

#pragma omp parallel reduction(+:n_out,n_total)  //private(temp,delta_z,z_avg,da_R,da,d_P,wt_f,theta,n1,n2,e_new_i,e_new_j)   //parallel code!
  {
    calc_temp ct,ctR;
    calc_temp_ini(ct,data_inf);
    calc_temp_ini(ctR,data_inf);
    gal *patch_gal;
    patch_gal=new gal[patch_gal_size_max];
    if(data_inf.data_sorted!=0)
      {
	int start=0,end=n_num;
	int startR=0,endR=n_num;

#pragma omp for schedule(dynamic,omp_dynamic_num)//collapse(2) //schedule(dynamic,10) nowait
        for (int i=0;i<data_inf.n_patch;i++)
          {
	    double start2=omp_get_wtime();
	    int patch=data_inf.patches[i];
	    start=0;end=n_num-1;
            startR=0;endR=nR_num-1;
	    //gal *patch_gal;
	    int n_patch_gal=0;

	    while (start<n_num)
	      {
		if(l[start].patch==patch)break;
		start++;
	      }
	    while(end>=0)
	      {
		if(l[end].patch==patch)break;
		end--;
	      }

	    while (startR<nR_num)
              {
                if(R[startR].patch==patch)break;
                startR++;
              }
            while(endR>=0)
              {
                if(R[endR].patch==patch)break;
                endR--;
              }

	    if ((start>=n_num||end<start)&&(startR>=nR_num||endR<startR)) continue;

	    else
	      {
		n_patch_gal=read_patch(patch,patch_gal_size_max,patch_gal,data_inf);
		cout<<"read patch "<<patch<<"  with gals=  "<<n_patch_gal<<"  "<<omp_get_wtime()-start2<<endl;
		ct.n_patch_gal+=n_patch_gal;
		do_wl_patch(l,patch_gal,ct,data_inf,n_num,
			    patch,n_patch_gal,start,end,n_out,n_total);
		do_wl_patch(R,patch_gal,ctR,data_inf,nR_num,
                            patch,n_patch_gal,startR,endR,nR_out,nR_total);
		//cout<<"done patch "<<patch<<"  "<<omp_get_wtime()-start2<<endl;
		  }
	  }
      }
    delete[] patch_gal;
   #pragma omp critical
    {
      comb_thread_bins(bjk,data_inf,ct);
      comb_thread_bins(bjkR,data_inf,ctR);
      // bool sub=false;
      // add_bins(bjk.b,ct.bjk.b,bjk.b,data_inf,sub);
      // add_bins(bjkR.b,ctR.bjk.b,bjkR.b,data_inf,sub);
      // data_inf.n_density+=ct.n_patch_gal;
      // if (data_inf.do_jk)
      //   {
      //     for(int i=0;i<data_inf.n_jk;i++)
      //       {
      //         add_bins(bjk.jk[i].b,ct.bjk.jk[i].b,bjk.jk[i].b,data_inf,sub);
	    //         add_bins(bjkR.jk[i].b,ctR.bjk.jk[i].b,bjkR.jk[i].b,data_inf,sub);
      //       }
      //   }
    }
  }

  if(data_inf.do_jk==1)
    {
      final_jk_bins(bjk,data_inf);
      final_jk_bins(bjkR,data_inf);
    }

  long n_in=n_total-n_out;
  long nR_in=nR_total-nR_out;
  cout<<"wl 1 ends: time:" <<omp_get_wtime()-start<<"  total pair= "<<n_total<<"     "<<"pairs out= "<<n_out<<"          pairs in="<<n_in<<"  source files used:"<<data_inf.n_density<<endl;
  cout<<"wl 1 ends: time:" <<omp_get_wtime()-start<<"  total pair= "<<nR_total<<"     "<<"pairs out= "<<nR_out<<"          pairs in="<<nR_in<<"  source files used:"<<data_inf.n_density<<endl;
  return n_in;
}
