#include <omp.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define _USE_MATH_DEFINES
#include <time.h>
#include "data_def.h"
#include "calcs.h"
#include "sky_calcs.h"
#include "PB_calcs.h"
#include "initialization.h"
#include "outp.h"
#include "corels.h"
#include <algorithm>    // std::min

using namespace std;

const double pi=M_PI;
const double r2d=180.0/pi;

void jk_initilize(data_info &data_inf)
{
  data_inf.n_jk=data_inf.n_jk_regions;
  if (data_inf.do_jk==0){
    return;}
  cout<<"jk_ini: "<<data_inf.n_jk_regions<<endl;
  data_inf.jk_regions=new int [data_inf.n_jk_regions];
  for (int i=0;i<data_inf.n_jk_regions;i++)
    {
       data_inf.jk_regions[i]=i;
    }
  
  if (data_inf.do_jk>1)
  {
    data_inf.n_jk=data_inf.n_jk_regions*(data_inf.n_jk_regions+1)/2;
  }
}

void seed()
{
  srand(124693517);
}

double jk_prob_PB()
{
  return rand() / (RAND_MAX + 1.);
}

void no_jk(calc_temp &ct)
{;}//do nothing. jk_prob already set to 0 in intialization

void jk_check_sky(calc_temp &ct)
{
  if (ct.gjk[0]==ct.gjk[1])
    ct.jk_prob[1]=0;
}

/*void jk_check_sky_wl(calc_temp &ct)
{
  ct.jk_prob[1]=0;
  }*/

void jk_check_PB(calc_temp &ct)
{
  if (ct.gjk[0]==ct.gjk[1])
    {
      ct.jk_prob[0]=1;
      ct.jk_prob[1]=0;
    }
}


void final_jk_bins(bin_jk &bjk, data_info &data_inf) //subtract jk region from full to get jk sample
{
  bool sub=true;
  for(int i=0;i<data_inf.n_jk;i++)
    {
      add_bins(bjk.b,bjk.jk[i].b,bjk.jk[i].b,data_inf,sub);
    }
}


void bin_jk_nd_ini(bin_jk &bjk,data_info &data_inf,int &njk_samples){
  bjk.b=new bin [data_inf.n_bins];
  bjk.jk=new jackknife[njk_samples];

  if(data_inf.lin_bin) // linear bins intialisation
    bin_lin(bjk.b,data_inf);
  else
    bin_log(bjk.b,data_inf);

    for (int i=0;i<njk_samples;i++)
      {
        //bjk.jk[i].region=i;
        bjk.jk[i].b=new bin[data_inf.n_bins];

        if (data_inf.lin_bin)
  	     bin_lin(bjk.jk[i].b,data_inf);
        else bin_log(bjk.jk[i].b,data_inf);
      }
}


void bins_ini_nd(bins &B,data_info &data_inf,int &njk_samples)
{
  if (data_inf.do_auto_corr){
    bin_jk_nd_ini(B.SS,data_inf,njk_samples);
    bin_jk_nd_ini(B.final_auto_jk,data_inf,njk_samples);
    bin_jk_nd_ini(B.final_auto,data_inf,njk_samples);
  }
  bin_jk_nd_ini(B.SD,data_inf,njk_samples);
  bin_jk_nd_ini(B.SRs,data_inf,njk_samples);
  bin_jk_nd_ini(B.SRd,data_inf,njk_samples);
  bin_jk_nd_ini(B.DRs,data_inf,njk_samples);
  bin_jk_nd_ini(B.RsRd,data_inf,njk_samples);
  bin_jk_nd_ini(B.RsRs,data_inf,njk_samples);

  bin_jk_nd_ini(B.final_cross_jk,data_inf,njk_samples);
  bin_jk_nd_ini(B.final_cross,data_inf,njk_samples);
}


void final_jk_bins_nd(bin_jk &bjk, data_info &data_inf)
{
  int n_jk_samples=1;

  for (int i=0;i<data_inf.do_jk;i++)
    n_jk_samples*=(data_inf.n_jk_regions-i);

  for (int i=1;i<data_inf.do_jk;i++)
    n_jk_samples/=i;
  bins bjk_nd;
  bins_ini_nd(bjk_nd,data_inf,n_jk_samples);

}

void corel_jk1(data_info &data_inf,calc_temp &ct, int i)
{
  //  cout<<"corel_jk1"<<ct.n1<<" "<<i<<" "<<ct.gjk[0]<<"  "<<ct.gjk[1]<<endl;
      ct.bjk.jk[i].b[ct.n1].num++;
      //ct.bjk.jk[i].b[ct.n1].wt_num+=ct.wt_f.val;

      //ct.bjk.jk[i].b[ct.n1].wt_num.err+=ct.wt_f.err;

      ct.bjk.jk[i].b[ct.n1].p_bin[ct.n2].num++;

      ct.bjk.jk[i].b[ct.n1].p_bin[ct.n2].wt_num+=ct.wt_f.val;
      //ct.bjk.jk[i].b[ct.n1].p_bin[ct.n2].wt_num.err+=ct.wt_f.err;
      if (data_inf.which_corr==0)
        return ;

      ct.bjk.jk[i].b[ct.n1].p_bin[ct.n2].data[1].val+=ct.ep.val; //e+
      //      ct.bjk.jk[i].b[ct.n1].p_bin[ct.n2].data[1].err+=ct.ep.err; //e+
      ct.bjk.jk[i].b[ct.n1].p_bin[ct.n2].data[2].val+=ct.ex.val; //eX

      //      ct.bjk.jk[i].b[ct.n1].p_bin[ct.n2].data[2].err+=ct.ex.err; //eX
      ct.bjk.jk[i].b[ct.n1].p_bin[ct.n2].data[3].val+=ct.epp.val; //e+

      //      ct.bjk.jk[i].b[ct.n1].p_bin[ct.n2].data[3].err+=ct.epp.err; //e+

      ct.bjk.jk[i].b[ct.n1].p_bin[ct.n2].data[4].val+=ct.exx.val; //eX

      //      ct.bjk.jk[i].b[ct.n1].p_bin[ct.n2].data[4].err+=ct.exx.err; //eX
      ct.bjk.jk[i].b[ct.n1].p_bin[ct.n2].data[5].val+=ct.epx.val; //e+

      //      ct.bjk.jk[i].b[ct.n1].p_bin[ct.n2].data[5].err+=ct.epx.err;

      ct.bjk.jk[i].b[ct.n1].p_bin[ct.n2].data[6].val+=ct.etheta.val; //eX
      //      ct.bjk.jk[i].b[ct.n1].p_bin[ct.n2].data[6].err+=ct.etheta.err; //eX
}

void indx_jk2(data_info &data_inf,int ij[2],int k){//XXX Need to verify this
  ij[0]=data_inf.n_jk_regions+floor(0.5-0.5*(sqrt(pow(2*data_inf.n_jk_regions+1,2)-8*k)));
  ij[1]=k-data_inf.n_jk_regions*ij[0]+ij[0]+ij[0]*(ij[0]-1.)/2;
  //
  // ij[0] = data_inf.n_jk_regions - 2 - floor(sqrt(-8*k +
  //                             4*data_inf.n_jk_regions*(data_inf.n_jk_regions-1)-7)/2.0 - 0.5);
  // ij[1] = k + ij[0] + 1 - data_inf.n_jk_regions*(data_inf.n_jk_regions-1)/2 +
  //                               (data_inf.n_jk_regions-ij[0])*((data_inf.n_jk_regions-ij[0])-1)/2;
}


int jk2_indx(data_info &data_inf,int ij[2])
{
  // int i=ij[0];
  // int j=ij[1];
  // if (i>j){
  //   i=ij[1];
  //   j=ij[0];
  // }
  // return data_inf.n_jk_regions*i-(i-1)*i/2+j-i;
  if (ij[0]>ij[1]){
    return data_inf.n_jk_regions*ij[1]-(ij[1]-1)*ij[1]/2+ij[0]-ij[1];
  }
  else{
    return data_inf.n_jk_regions*ij[0]-(ij[0]-1)*ij[0]/2+ij[1]-ij[0];
  }
}

int corel_jk(gal &g1,gal &g2,data_info &data_inf,calc_temp &ct) // correlations b/w 2 galaxies
{
  ct.fail=ct.corel(g1,g2,data_inf,ct);

  if (ct.fail!=0 || data_inf.do_jk==0)
    {
      return ct.fail;
    }

    if (data_inf.do_jk==1){
      ct.jk_check(ct);
      for (int i=0;i<2;i++){
        if(ct.jk_prob[i]<0.5)
          continue;
        corel_jk1(data_inf,ct,ct.gjk[i]);
      }
    }
    else{
      int i=jk2_indx(data_inf,ct.gjk);
      corel_jk1(data_inf,ct,i);
    }
    return 0;
}
