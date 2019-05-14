#include "HBHS.h"
#include "THDM.h"
#include "DecayTable.h"
#include <iostream>

static bool HB_initialized = false;
static bool HS_initialized = false;

#if defined HiggsBounds


void HB_init() {
  int nH0=3;
  int nHp=1;
  int hbflag=3;	
  
  if (HB_initialized) return;
  
  printf("\nInitializing HiggsBounds... ");
   
   // Third argument is HB analysis setting: 1='onlyL', 2='onlyH' 3='LandH'
    initialize_higgsbounds_int_(&nH0, &nHp, &hbflag);
//  initialize_higgsbounds_(&nH0,&nHp,whichexpt);

  printf("Please cite the relevant publications when using Higgs mass limits.\n");

  HB_initialized = true;

}


void HB_set_input(THDM model) {

	HB_set_input_effC(model);

}


void HB_set_input_effC(THDM model) {

	if (!HB_initialized) {
		cout << "WARNING: HiggsBounds must be initialized with HB_init() before usage" << endl;
		return;
	}

    bool debug =false;

    double Mh[3];
    double GammaTotal[3]; 
    double g2hjss_s[3];
    double g2hjss_p[3];
    double g2hjcc_s[3];
    double g2hjcc_p[3];
    double g2hjbb_s[3];
    double g2hjbb_p[3]; 
    double g2hjtt_s[3];
    double g2hjtt_p[3];
    double g2hjmumu_s[3];
    double g2hjmumu_p[3];
    double g2hjtautau_s[3];
    double g2hjtautau_p[3]; 
    double g2hjWW[3];
    double g2hjZZ[3];
    double g2hjZga[3];
    double g2hjgaga[3];
    double g2hjgg[3];
    double g2hjggZ[3]; 
    double g2hjhiZ[3][3]; 
    double BR_hjinvisible[3];
    double BR_hjhihi[3][3];
    
    double MHp[1];
	double MHplusGammaTot[1];
	double CS_lep_HpjHmi_ratio[1];
	double BR_tWpb[1];
	double BR_tHpjb[1];
	double BR_Hpjcs[1];
	double BR_Hpjcb[1];
	double BR_Hptaunu[1];
    
    double sba, tanb, lam6, lam7, m122;
    model.get_param_phys(Mh[0], Mh[1], Mh[2], MHp[0], sba, lam6, lam7, m122, tanb); 
    
#if defined debug
    printf("Masses for HB/HS: %16.8E %16.8E %16.8E\n", Mh[0], Mh[1], Mh[2]);
#endif

    THDM sm_like;
    DecayTable table(model), sm_table(sm_like);
 
    SM sm = model.get_SM();

  	double g=sm.get_g();
  	double costw=sm.get_costw();
  	double mt = sm.get_umass_pole(3);
 
    complex <double> c,cs,cp,c_sm,cs_sm,cp_sm;
    
    for (int h=1;h<=3;h++) {
      double mh = Mh[h-1];
//      sm_like.set_param_phys(Mh[h-1], Mh[h-1]*10., Mh[h-1]*10., Mh[h-1]*10., 1.0, 0., 0., 0., 1.0);
      sm_like.set_param_sm(Mh[h-1]);
      sm_like.set_yukawas_type(1);
      
      table.set_model(model);
      sm_table.set_model(sm_like);
    
      model.get_coupling_hdd(h,2,2,cs,cp);
	  sm_like.get_coupling_hdd(1,2,2,cs_sm,cp_sm);
      g2hjss_s[h-1] = pow(abs(cs/cs_sm),2);
      g2hjss_p[h-1] = pow(abs(cp/cs_sm),2);
	  if (debug) printf("%2d %5s %16.8E %16.8E\n", h, "ss", g2hjss_s[h-1], g2hjss_p[h-1]);

      model.get_coupling_hdd(h,3,3,cs,cp);
	  sm_like.get_coupling_hdd(1,3,3,cs_sm,cp_sm);
      g2hjbb_s[h-1] = pow(abs(cs/cs_sm),2);
      g2hjbb_p[h-1] = pow(abs(cp/cs_sm),2);
      if (debug) printf("%2d %5s %16.8E %16.8E\n", h, "bb", g2hjbb_s[h-1], g2hjbb_p[h-1]);

      model.get_coupling_huu(h,2,2,cs,cp);
	  sm_like.get_coupling_huu(1,2,2,cs_sm,cp_sm);
      g2hjcc_s[h-1] = pow(abs(cs/cs_sm),2);
      g2hjcc_p[h-1] = pow(abs(cp/cs_sm),2);
      if (debug) printf("%2d %5s %16.8E %16.8E\n", h, "cc", g2hjcc_s[h-1], g2hjcc_p[h-1]);

      model.get_coupling_huu(h,3,3,cs,cp);
	  sm_like.get_coupling_huu(1,3,3,cs_sm,cp_sm);
      g2hjtt_s[h-1] = pow(abs(cs/cs_sm),2);
      g2hjtt_p[h-1] = pow(abs(cp/cs_sm),2);
      if (debug) printf("%2d %5s %16.8E %16.8E\n", h, "tt", g2hjtt_s[h-1], g2hjtt_p[h-1]);

      model.get_coupling_hll(h,2,2,cs,cp);
	  sm_like.get_coupling_hll(1,2,2,cs_sm,cp_sm);
      g2hjmumu_s[h-1] = pow(abs(cs/cs_sm),2);
      g2hjmumu_p[h-1] = pow(abs(cp/cs_sm),2);
      if (debug) printf("%2d %5s %16.8E %16.8E\n", h, "mumu", g2hjmumu_s[h-1], g2hjmumu_p[h-1]);

      model.get_coupling_hll(h,3,3,cs,cp);
	  sm_like.get_coupling_hll(1,3,3,cs_sm,cp_sm);
      g2hjtautau_s[h-1] = pow(abs(cs/cs_sm),2);
      g2hjtautau_p[h-1] = pow(abs(cp/cs_sm),2);
      if (debug) printf("%2d %5s %16.8E %16.8E\n", h, "tata", g2hjtautau_s[h-1], g2hjtautau_p[h-1]);

      model.get_coupling_vvh(2, 2, h, c);
      sm_like.get_coupling_vvh(2, 2, 1, c_sm);
      g2hjZZ[h-1] = pow(abs(c/c_sm),2);
      if (debug) printf("%2d %5s %16.8E\n", h, "ZZ", g2hjZZ[h-1]);

      model.get_coupling_vvh(3, 3, h, c);
      sm_like.get_coupling_vvh(3, 3, 1, c_sm);
      g2hjWW[h-1] = pow(abs(c/c_sm),2);
      if (debug) printf("%2d %5s %16.8E\n", h, "WW", g2hjWW[h-1]);
    
      double hgaga = table.get_gamma_hgaga(h);
      double hgaga_sm = sm_table.get_gamma_hgaga(1);
      g2hjgaga[h-1] = hgaga/hgaga_sm;
      if (debug) printf("%2d %5s %16.8E\n", h, "gaga",g2hjgaga[h-1]);

      double hZga = table.get_gamma_hZga(h);
      double hZga_sm = sm_table.get_gamma_hZga(1);
      g2hjZga[h-1] = hZga/hZga_sm;
      if (debug) printf("%2d %5s %16.8E\n", h, "Zga", g2hjZga[h-1]);
    
      double hgg = table.get_gamma_hgg(h);
      double hgg_sm = sm_table.get_gamma_hgg(1);
      g2hjgg[h-1] = hgg/hgg_sm;
      g2hjggZ[h-1] = 0.;
      if (debug) printf("%2d %5s %16.8E\n", h, "gg", g2hjgg[h-1]);    

      if ((h<=2)&&(mh>=90.)) {    
		  GammaTotal[h-1] =  smgamma_h_(&mh);
		  GammaTotal[h-1] += smgamma_h_(&mh)*(g2hjWW[h-1] - 1.)*smbr_hww_(&mh);
		  GammaTotal[h-1] += smgamma_h_(&mh)*(g2hjZZ[h-1] - 1.)*smbr_hzz_(&mh);
		  GammaTotal[h-1] += smgamma_h_(&mh)*(g2hjgg[h-1] - 1.)*smbr_hgg_(&mh);
		  GammaTotal[h-1] += smgamma_h_(&mh)*(g2hjtt_s[h-1] + g2hjtt_p[h-1]/(1.-2.*pow(mt,2)/pow(mh,2)) - 1.)*smbr_htoptop_(&mh);
		  GammaTotal[h-1] += smgamma_h_(&mh)*(g2hjbb_s[h-1] + g2hjbb_p[h-1] - 1.)*smbr_hbb_(&mh);
		  GammaTotal[h-1] += smgamma_h_(&mh)*(g2hjtautau_s[h-1] + g2hjtautau_p[h-1] - 1.)*smbr_htautau_(&mh);
		  GammaTotal[h-1] += smgamma_h_(&mh)*(g2hjmumu_s[h-1] + g2hjmumu_p[h-1] - 1.)*smbr_hmumu_(&mh);
		  GammaTotal[h-1] += smgamma_h_(&mh)*(g2hjss_s[h-1] + g2hjss_p[h-1] - 1.)*smbr_hss_(&mh);
		  GammaTotal[h-1] += smgamma_h_(&mh)*(g2hjcc_s[h-1] + g2hjcc_p[h-1] - 1.)*smbr_hcc_(&mh);
		  GammaTotal[h-1] += smgamma_h_(&mh)*(g2hjZga[h-1]-1)*smbr_hzgam_(&mh);
		  GammaTotal[h-1] += smgamma_h_(&mh)*(g2hjgaga[h-1]-1)*smbr_hgamgam_(&mh);
	  
		  for (int i=1; i<=4; i++) {
			  GammaTotal[h-1] += table.get_gamma_hhh(h,i,i);
		  }

		  for (int i=1; i<=4; i++) {
			for (int j=1; j<=3; j++) {
			  GammaTotal[h-1] += table.get_gamma_hvh(h,j,i);
			}
		  }
	  
	  } else {
	  
	    GammaTotal[h-1] = table.get_gammatot_h(h);
	  }

      
      if (debug) printf("gtot %16.8E %16.8E %16.8E %16.8E\n",  GammaTotal[h-1], table.get_gammatot_h(h), HB_get_gammah(Mh[h-1]), sm_table.get_gammatot_h(1));


    }
    
    
    

//      GammaTotal[h-1] = GammaTotal[h-1]*HB_get_gammah(Mh[h-1])/sm_table.get_gammatot_h(1);
//      	printf("gtot %16.8E %16.8E %16.8E\n",  GammaTotal[h-1], sm_table.get_gammatot_h(1), HB_get_gammah(Mh[h-1]));
  	 
  	 
  	 for (int j=1;j<=3;j++) {  	 
      for (int i=1;i<=3;i++) {
       BR_hjhihi[i-1][j-1]=table.get_gamma_hhh(j,i,i)/GammaTotal[j-1];
       model.get_coupling_vhh(2,j,i,c);
       g2hjhiZ[i-1][j-1]=pow(abs(c)/(g/2./costw),2);
       if (debug) printf("%2d %2d hihjZ %16.8E\n", j, i, g2hjhiZ[i-1][j-1]);
       if (debug) printf("%2d %2d hj->hihi %16.8E\n", j, i, BR_hjhihi[i-1][j-1]);
      }
     }

    higgsbounds_neutral_input_effc_(Mh,GammaTotal,
    g2hjss_s, g2hjss_p,
    g2hjcc_s, g2hjcc_p,
    g2hjbb_s, g2hjbb_p, 
    g2hjtt_s, g2hjtt_p,
    g2hjmumu_s,g2hjmumu_p,
    g2hjtautau_s, g2hjtautau_p, 
    g2hjWW, g2hjZZ, g2hjZga, g2hjgaga,
    g2hjgg, g2hjggZ,g2hjhiZ, 
    BR_hjinvisible,
    BR_hjhihi);
    
    CS_lep_HpjHmi_ratio[0] = 1.;
	BR_tWpb[0] = sm.get_gamma_top()/table.get_gammatot_top();
	BR_tHpjb[0]=table.get_gamma_uhd(3,4,3)/table.get_gammatot_top();

	BR_Hpjcs[0] = table.get_gamma_hdu(4,2,2)/table.get_gammatot_h(4);
	BR_Hpjcb[0] = table.get_gamma_hdu(4,3,2)/table.get_gammatot_h(4);
	BR_Hptaunu[0] = table.get_gamma_hln(4,3,3)/table.get_gammatot_h(4);

	higgsbounds_charged_input_(MHp,
		MHplusGammaTot,
		CS_lep_HpjHmi_ratio,
		BR_tWpb,
		BR_tHpjb,
		BR_Hpjcs,
		BR_Hpjcb,
		BR_Hptaunu);

}


double HB_get_gammah(double m) {
 
   return smgamma_h_(&m);
}


void HB_run_full(int hbres[], int hbchan[], double hbobs[], int hbcomb[]) {

//    printf("Running HB full\n");
	if (!HB_initialized) {
		cout << "WARNING: HiggsBounds must be initialized with HB_init() before usage" << endl;
		return;
	}

	run_higgsbounds_full_(hbres,hbchan,hbobs,hbcomb);
}

void HS_init() {

 int nH0 = 3;
 int nHp = 1;
 
 printf("\nInitializing HiggsSignals... ");

 initialize_higgssignals_latestresults_(&nH0,&nHp);
 HS_initialized = true;  
}

void HS_finish() {

  finish_higgssignals_();
}

void HB_finish() {

  finish_higgsbounds_();
}


void HS_run(double *csqmu, double *csqmh, double *csqtot, int *nobs, double *pval) {
	if (!HS_initialized) {
		cout << "WARNING: HiggsSignals must be initialized with HS_init() before usage" << endl;
		return;
	}

	int mode = 1;

  run_higgssignals_(&mode, csqmu, csqmh, csqtot, nobs, pval);
}


void HS_set_pdf(int pdf) {

  setup_pdf_(&pdf);

}

void HS_set_assignment_range(double range) {
  setup_assignmentrange_(&range);
}

void HS_setup_assignment_range_massobservables(double range) {
  setup_assignmentrange_massobservables_(&range);
}


void HS_set_rate_uncertainties(double dCS[], double dBR[]) {
  setup_rate_uncertainties_(dCS, dBR);

}

void HS_set_mass_uncertainties(double dMh[]) {
	higgssignals_neutral_input_massuncertainty_(dMh);
}

void HS_get_Rvalues(int i, int collider, double *R_H_WW, double *R_H_ZZ, double *R_H_gaga, double *R_H_tautau, double *R_H_bb, double *R_VH_bb) {
  get_rvalues_(&i, &collider, R_H_WW, R_H_ZZ, R_H_gaga, R_H_tautau, R_H_bb, R_VH_bb);

}

void HS_set_output_level(int level) {

  setup_output_level_(&level);

}

#endif

