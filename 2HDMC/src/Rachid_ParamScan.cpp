#include "THDM.h"
#include "SM.h"
#include "HBHS.h"
#include "Constraints.h"
#include "DecayTable.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>

using namespace std;

int main(int argc, char* argv[]) {

  if ( argc != 15)
  {
	 printf("ParameterScan_Physic_MultiDim usage:\n");
	 printf("ParameterScan_Physic_MultiDim <tag> <BitWriteLHA> <mh> <mH>  <mA> <tanb> <mC> <sinba> <l6> <l7> <m12_2> <yukawa-type>");
  }

  double mh_in     = (double)atof(argv[1]);
  double mH_in     = (double)atof(argv[2]);
  double mA_in    = (double)atof(argv[3]);
  double tanb_in   = (double)atof(argv[4]);
  double mC_in     = (double)atof(argv[5]);
  double sba_in     = (double)atof(argv[6]);
  double l6_in     = (double)atof(argv[7]);
  double l7_in     = (double)atof(argv[8]);
  double m12_2_in     = (double)atof(argv[9]);
  int yt_in        = (int)atoi(argv[10]);
  std::string file_tag  = argv[11];

  printf("Inside ParameterScan_Hybrid_MultiDim\n");
  printf("mh:       %8.4f\n", mh_in);
  printf("mH:       %8.4f\n", mH_in);
  printf("mA:       %8.4f\n", mA_in);
  printf("tan(b):   %8.4f\n", tanb_in);
  printf("mC:       %8.4f\n", mC_in);
  printf("sin(b-a): %8.4f\n", sba_in);
  printf("l6:       %8.4f\n", l6_in);
  printf("l7:       %8.4f\n", l7_in);
  printf("m12_2:       %8.4f\n", m12_2_in);

  // Reference SM Higgs mass for EW precision observables
  double mh_ref = 125.;
  
  // Create SM and set parameters
  SM sm;
  sm.set_qmass_pole(6, 172.5);		
  sm.set_qmass_pole(5, 4.75);		
  sm.set_qmass_pole(4, 1.42);	
  sm.set_lmass_pole(3, 1.77684);	
  sm.set_alpha(1./127.934);
  sm.set_alpha0(1./137.0359997);
  sm.set_alpha_s(0.119);
  sm.set_MZ(91.15349);
  sm.set_MW(80.36951);
  sm.set_gamma_Z(2.49581);
  sm.set_gamma_W(2.08856);
  sm.set_GF(1.16637E-5);

  // Create 2HDM and set SM parameters
  THDM model;
  model.set_SM(sm);

  // Parameter set validation check
  bool pset = model.set_param_phys(mh_in,mH_in,mA_in,mC_in,sba_in,l6_in,l7_in,m12_2_in,tanb_in);
  
  if (!pset) {
    cerr << "The specified parameters are not valid\n";
    return -1;
  }

  HB_init();
  HS_init();

  // Set Yukawa couplings
  model.set_yukawas_type(yt_in);


  // Prepare to calculate observables
  Constraints constr(model);
  bool BitAllowedStability      = constr.check_stability();
  bool BitAllowedUnitarity      = constr.check_unitarity();
  bool BitAllowedPerturbativity = constr.check_perturbativity();
  double S,T,U,V,W,X;   

  constr.oblique_param(mh_ref,S,T,U,V,W,X);

  
  //////////////////////////////////
  // - HiggsBounds/HiggsSignals - //
  //////////////////////////////////

// See HiggsSignals manual for more information
  int mass_pdf = 2;
  HS_set_pdf(mass_pdf);
  HS_setup_assignment_range_massobservables(2.);
  HS_set_output_level(0);

// Share couplings of 2HDM model with HiggsBounds/HiggsSignals
  HB_set_input_effC(model);
  
  int hbres[6];
  int hbchan[6];
  double hbobs[6];
  int hbcomb[6];  

// Run HiggsBounds 'full', i.e. with each Higgs result separately  
  HB_run_full(hbres, hbchan, hbobs, hbcomb);
  printf("\nHiggsBounds results (full):\n");
  printf("  Higgs  res  chan       ratio        ncomb\n");
  for (int i=1;i<=4;i++) {
    printf("%5d %5d %6d %16.8E %5d   %s\n", i, hbres[i],hbchan[i],hbobs[i],hbcomb[i],hbobs[i]<1 ? "Allowed" : "Excluded");
  }
  printf("------------------------------------------------------------\n");
  printf("  TOT %5d %6d %16.8E %5d   %s\n", hbres[0],hbchan[0],hbobs[0],hbcomb[0],hbobs[0]<1 ? "ALLOWED" : "EXCLUDED");

  double tot_hbobs = hbobs[0];
  double csqmu;
  double csqmh;
  double csqtot;
  int nobs;
  double pval;
  
  double dMh[3]={0., 0., 0.,};
  HS_set_mass_uncertainties(dMh);
  HS_run(&csqmu, &csqmh, &csqtot, &nobs, &pval);
  printf("\nHiggsSignals results:\n");
  printf(" Chi^2 from rates: %16.8E\n", csqmu);
  printf("  Chi^2 from mass: %16.8E\n", csqmh);
  printf("      Total chi^2: %16.8E\n", csqtot);
  printf("    # observables: %16d\n\n", nobs);


  
  DecayTable table(model);
       double BrhWW = table.get_gamma_hvv(1,3)/table.get_gammatot_h(1);
  double BrhZZ = table.get_gamma_hvv(1,2)/table.get_gammatot_h(1);
  double Brhtautau = table.get_gamma_hll(1,3,3)/table.get_gammatot_h(1);
  double Brhss = table.get_gamma_hdd(1,2,2)/table.get_gammatot_h(1);
  double Brhbb = table.get_gamma_hdd(1,3,3)/table.get_gammatot_h(1);
  double Brhuu = table.get_gamma_huu(1,1,1)/table.get_gammatot_h(1);
  double Brhcc = table.get_gamma_huu(1,2,2)/table.get_gammatot_h(1);
  double Brhgg = table.get_gamma_hgg(1)/table.get_gammatot_h(1);
  double Brhgaga = table.get_gamma_hgaga(1)/table.get_gammatot_h(1);

// H decays #####################################################

  double BrHWW = table.get_gamma_hvv(2,3)/table.get_gammatot_h(2);
  double BrHZZ = table.get_gamma_hvv(2,2)/table.get_gammatot_h(2);
  double BrHtautau = table.get_gamma_hll(2,3,3)/table.get_gammatot_h(2);
  double BrHss = table.get_gamma_hdd(2,2,2)/table.get_gammatot_h(2);
  double BrHbb = table.get_gamma_hdd(2,3,3)/table.get_gammatot_h(2);
  double BrHuu = table.get_gamma_huu(2,1,1)/table.get_gammatot_h(2);
  double BrHcc = table.get_gamma_huu(2,2,2)/table.get_gammatot_h(2);
  double BrHtt = table.get_gamma_huu(2,3,3)/table.get_gammatot_h(2);
  double BrHHpHm = table.get_gamma_hhh(2,4,4)/table.get_gammatot_h(2);
  double BrHhh = table.get_gamma_hhh(2,1,1)/table.get_gammatot_h(2);
  double BrHAA = table.get_gamma_hhh(2,3,3)/table.get_gammatot_h(2);
  double BrHgg = table.get_gamma_hgg(2)/table.get_gammatot_h(2);
  double BrHgaga = table.get_gamma_hgaga(2)/table.get_gammatot_h(2);
  double BrHZA = table.get_gamma_hvh(2,2,3)/table.get_gammatot_h(2);

// A decays #####################################################

  double BrAtautau = table.get_gamma_hll(3,3,3)/table.get_gammatot_h(3);
  double BrAss = table.get_gamma_hdd(3,2,2)/table.get_gammatot_h(3);
  double BrAbb = table.get_gamma_hdd(3,3,3)/table.get_gammatot_h(3);
  double BrAuu = table.get_gamma_huu(3,1,1)/table.get_gammatot_h(3);
  double BrAcc = table.get_gamma_huu(3,2,2)/table.get_gammatot_h(3);
  double BrAtt = table.get_gamma_huu(3,3,3)/table.get_gammatot_h(3);
  double BrAgg = table.get_gamma_hgg(3)/table.get_gammatot_h(3);
  double BrAgaga = table.get_gamma_hgaga(3)/table.get_gammatot_h(3);
  double BrAZh = table.get_gamma_hvh(3,2,1)/table.get_gammatot_h(3);
  double BrAZH = table.get_gamma_hvh(3,2,2)/table.get_gammatot_h(3);


// Hp decays #####################################################
  double BrHpmunu = table.get_gamma_hln(4,2,2)/table.get_gammatot_h(4); 
  double BrHptanu = table.get_gamma_hln(4,3,3)/table.get_gammatot_h(4); 
  double BrHptb = table.get_gamma_hdu(4,3,3)/table.get_gammatot_h(4);
  double BrHpcs = table.get_gamma_hdu(4,2,2)/table.get_gammatot_h(4);
  double BrHpcb = table.get_gamma_hdu(4,3,2)/table.get_gammatot_h(4);
  double BrHpWh = table.get_gamma_hvh(4,3,1)/table.get_gammatot_h(4);
  double BrHpWH = table.get_gamma_hvh(4,3,2)/table.get_gammatot_h(4);
  double BrHpWA = table.get_gamma_hvh(4,3,3)/table.get_gammatot_h(4);
//#######################################################################
       
  std::string filename_param_chisq = file_tag;

  std::ofstream file_param_chisq;
  file_param_chisq.open(filename_param_chisq.c_str(), std::ios_base::app);

  std::stringstream line;
   line << std::setprecision(6) << std::fixed << mh_in << " ";    // 1
   line << std::setprecision(6) << std::fixed << mH_in << " ";    //2
   line << std::setprecision(6) << std::fixed << mA_in << " ";    //3
   line << std::setprecision(6) << std::fixed << mC_in << " ";    //4
   line << std::setprecision(6) << std::fixed << sba_in << " ";   //5
   line << std::setprecision(6) << std::fixed << tanb_in << " ";  //6
   line << std::setprecision(6) << std::fixed << m12_2_in << " "; //7
   line << std::setprecision(6) << std::fixed << tot_hbobs << " ";//8
   line << std::setprecision(6) << std::fixed << csqtot << " ";   //9  
   line << std::setprecision(6) << std::fixed << pval << " ";
   line << std::setprecision(1) << std::fixed << BitAllowedStability << " "; //10
   line << std::setprecision(1) << std::fixed << BitAllowedUnitarity << " ";  //11
   line << std::setprecision(1) << std::fixed << BitAllowedPerturbativity << " ";  //12
   //############  h decays #####################################
  line << std::setprecision(6) << std::fixed << BrhWW << " ";
  line << std::setprecision(6) << std::fixed << BrhZZ << " ";
  line << std::setprecision(6) << std::fixed << Brhtautau << " ";
  line << std::setprecision(6) << std::fixed << Brhss << " ";
  line << std::setprecision(6) << std::fixed << Brhbb << " ";
  line << std::setprecision(6) << std::fixed << Brhuu << " ";
  line << std::setprecision(6) << std::fixed << Brhcc << " ";
  line << std::setprecision(6) << std::fixed << Brhgg << " ";
  line << std::setprecision(6) << std::fixed << Brhgaga << " ";

//########### H decays  #########################################

  line << std::setprecision(6) << std::fixed << BrHWW << " ";
  line << std::setprecision(6) << std::fixed << BrHZZ << " ";
  line << std::setprecision(6) << std::fixed << BrHtautau << " ";
  line << std::setprecision(6) << std::fixed << BrHss << " ";
  line << std::setprecision(6) << std::fixed << BrHbb << " ";
  line << std::setprecision(6) << std::fixed << BrHuu << " ";
  line << std::setprecision(6) << std::fixed << BrHcc << " ";
  line << std::setprecision(6) << std::fixed << BrHtt << " ";
  line << std::setprecision(6) << std::fixed << BrHHpHm << " ";
  line << std::setprecision(6) << std::fixed << BrHhh << " ";
  line << std::setprecision(6) << std::fixed << BrHAA << " ";
  line << std::setprecision(6) << std::fixed << BrHgg << " ";
  line << std::setprecision(6) << std::fixed << BrHgaga << " ";
  line << std::setprecision(6) << std::fixed << BrHZA << " ";

//##########  A decays  ###########################################

  line << std::setprecision(6) << std::fixed << BrAtautau << " ";
  line << std::setprecision(6) << std::fixed << BrAss << " ";
  line << std::setprecision(6) << std::fixed << BrAbb << " ";
  line << std::setprecision(6) << std::fixed << BrAuu << " ";
  line << std::setprecision(6) << std::fixed << BrAcc << " ";
  line << std::setprecision(6) << std::fixed << BrAtt << " ";
  line << std::setprecision(6) << std::fixed << BrAgg << " ";
  line << std::setprecision(6) << std::fixed << BrAgaga << " ";
  line << std::setprecision(6) << std::fixed << BrAZh << " ";
  line << std::setprecision(6) << std::fixed << BrAZH << " ";
//#########  H+ decays ################################################

  line << std::setprecision(6) << std::fixed << BrHptb << " ";
  line << std::setprecision(6) << std::fixed << BrHpcs << " ";
  line << std::setprecision(6) << std::fixed << BrHpcb << " ";
  line << std::setprecision(6) << std::fixed << BrHpWh << " ";
  line << std::setprecision(6) << std::fixed << BrHpWA << " ";
  line << std::setprecision(6) << std::fixed << BrHpWH << " ";
  line << std::setprecision(6) << std::fixed << BrHpmunu << " ";
  line << std::setprecision(6) << std::fixed << BrHptanu << " ";
  line << std::setprecision(6) << std::fixed << S << " ";
  line << std::setprecision(6) << std::fixed << T << " ";
  line << std::setprecision(6) << std::fixed << U << "\n";
  file_param_chisq << line.rdbuf();
  
  
  HB_finish();
  HS_finish();
}
