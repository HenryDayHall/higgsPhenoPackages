/*******************
 * Quick demo program
 * ****************/
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

    // Set parameters of the 2HDM in the 'physical' basis
    // Some parameters stay the same
    double mh       = 125.;
    double mH       = (double)atof(argv[1]); // 100 1000
    double mA       = (double)atof(argv[2]);  //100 1000
    double mC       = (double)atof(argv[3]); //100 1000
    double sba      = (double)atof(argv[4]);  //0.9 1
    double lambda_6 = 0.;  // 
    double lambda_7 = 0.;
    double tb       = (double)atof(argv[5]);  //0 25
    int yt_in = (int)atoi(argv[6]); //Type



    double m12_2    = mA*mA*sin(atan(tb))*cos(atan(tb));

    printf("Using parameters\n");
    printf("mh:       %8.4f\n", mh);
    printf("mH:       %8.4f\n", mH);
    printf("mA:       %8.4f\n", mA);
    printf("tan(b):   %8.4f\n", tb);
    printf("mC:       %8.4f\n", mC);
    printf("sin(b-a): %8.4f\n", sba);
    printf("lambda_6:       %8.4f\n", lambda_6);
    printf("lambda_7:       %8.4f\n", lambda_7);
    printf("m12_2:       %8.4f\n", m12_2);


    

    bool pset = model.set_param_phys(mh,mH,mA,mC,sba,lambda_6,lambda_7,m12_2,tb);
    if (!pset) {
      cerr << "The specified parameters are not valid\n";
      return -1;
    }

    // Set Yukawa couplings to model type
    model.set_yukawas_type(yt_in);
    string sep = " ";
    


    // Prepare to calculate observables
    Constraints constr(model);
    bool BitAllowedStability      = constr.check_stability();
    bool BitAllowedUnitarity      = constr.check_unitarity();
    bool BitAllowedPerturbativity = constr.check_perturbativity();
    double S,T,U,V,W,X;   

    constr.oblique_param(mh_ref,S,T,U,V,W,X);

    
  
  DecayTable table(model);
       double BrhWW = table.get_gamma_hvv(1,3)/table.get_gammatot_h(1);
  double BrhZZ = table.get_gamma_hvv(1,2)/table.get_gammatot_h(1);
  double Brhtautau = table.get_gamma_hll(1,3,3)/table.get_gammatot_h(1);

  std::cout << "Stability = "<< std::fixed << BitAllowedStability << "\n"; //10
  std::cout << "Unitarity = " << std::fixed << BitAllowedUnitarity << "\n";  //11
  std::cout << "Perturbativity = " << std::fixed << BitAllowedPerturbativity << "\n";  //12
   //############  h decays #####################################
  std::cout << "BrhWW = "<< std::fixed << BrhWW << "\n";
  std::cout << "BrhZZ = " << std::fixed << BrhZZ << "\n";
}
