// TODO: This documentation is outdated, from this line:
//     input_file >> timestep >> temp1 >> rad1 >> cpTemp >> iLeak >> diLeak
// I deduce that the profile file must contain:
//     duration, temperature, fluence, cooling pipe temperature, leakage current, std leakage current
// which seems to correspond to what this script does:
// https://gitlab.cern.ch/dbrzhech/PixelMonitoring/-/blob/b795aaee606c38a91ed9158f1a0addc43f511a7d/getPLC_Tpipe.py
// Which, unlike the names says, prepares a file with this input.
// However it does so for FPix, and does not run anymore.
//
/*
Compile Code with: g++ data_processing.cpp -I /path/to/boost/boost_1_62_0 -Wall -std=c++11 -o output `root-config --cflags` `root-config --libs`

Requirements: Root (for plotting) and Boost (for Date handling) have to be installed

Call program with 1 Arguments:  1. name of the file with the temperature/radiation profile (you may want to uncomment the gA,Y,C to be read in again as additional parameters, around line 686)

Radiation and temperature profile have to fulfil the following scheme: duration (in seconds!!!!!!!!!) / temperature in K / irradiation in neq/cm2s (ALL have to be integer !)

Attention:  timesteps in the profile file HAVE to be devideable by global timestep ->Integer result required
            profile files may not have additional (even/especially empty) lines after the last line

Attention 2: the value for the global_layer_conversion has to be changed in order to fit to another layer (parameter transforms luminosity in neq) eg: B-Layer: 2.5e12, IBL 6e12

Attention 3: to change from one detector type to another, change mainly 3 things: thickness and global_layer_conversion (global values) and Ndonor
*/


#include <iostream>			// Basic IO
#include <fstream>			// Read and write files
#include <cstdlib>			// convert arguments of function call into int, float,...
#include <string>
#include <vector>
#include <math.h>           		// for fabs()
#include <sstream>         		// to get string into stream
//#include <TH2F.h>           		// root stuff for file reading
#include <TFile.h>          		// more root stuff
#include <TH2D.h>           		// root stuff for file reading
#include <TCanvas.h>
#include <TROOT.h>
#include <TGraphErrors.h>
#include <TVector.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TF1.h>
#include <TMath.h>
#include <TLine.h>
#include <TH1.h>
#include <cmath>
#include "TDatime.h"
#include <time.h>
#include <ctime>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

using namespace std;
// using namespace boost::gregorian;
// using namespace boost::posix_time;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool use_CMS_paper_alpha0 = false; //if false: temperature averaging for leakage current calculation a la Moll (theta function for time scale change), if true: use CMS approach of temperature average (attention: constants can only be changed directly in the code in this case)
// begin date: 2017-05-23 14:32:22.210863
bool debug=false;                                    // additional debugging console outputs
bool plot_pdf=false;                                  // plots are by default saved as root file, but can also be directly exported as pdf
bool overwrite_files=true;                           // false = append internal clock identifier to file name, false = overwrite pdfs, append root files // files are by default created with a unique name depending on the internal clock and you can switch it of here such that pdf files are overwritten and the root file is extended (old plots will stay in the root file as well)
double timestep=1;                                   // step size is 1 second, do not change this!
double donorremovalfraction=0.99;                    // fraction of donors which can be removed by donor removal
const double userTrefC=0.;
double userTref=273.15+userTrefC;                          // set a reference temperature for the volume corrected leakage current plot (it will only effect this one plot!) Now: implemented!
double bandGap=1.21;                                 // eV used for scaling temperatures

// string output_file_name = "simulation_results";      // set unique file name for each simulation
string startTime = "2017-05-23 14:32:22";
// date d(2017,May,23);                                  // IBL     //set a date for the plots to begin (to be correct it has to be equal to the beginning of your (!) temp/irr profile)
boost::posix_time::ptime t1(boost::posix_time::time_from_string(startTime));
//date d(2011,Feb,11);                               // PIXEL


const double Ndonor_0 = 1.7e12;                      // IBL       // initial donor concentration (right now the code only works for originally n-type bulk material!)
//const double Ndonor_0 = 1.4e12;                    // B-Layer Layer1/2 Disks

// double thickness=200;                                // IBL       // sensor thickness
//double thickness=250;                              // B-Layer Layer1/2 Disks
double thickness=285;                                // CMS sensor thickness
//double global_layer_conversion=0.92e12;            // Layer 1    //conversion factor from luminosity to neq (1fb-1=2.3e12neq/cm2) - is used only for computation of total luminosity, will be wrong if there are multiple conversion factor in the original profiles as for example due to different center of mass energies
//double global_layer_conversion=1.1e12;             // L1 average
//double global_layer_conversion=0.571e12;           // Layer 2
//double global_layer_conversion=0.7e12;             // Layer 2 average
//double global_layer_conversion=0.582e12;           // Disks
double global_layer_conversion=6.262e12;             // IBL
//double global_layer_conversion=2.8e12;//2.929e12;  // B-Layer

int limit_file_size_input = 0;                       // how many lines should be read from the profile, 0 for everything

// TODO: Understand this scaling, why is it "temporary"?
float DoseRateScaling=0.0859955711871/8.9275E-02; //temporary scaling for bpix
// float DoseRateScaling=0.0859955711871/8.9275E-02;                            // this parameter is multiplied to all fluence rates from the input profile

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Do not change things below this line without thinking carefully!

//Definition of Hamburg model depletion voltage functions

TF1 *Naccept_rev_TF1 = new TF1("Naccept_rev_TF1","[0]*[1]*(1-TMath::Exp(-[2]*x))/[2] + [3]*TMath::Exp(-[2]*x)",0,10000000);
TF1 *Naccept_const_TF1 = new TF1("Naccept_const_TF1","[0] * [1]*x",0,10000000);
TF1 *Nneutrals_rev_TF1 = new TF1("Ndonor_rev_TF1","[0]*[1]*(1-TMath::Exp(-[2]*x))/[2] + [3]*TMath::Exp(-[2]*x)",0,10000000);
TF1 *Ndonor_const_TF1 = new TF1("Ndonor_const_TF1","-[0]*(1-TMath::Exp(-[1]*[2]*x))",0,10000000);

TF1 *Naccept_rev_TF1_approx = new TF1("Naccept_rev1_TF1","[0]*[1]*x+[3]*TMath::Exp(-[2]*x)",0,10000000);
TF1 *Ndonor_neutrals_rev_TF1_approx = new TF1("Ndonor_rev1_TF1","[0]*[1]*x+[3]*TMath::Exp(-[2]*x)",0,10000000);

TF1 *Ndonor_TF1 = new TF1("Ndonor_TF1","[0]*[1]/[2] * ([2]*x+TMath::Exp(-[2]*x)-1) +[3]*(1-TMath::Exp(-[2]*x)) ",0,10000000);
TF1 *Ndonor_g1_TF1 = new TF1("Ndonor_TF1","[0]/[1] * ([1]*x+TMath::Exp(-[1]*x)-1) +[2]*(1-TMath::Exp(-[1]*x)) ",0,10000000);

// TODO: Update these paths
TFile *fin = TFile::Open("../PixelMonitoring/FLUKA/fluence_field.root");
TH2F *fluence = (TH2F*)fin->Get("fluence_allpart_6500GeV_phase1");

struct DataElement
{
    int duration;
    float temperature;
    long int dose_rate;
    float coolPipeTemp;
    float iLeak_data;
    float diLeak_data;
};

double deltaT(double flux){
  const double scale=3/1e8;
  //const double scale=3/2.5e7;
  return scale*flux;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//calculate depletion voltage from given effective doping concentration

//some problem with memory !!!!!!

double NtoV_function(double doping_conc)
{
    return 1.6021766208e-13/(2*11.68*8.854187817)*fabs(doping_conc)*thickness*thickness;                // q/2*epsilon*epsilon_0 * Neff*D^2
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double getPhi_eq(double r0, double z0, double phi_aver, double Fluka_aver){
  double phi_eq = fluence->GetBinContent(fluence->GetXaxis()->FindBin(r0),fluence->GetYaxis()->FindBin(z0));
  phi_eq = phi_aver*phi_eq/Fluka_aver;
  return phi_eq;
}

double getFluka_aver(double r0, double z0){
  const double R = 3.2;
  const double dr = 0.1;
  double r = r0-R+dr/2.;
  double phi_eq = 0;
  while(r <= r0+R){
    phi_eq += fluence->GetBinContent(fluence->GetXaxis()->FindBin(r0),fluence->GetYaxis()->FindBin(z0));
    r += dr;
  }
  return phi_eq*dr/(2*R);
}

class Annealing_constants {                                         // Class to store a set of annealing specific constants and compute temperature dependent values
private:
    double gA;
    double gY;
    double gC;
    double ka;
    double ky1;
    double Ea;
    double Ey;
    double cc;

public:
    Annealing_constants(double,double,double,double,double,double,double,double);
    double get_gA(int) const;
    double get_gY(int) const;
    double get_gC(int) const;
    double get_ka(int) const;
    double get_ky1(int) const;
    double get_Ea() const;
    double get_Ey() const;
    double get_cc() const;
};

Annealing_constants::Annealing_constants(double a ,double b ,double c ,double d ,double e ,double f ,double g , double h)
{
    gA=a ;
    gY=b ;
    gC=c ;
    ka=d ;
    ky1=e ;
    Ea=f ;
    Ey=g ;
    cc=h ;
}

double Annealing_constants::get_gA(int T) const
{
    return gA;
}

double Annealing_constants::get_gY(int T) const
{
    return gY;
}

double Annealing_constants::get_gC(int T) const
{
    return gC;
}

double Annealing_constants::get_ka(int T) const
{
    return ka*exp(-Ea/(8.6173303e-5*(double)T));
    //return ka;
}

double Annealing_constants::get_ky1(int T) const
{
    return ky1*exp(-Ey/(8.6173303e-5*(double)T));
    //return ky1;
}

double Annealing_constants::get_Ea() const
{
    return Ea;
}

double Annealing_constants::get_Ey() const
{
    return Ey;
}

double Annealing_constants::get_cc() const
{
    return cc;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// leakage current constants will be stored in this struct

struct leakage_current_consts
{
  double alpha_1;
  double alpha_0_star;
  double beta;
  double k01;
  double E1;
  double E1_star;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// all atributes of a sensor (like irradiation, doping concentration, leakage currents, etc) will be stored in this class

class Sensor                    //Class to handle Sensors
{
private:							//Content of a sensor
	double Nacceptor;
	double Ndonor;
  double Ndonor_const;
  double Nacceptor_reversible;
  double Nneutral_reversible;
  double Nacceptor_stable_constdamage;
  double Ndonor_stable_donorremoval;
  double Nacceptor_stable_reverseannealing;

  double Nbenef_anneal_g1;
  double Nrevers_anneal_g1;
  double Nnadefects_g1;
  double Nconstdamage_g1;

	int Temperature;
  vector<vector<double>> leakage_current_alpha;
  vector<vector<double>> time_history;
  vector<double> G_i;
  vector<double> leakage_current;
  double volume;
  vector<double> alpha_vec;
  vector<double> powerconsumption;
  vector<double> inverse_weighted_temp;


public:
	Sensor(double, double, int, double, double);		    //Constructor
  Sensor();
	double get_Neff() const;
  void set_Nacceptor(double);
  void set_Ndonor(double);
  double get_Nacceptor() const;
  double get_Ndonor() const;

  void set_temp(int);
  double get_temp() const;
  void irradiate(leakage_current_consts,Annealing_constants, long int, float, long double);
  vector<double> get_G_i() const;
  vector<double> get_leakage_current() const;
  vector<double> get_alpha_vec() const;
  vector<double> get_powerconsumption() const;

  double get_Ndonor_const() const;
  double get_Nbenef_anneal_g1() const;
  double get_Nrevers_anneal_g1() const;
  double get_Nnadefects_g1() const;
  double get_Nconstdamage_g1() const;
  // double get_Ndonor_const() const;
  double get_Ndonor_stable_donorremoval() const;
};

Sensor::Sensor( double a, double b, int c, double d, double e)	//Constructordefinition
{
  Nacceptor=a;
  Ndonor=b;
  Temperature=c;
  Nacceptor_reversible=d;
  Nneutral_reversible=e;
}

Sensor::Sensor()	//Constructordefinition
{
  Nacceptor=0;
  Ndonor=Ndonor_0;
  Temperature=0;
  Nneutral_reversible=0;
  Nacceptor_reversible=0;
  Nnadefects_g1=0;
  // Nneutral_reversible=0;
  // volume= 0.135;
  volume=0.0285*6.48*1.62;
  Nacceptor_stable_constdamage=0;
  Ndonor_stable_donorremoval=donorremovalfraction*fabs(Nacceptor-Ndonor);
  Ndonor_const=Ndonor-Ndonor_stable_donorremoval;
  Nacceptor_stable_reverseannealing=0;
}

double Sensor::get_Neff() const
{
	return Nacceptor-Ndonor;
}

double Sensor::get_temp() const
{
  return Temperature;
}

void Sensor::set_temp(int value)
{
  Temperature=value;
  return;
}

double Sensor::get_Nacceptor() const
{
  return Nacceptor;
}

double Sensor::get_Ndonor() const
{
  return Ndonor;
}

double Sensor::get_Ndonor_const() const
{
  return Ndonor_const;
}

void Sensor::set_Nacceptor(double new_value)
{
  Nacceptor=new_value;
  return;
}

void Sensor::set_Ndonor(double new_value)
{
  Ndonor=new_value;
  return;
}

vector<double> Sensor::get_G_i() const
{
  return G_i;
}

vector<double> Sensor::get_leakage_current() const
{
  return leakage_current;
}

vector<double> Sensor::get_alpha_vec() const
{
  return alpha_vec;
}

vector<double> Sensor::get_powerconsumption() const
{
  return powerconsumption;
}

double Sensor::get_Nbenef_anneal_g1() const{
  return Nbenef_anneal_g1;
}
double Sensor::get_Nrevers_anneal_g1() const{
  return Nrevers_anneal_g1;
}
double Sensor::get_Nnadefects_g1() const{
  return Nnadefects_g1;
}

double Sensor::get_Nconstdamage_g1() const{
  return Nconstdamage_g1;
}

double Sensor::get_Ndonor_stable_donorremoval() const{
  return Ndonor_stable_donorremoval;
}

void Sensor::irradiate(leakage_current_consts leconsts, Annealing_constants constants,long int phi, float t, long double totalDose)
{
  //t=t*3600;                                                                     //conversion from hours to seconds
  double a=1e-30;
  //calculating the effective doping concentration
  // debug=1;
  if(debug) cout << t<<" "<< constants.get_gA(Temperature)<<" "<<phi<<" "<<constants.get_ka(Temperature)<<" "<<Nacceptor<<" "<< constants.get_gC(Temperature)<<" "<< constants.get_gY(Temperature)<<" "<< constants.get_ky1(Temperature) <<" "<<Ndonor <<endl;
  
  //
  // phi/= 4.;
  //
  double R0 = 7.81;
  double R = 6.4/2.;
  double dr = 0.1;
  double r = R0 - R + dr/2+30*dr;
  const double z0 = -32.;
  double phi_aver = phi;
  // double Fluka_aver = getFluka_aver(R0, z0);
  double Fluka_aver = 0.07340012983430498;
  double phiR = 0;
  double phi_aver_check = 0;
  // phi = getPhi_eq(r,z0,phi_aver,Fluka_aver);
  // phi *= (0.0859955711871/8.9275E-02); //(why 8.9275E-02 ???)

  Naccept_rev_TF1->SetParameter(0,constants.get_gA(Temperature));
  Naccept_rev_TF1->SetParameter(1,phi);
  Naccept_rev_TF1->SetParameter(2,constants.get_ka(Temperature));
  Naccept_rev_TF1->SetParameter(3,Nacceptor_reversible);

  Naccept_rev_TF1_approx->SetParameter(0,constants.get_gA(Temperature));
  Naccept_rev_TF1_approx->SetParameter(1,phi);
  Naccept_rev_TF1_approx->SetParameter(2,constants.get_ka(Temperature));
  Naccept_rev_TF1_approx->SetParameter(3,Nacceptor_reversible);

  Naccept_const_TF1->SetParameter(0,constants.get_gC(Temperature));
  Naccept_const_TF1->SetParameter(1,phi);

  Nneutrals_rev_TF1->SetParameter(0,constants.get_gY(Temperature));
  Nneutrals_rev_TF1->SetParameter(1,phi);
  Nneutrals_rev_TF1->SetParameter(2,constants.get_ky1(Temperature));
  Nneutrals_rev_TF1->SetParameter(3,Nneutral_reversible);

  Ndonor_neutrals_rev_TF1_approx->SetParameter(0,constants.get_gY(Temperature));
  Ndonor_neutrals_rev_TF1_approx->SetParameter(1,phi);
  Ndonor_neutrals_rev_TF1_approx->SetParameter(2,constants.get_ky1(Temperature));
  Ndonor_neutrals_rev_TF1_approx->SetParameter(3,Nneutral_reversible);

  Ndonor_const_TF1->SetParameter(0,Ndonor_stable_donorremoval);
  Ndonor_const_TF1->SetParameter(1,constants.get_cc());
  Ndonor_const_TF1->SetParameter(2,phi);

  Ndonor_TF1->SetParameter(0,constants.get_gY(Temperature));
  Ndonor_TF1->SetParameter(1,phi);
  Ndonor_TF1->SetParameter(2,constants.get_ky1(Temperature));
  Ndonor_TF1->SetParameter(3,Nneutral_reversible);


  if(constants.get_ka(Temperature)>a && constants.get_ky1(Temperature)>a)
  {
    Nacceptor_reversible             =  Naccept_rev_TF1->Eval(t);
    Nacceptor_stable_constdamage     +=  Naccept_const_TF1->Eval(t);
    Nneutral_reversible              =  Nneutrals_rev_TF1->Eval(t);
    Ndonor_stable_donorremoval       += Ndonor_const_TF1->Eval(t);
    Nacceptor_stable_reverseannealing   +=  Ndonor_TF1->Eval(t);
  }
  else
  {
    cout << "Potential numerical problem due to ultra low temp and therby caused very small ky1 and ka values. Using approach in order to perform calculation. In general, no problem!"<<endl;

    Nacceptor_reversible           =  Naccept_rev_TF1_approx->Eval(t);
    Nacceptor_stable_constdamage   +=  Naccept_const_TF1->Eval(t);
    Nneutral_reversible            =  Ndonor_neutrals_rev_TF1_approx->Eval(t);
    Ndonor_stable_donorremoval     += Ndonor_const_TF1->Eval(t);
    Nacceptor_stable_reverseannealing +=  Ndonor_TF1->Eval(t);
  }

  Nconstdamage_g1   = Nacceptor_stable_constdamage/constants.get_gC(Temperature);
  Nbenef_anneal_g1  = Nacceptor_reversible/constants.get_gA(Temperature);
  Ndonor_g1_TF1->SetParameter(0,phi);
  Ndonor_g1_TF1->SetParameter(1,constants.get_ky1(Temperature));
  Ndonor_g1_TF1->SetParameter(2,Nnadefects_g1);
  Nrevers_anneal_g1 += Ndonor_g1_TF1->Eval(t);
  Nnadefects_g1     = Nneutral_reversible/constants.get_gY(Temperature);

  // cout << "ka = " << constants.get_ka(Temperature) << "; ky1 = " << constants.get_ky1(Temperature) << ": " << NtoV_function(Nacceptor_stable_reverseannealing) << " = " << NtoV_function(Nacceptor_stable_reverseannealing_g) << "+" << NtoV_function(Nacceptor_stable_reverseannealing_no_g) <<  endl;
  Nacceptor =  Nacceptor_reversible + Nacceptor_stable_constdamage + Nacceptor_stable_reverseannealing;
  Ndonor    =  Ndonor_const         + Ndonor_stable_donorremoval;


//calculating the leackage current in the following part

  vector<double> tmp(3);
  vector<double> tmp2(2);
  vector<double> tmp3(2);

  tmp3.at(0)=0;
  tmp3.at(1)=0;

//////////////////////////////

  double G_i_tmp=0;   // G_i_tmp needs to be outside of the if clause to be existend also afterwards when it is pushed back in the G_i vector, can be put back to pos1 after finished with the if clauses

  if(use_CMS_paper_alpha0)
    {
      double CMS_alpha0_c1=-8.9e-17;
      double CMS_alpha0_c2=4.6e-14;
      double CMS_beta=2.9e-18;

      double inverse_weighted_temp_tmp=0;
      inverse_weighted_temp_tmp = t/(double)Temperature;
      inverse_weighted_temp.push_back(inverse_weighted_temp_tmp);

      tmp.at(0)=phi*t*leconsts.alpha_1;
      tmp.at(1)=phi*t;                                                          //needs to be updated at every itaration in the while loop since the averaged temperature needs to be recalculated
      tmp.at(2)=phi*t*CMS_beta;

      leakage_current_alpha.push_back(tmp);

      tmp2.at(0)=t*leconsts.k01*exp(-leconsts.E1/(8.6173303e-5*(double)Temperature));             // time comes in seconds and k01 is seconds-1 so units are fine here
      tmp2.at(1)=t / 60.0;                                                      // t/60 to convert time from seconds to minutes since this is the proposed unit here

      if(time_history.size() > 0.1) tmp3=time_history.at(time_history.size()-1);

      tmp2.at(0)+=tmp3.at(0);
      tmp2.at(1)+=tmp3.at(1);

      time_history.push_back(tmp2);


      int i=0;
      double temperature_average_i=0;

      while(i<leakage_current_alpha.size())
        {
          temperature_average_i=0;

          for(int j=i; j<leakage_current_alpha.size(); j++)                     //calculate the average inverse temperature weighted with the individual time from the moment of irradiation (i) to now (leakage_current_alpha.size())
            {
              temperature_average_i += inverse_weighted_temp.at(j);
            }

          if(i>=1)                                                              //for i=0 the time is not the difference but just the total time
          {
            temperature_average_i /= (double)(time_history.back().at(1)*60.0-time_history.at(i-1).at(1)*60.0);
            G_i_tmp += leakage_current_alpha.at(i).at(0) * exp(-(time_history.back().at(0)-time_history.at(i-1).at(0))) + leakage_current_alpha.at(i).at(1)*(CMS_alpha0_c1+CMS_alpha0_c2*temperature_average_i) - leakage_current_alpha.at(i).at(2) * log(time_history.back().at(1)-time_history.at(i-1).at(1) )  ;
          }
          else
          {
            temperature_average_i /= (double)(time_history.back().at(1)*60.0);
            G_i_tmp += leakage_current_alpha.at(i).at(0) * exp(-time_history.back().at(0)) + leakage_current_alpha.at(i).at(1)*(CMS_alpha0_c1+CMS_alpha0_c2*temperature_average_i) - leakage_current_alpha.at(i).at(2) * log(time_history.back().at(1)) ;
          }

          i++;
        }
    }

  else
  {
      tmp.at(0)=phi*t*leconsts.alpha_1;
      tmp.at(1)=phi*t*leconsts.alpha_0_star;            //temperature dependence of alpha0 is in the theta function of the time calculation - high temp = longer times and vice versa
      tmp.at(2)=phi*t*leconsts.beta;

      leakage_current_alpha.push_back(tmp);

      tmp2.at(0)=t*leconsts.k01*exp(-leconsts.E1/(8.6173303e-5*(double)Temperature));             // time comes in seconds and k01 is seconds-1 so units are fine here
      tmp2.at(1)=t / 60. * exp(-leconsts.E1_star*(1.0/(double)Temperature-1.0/293.15)/(8.6173303e-5));  //T-ref is hardcoded and fixed, it cannot be changed in a trivial way since the other parameters, especially alpha 0 star were evaluated at this reference temperature!!!    // t/60 to convert time from seconds to minutes since this is the proposed unit here

      if(time_history.size() > 0.1) tmp3=time_history.at(time_history.size()-1);

      tmp2.at(0)+=tmp3.at(0);
      tmp2.at(1)+=tmp3.at(1);

      time_history.push_back(tmp2);

      //pos1
      int i=0;

      while(i<leakage_current_alpha.size())
      {
        if(i>=1) G_i_tmp += leakage_current_alpha.at(i).at(0) * exp(-(time_history.back().at(0)-time_history.at(i-1).at(0))) + leakage_current_alpha.at(i).at(1) - leakage_current_alpha.at(i).at(2) * log(time_history.back().at(1)-time_history.at(i-1).at(1)) ;
        else G_i_tmp += leakage_current_alpha.at(i).at(0) * exp(-time_history.back().at(0)) + leakage_current_alpha.at(i).at(1) - leakage_current_alpha.at(i).at(2) * log(time_history.back().at(1)) ;
        i++;
      }
  }

//////////////////////////////


  G_i.push_back(G_i_tmp);
  // if(G_i_tmp/(totalDose*1.0e6+(double)phi+0.001)>5e-16)alpha_vec.push_back(1e-17);    //totalDose was devided by 1e6 to fit into a long double
  // else alpha_vec.push_back(G_i_tmp/(totalDose*1.0e6+(double)phi+0.001));

  if(G_i_tmp/(totalDose+(double)phi+0.001)>5e-16)alpha_vec.push_back(1e-17);
  else alpha_vec.push_back(G_i_tmp/(totalDose+(double)phi+0.001));

  leakage_current.push_back(G_i_tmp*1000.0*thickness*1.e-4);                                      //insert volume here for current per module
  powerconsumption.push_back(leakage_current.back()*NtoV_function(get_Neff()));

  return;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//function to read in the temp/rad profile

double iLeak_scale(double iLeak, double T, double Tref, double bandGap){
    // cout << (Tref-T) << ": " << iLeak << " ~ " << iLeak*(T*T/(Tref*Tref)*exp(-(bandGap/(2*8.6173303e-5))*(1/T - 1/Tref))) << endl;
    return iLeak*(T*T/(Tref*Tref)*exp(-(bandGap/(2*8.6173303e-5))*(1/T - 1/Tref)));
}

vector<DataElement> get_Profile(string filename)                                 // reads in a file of format (timestep/temperature/radiation_dose_rate) and gives back a vector of DataElements containing each the timestep/temperature/dose_rate information of one line of the input file
{
    ifstream input_file(filename.c_str());                                              // open stream to data file
    vector<DataElement> Temperature_Profile_vector;
    DataElement tmp_data;
    int lineNumber=0;
    if (input_file.is_open())                                                           // check whether file opening was succesfull
    {
        if(debug) cout<< "Reading temperature/radiation file: " << filename << endl;

        int timestep;
        float temp1;
        long int rad1;
        float cpTemp;
        float iLeak;
        float diLeak;
        float fillNum;
        while(true)
        {
            input_file >> timestep >> temp1 >> rad1 >> cpTemp >> iLeak >> diLeak;//; >> fillNum;

            //if(temp1<294 && temp1>290) //TODO: remove!
            //{
            //  cout << "modifying temperature from: " << temp1 << endl;
            //  //temp1=temp1+5;
            //}
            //else //temp1=temp1+5;
            //if (rad1!=0)temp1=temp1+5;

            tmp_data.duration=timestep;
            tmp_data.dose_rate=rad1*DoseRateScaling;
            tmp_data.temperature=cpTemp+deltaT(tmp_data.dose_rate);
	    //tmp_data.temperature=temp1;
            tmp_data.coolPipeTemp=cpTemp;
            tmp_data.iLeak_data=iLeak_scale(iLeak, userTref, tmp_data.temperature, bandGap);
            tmp_data.diLeak_data=iLeak_scale(diLeak, userTref, tmp_data.temperature, bandGap);

            Temperature_Profile_vector.push_back(tmp_data);

            if(debug) cout << "Line read: " << timestep << " " << temp1 << " " << rad1  << endl;

            if(input_file.eof())break;                          	// if end of file is reached: brake the loop

            limit_file_size_input--;
            if(limit_file_size_input==0)break;
        }
    }

    else
    {
        if(debug) cout << "Error opening the temperature/radiation file!"<<endl;
        return Temperature_Profile_vector;
    }

    input_file.close();

    return Temperature_Profile_vector;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//function to plot different informations in root and pdf files

void plot_vectors(vector<double> vecx, vector<double> vecy, string yName, string xName, long double totalDose, string plotname, int mode, string output_file_name)
{
    string output_file_name_pdf=output_file_name;
    output_file_name_pdf.insert(output_file_name_pdf.size(),"_");
    output_file_name_pdf.insert(output_file_name_pdf.size(),plotname);
    output_file_name_pdf.insert(output_file_name_pdf.size(),".pdf");

    string output_file_name_root=output_file_name;
    output_file_name_root.insert(output_file_name_root.size(),".root");


    TFile *file = new TFile(output_file_name_root.c_str(),"UPDATE");
    //file->cd();

    if(mode==1 || mode==2 || mode==66 || mode==8 || mode==888 || mode==999) //in case of plot as a function of time, convert the days from simulation into a more handy date
    {
        for(int k=0; k<vecx.size();k++)
        {
            int day=vecx.at(k);
            //--------------------------------------
            int seconds = (vecx.at(k) - day)*86400;

            int hours = seconds/3600;

            seconds-=hours*3600;

            int minutes = seconds/60;

            seconds-=minutes*60;

            // -------------------------------------
            boost::posix_time::ptime t2 = t1 + boost::posix_time::hours(hours+day*24)
                                             + boost::posix_time::minutes(minutes)
                                             + boost::posix_time::seconds(seconds);
            // date d2 = d + days(day);
            boost::gregorian::date d2 = t2.date();
            TDatime da (d2.year(), d2.month(), d2.day(), t2.time_of_day().hours(), t2.time_of_day().minutes(), t2.time_of_day().seconds());
            // TDatime da (d2.year(), d2.month(), d2.day(),0,0,0);

            vecx.at(k)=da.Convert();
        }
    }
    else  //otherwise the axis is fluence and we need to reconvert it since for safty reasons it was stored devided by 1e6 in order not to overflow the long double
    {
        // for(int k=0; k<vecx.size();k++)
        // {
        //    vecx.at(k)*=1.0e6;
        // }
    }

    if(mode==8 || mode==9) //in case of 8, change from unit area to unit surface for the leakage current (IMPLEMENT USER TREF HERE)
    {
        for(int k=0; k<vecy.size();k++)
        {
            vecy.at(k)/=thickness*1.0e-4; //convert from surface to volume normalisation
            //vecy.at(k)*=(userTref*userTref/(294.15*294.15)*exp(-6500*(1/userTref - 1/294.15))); //scale to user defined temperature
            vecy.at(k)*=(userTref*userTref/(294.15*294.15)*exp(-(bandGap/(2*8.6173303e-5))*(1/userTref - 1/294.15))); //scale to user defined temperature
        }
    }

    if(mode==888 || mode==889) //in case of 8, change from unit area to unit surface for the leakage current (IMPLEMENT USER TREF HERE)
    {
      // if(mode==888){
      //   for(int k=0; k<vecx.size();k++)
      //   {
      //       int day=vecx.at(k);
      //       //--------------------------------------
      //       int seconds = (vecx.at(k) - day)*86400;
      //
      //       int hours = seconds/3600;
      //
      //       seconds-=hours*3600;
      //
      //       int minutes = seconds/60;
      //
      //       seconds-=minutes*60;
      //       // -------------------------------------
      //       date d2 = d + days(day);
      //
      //       TDatime da (d2.year(), d2.month(), d2.day(), hours, minutes, seconds);
      //
      //       // TDatime da (d2.year(), d2.month(), d2.day(),0,0,0);
      //
      //       vecx.at(k)=da.Convert();
      //   }
      // }
      for(int k=0; k<vecy.size();k++)
      {
          // cout << k << ": " << vecx.at(k) << " -> " << vecy.at(k) << endl;
          vecy.at(k)*=6.48*1.62; //convert from surface to volume normalisation
          //vecy.at(k)*=(userTref*userTref/(294.15*294.15)*exp(-6500*(1/userTref - 1/294.15))); //scale to user defined temperature
          vecy.at(k)*=(userTref*userTref/(294.15*294.15)*exp(-(bandGap/(2*8.6173303e-5))*(1/userTref - 1/294.15))); //scale to user defined temperature
      }
    }
    if(mode==999 || mode==9999) //in case of 8, change from unit area to unit surface for the leakage current (IMPLEMENT USER TREF HERE)
    {
      // if(mode==999){
      //   for(int k=0; k<vecx.size();k++)
      //   {
      //       int day=vecx.at(k);
      //       //--------------------------------------
      //       int seconds = (vecx.at(k) - day)*86400;
      //
      //       int hours = seconds/3600;
      //
      //       seconds-=hours*3600;
      //
      //       int minutes = seconds/60;
      //
      //       seconds-=minutes*60;
      //       // -------------------------------------
      //       date d2 = d + days(day);
      //
      //       TDatime da (d2.year(), d2.month(), d2.day(), hours, minutes, seconds);
      //
      //       // TDatime da (d2.year(), d2.month(), d2.day(),0,0,0);
      //
      //       vecx.at(k)=da.Convert();
      //   }
      // }
      // for(int k=0; k<vecy.size();k++)
      // {
      //     // vecy.at(k)*=6.48*1.62; //convert from surface to volume normalisation
      //     //vecy.at(k)*=(userTref*userTref/(294.15*294.15)*exp(-6500*(1/userTref - 1/294.15))); //scale to user defined temperature
      //     vecy.at(k)*=(userTref*userTref/(294.15*294.15)*exp(-(bandGap/(2*8.6173303e-5))*(1/userTref - 1/294.15))); //scale to user defined temperature
      // }
    }
    const TVectorD t_timevector(vecx.size(),&vecx[0]);
    const TVectorD t_vdepvector(vecy.size(),&vecy[0]);
    const TVectorD t_time_error(vecx.size());
    const TVectorD t_vdep_error(vecy.size());

    TCanvas * c1 = new TCanvas("TCanvas_name","Title",0,0,1024,668);
    gStyle->SetOptTitle(0);

    TGraphErrors *gr = new TGraphErrors(t_timevector,t_vdepvector,t_time_error,t_vdep_error);

    if(mode==1 || mode==2 || mode==66 || mode==8 || mode==888 || mode==999) //in case of plot as a function of time, convert the days from simulation into a more handy date
    {
        gr->GetXaxis()->SetRangeUser(1488326400.,1577836799.);
        gr->GetXaxis()->SetTimeDisplay(1);
        gr->GetXaxis()->SetNdivisions(6,2,0);
        gr->GetXaxis()->SetTimeFormat("%d/%m/%Y");
        gr->GetXaxis()->SetTimeOffset(0,"gmt");
    }

    gr->GetXaxis()->SetTitle(xName.c_str());
    gr->GetYaxis()->SetTitle(yName.c_str());

    // if(mode!=66)
    // {
    //   gr->SetName(plotname.c_str());
    //   gr->Draw("AP");
    //   if(plot_pdf) c1->Print(output_file_name_pdf.c_str());
    //   gr->Write();
    // }

    gr->SetName(plotname.c_str());
    gr->Draw("AP");
    if(plot_pdf) c1->Print(output_file_name_pdf.c_str());
    gr->Write();

    if(mode==1)
    {
        cout << "Final " << plotname << " is: " << vecy.at(vecy.size()-1)<< yName <<endl;
        // cout << "Total collected fluence is: " << totalDose*1.0e6/global_layer_conversion << "fb-1 (only one conversion factor supported - if more than one was used for the profile, this value will be wrong)" << endl;
        // cout << "Total collected dose is: " << totalDose*1.0e6 << "neq/cm2" << endl;
        cout << "Total collected fluence is: " << totalDose/global_layer_conversion << "fb-1 (only one conversion factor supported - if more than one was used for the profile, this value will be wrong)" << endl;
        cout << "Total collected dose is: " << totalDose << "neq/cm2" << endl;
    }
    else
    {
      // cout << "Final " << plotname << " is: " << vecy.at(vecy.size()-1)<<endl;
    }

    file->Close();
    delete c1;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {

    if(argc!=3)                                                 		// Check for correct number of arguments
    {
        cout << "Wrong number of arguments! There are: " << argc << " arguments." << endl;
        return -1;
    }

    string output_file_name = argv[2];

    if(!overwrite_files)
    {
        time_t current_time = time(nullptr);
        // output_file_name.insert(output_file_name.size(),to_string(current_time));
    }

    string profilinput_file_name = argv[1];

    leakage_current_consts LeCo = { 1.23e-17,   //alpha_1
                                    7.07e-17,   //alpha_0_star
                                    3.3e-18,    //xi
                                    1.2e13,     //k01
                                    1.11,       //E1
                                    1.3};       //E1_star

    Annealing_constants AnCo( 1.4e-2,           //atof(argv[2]),       //gA
                              6.7e-2,           //atof(argv[3]),       //gY
                              // 0.7e-2,           //atof(argv[4]),       //gC
                              0.0061566116294,           //atof(argv[4]),       //gC
                              2.4e13,           //ka
                              7.4e14,           //ky
                              1.09,             //Ea
                              1.325,            //Ey
                              6.4118e-14);      //cc

    //DoseRateScaling = atof(argv[5]) ;

fluence_2018_begin = nconst_begin_2018/gC_lin
    vector<DataElement> temp_rad_profile = get_Profile(profilinput_file_name);    //read in the input data profile, a vector containing the individual data elements (duration, temperature, dose rate) is returned

    int max_steps=temp_rad_profile.size();
    // int max_steps=100;
    vector<double> V_dep_vector;
    vector<double> Neff_vector;
    vector<double> Ndonor_vector;
    vector<double> Nacceptor_vector;
    vector<double> Temperature_vector;
    vector<double> cpTemperature_vector;

    vector<double> time_vector;
    vector<double> fluence_vector;
    vector<double> flux_vector;

    vector<double> time_vector_data;
    vector<double> fluence_vector_data;
    vector<double> iLeak_data;
    vector<double> diLeak_data;

    vector<double> N_benef_anneal_g1_vec;
    vector<double> N_revers_anneal_g1_vec;
    vector<double> N_nadefects_g1_vec;
    vector<double> N_constdamage_g1_vec;
    vector<double> N_donor_vec;
    // vector<double> N_donor_stable_donorremoval;

    Sensor *sensor = new Sensor();

    long double time=0.0;
    long double totalDose=0.;

    // main loop where all irradiation steps happen

    cout << "Profile succesfully read. Length of the profile is: " << max_steps << endl;

    // TODO: Why is this code commented out?
    // double gA_float_min = 0.4e-2;
    // double gA_float_max = 1.4e-2;
    // double gA_float = gA_float_min;
    // double gA_step = 0.1e-2;
    // int gA_Nsteps = int((gA_float_max-gA_float_min)/gA_step);
    // double gC_float_min = 0.6e-2;
    // double gC_float_max = 2.0e-2;
    // double gC_float = gC_float_min;
    // double gC_step = 0.1e-2;
    // int gC_Nsteps = int((gC_float_max-gC_float_min)/gC_step);
    // TString gA_str;
    // TString gC_str;
    // TString gAC_str;
    // bool adjust_gAgC = true;

    // if (adjust_gAgC){
    //   for(int i=0; i < gA_Nsteps; i++){
    //     gC_float = gC_float_min;
    //     for(int j=0; j < gC_Nsteps; j++){
    //       Annealing_constants AnCo( gA_float,           //atof(argv[2]),       //gA
    //                                 1.6e-2,           //atof(argv[3]),       //gY
    //                                 gC_float,           //atof(argv[4]),       //gC
    //                                 2.4e13,           //ka
    //                                 7.4e14,           //ky
    //                                 1.09,             //Ea
    //                                 1.325,            //Ey
    //                                 6.4118e-14);      //cc
    //
    //       for (int t=0; t<max_steps; t++)   // iterate through the profile and irradiate the sensor at each step
    //       {
    //           sensor->set_temp(temp_rad_profile.at(t).temperature);
    //           sensor->irradiate(LeCo,AnCo,temp_rad_profile.at(t).dose_rate,temp_rad_profile.at(t).duration,totalDose);
    //
    //           Temperature_vector.push_back(temp_rad_profile.at(t).temperature);
    //           cpTemperature_vector.push_back(temp_rad_profile.at(t).coolPipeTemp);
    //           Neff_vector.push_back(sensor->get_Neff());
    //           V_dep_vector.push_back(NtoV_function(sensor->get_Neff()));
    //           // Ndonor_vector.push_back(sensor->get_Ndonor());
    //           // Nacceptor_vector.push_back(sensor->get_Nacceptor());
    //           // totalDose+=temp_rad_profile.at(t).dose_rate/1.0e3*temp_rad_profile.at(t).duration/1.0e3;  // time (seconds) * dose rate (neq/cm2/s)
    //           totalDose+=double(temp_rad_profile.at(t).dose_rate)*double(temp_rad_profile.at(t).duration);
    //           time+=temp_rad_profile.at(t).duration;
    //           // if(temp_rad_profile.at(t).iLeak_data > 0.){
    //           //   time_vector_data.push_back(time/(24.0*3600.0));
    //           //   fluence_vector_data.push_back(totalDose);
    //           //   iLeak_data.push_back(temp_rad_profile.at(t).iLeak_data/1000.);
    //           //   diLeak_data.push_back(temp_rad_profile.at(t).diLeak_data/1000.);
    //           //   flux_vector.push_back(temp_rad_profile.at(t).dose_rate);
    //           // }
    //           time_vector.push_back(time/(24.0*3600.0));  // time vector in days
    //           fluence_vector.push_back(totalDose);
    //           // cout << t << ": " << fluence_vector.at(t) << endl;
    //           // fluence_vector.push_back(totalDose*1.0e6);
    //           // if(t%(int)(max_steps/100.)==0) cout << (int)((int)t*100./max_steps) << " percent done..." << endl;
    //       }
    //       gA_str = TString("_")+TString::Itoa(int(gA_float*1e2),10)+TString("p")+TString::Itoa(int((gA_float*1e2-int(gA_float*1e2))*10),10);
    //       // gA_str.Replace(1,1,"p");
    //       gC_str = TString("_")+TString::Itoa(int(gC_float*1e2),10)+TString("p")+TString::Itoa(int((gC_float*1e2-int(gC_float*1e2))*10),10);
    //       // gC_str.Replace(1,1,"p");
    //       plot_vectors(time_vector,V_dep_vector,"U_{dep} [V]","Date [days]",totalDose,(TString("U_dep")+gA_str+gC_str).Data(),66,output_file_name);
    //       sensor = new Sensor();
    //       V_dep_vector.clear();
    //       Neff_vector.clear();
    //       Ndonor_vector.clear();
    //       Nacceptor_vector.clear();
    //       Temperature_vector.clear();
    //       cpTemperature_vector.clear();
    //       time_vector.clear();
    //       fluence_vector.clear();
    //       flux_vector.clear();
    //       time_vector_data.clear();
    //       fluence_vector_data.clear();
    //       iLeak_data.clear();
    //       diLeak_data.clear();
    //       totalDose = 0;
    //       time = 0;
    //       gC_float += gC_step;
    //     }
    //     gA_float += gA_step;
    //   }
    // }

    for (int t=0; t<max_steps; t++)   // iterate through the profile and irradiate the sensor at each step
    {
        sensor->set_temp(temp_rad_profile.at(t).temperature);
        sensor->irradiate(LeCo,AnCo,temp_rad_profile.at(t).dose_rate,temp_rad_profile.at(t).duration,totalDose);

        Temperature_vector.push_back(temp_rad_profile.at(t).temperature);
        cpTemperature_vector.push_back(temp_rad_profile.at(t).coolPipeTemp);
        Neff_vector.push_back(sensor->get_Neff());
        V_dep_vector.push_back(NtoV_function(sensor->get_Neff()));
        Ndonor_vector.push_back(sensor->get_Ndonor());
        Nacceptor_vector.push_back(sensor->get_Nacceptor());
        // cout << AnCo.get_gY(temp_rad_profile.at(t).temperature);

        N_benef_anneal_g1_vec.push_back(sensor->get_Nbenef_anneal_g1());
        N_revers_anneal_g1_vec.push_back(sensor->get_Nrevers_anneal_g1());
        N_nadefects_g1_vec.push_back(sensor->get_Nnadefects_g1());
        N_constdamage_g1_vec.push_back(sensor->get_Nconstdamage_g1());
        N_donor_vec.push_back(sensor->get_Ndonor());
        // N_donor_stable_donorremoval.push_back(NtoV_function(sensor->get_Ndonor_stable_donorremoval()));
        // 14722058.0
        // 16450058.0
        // 45480458.0
        // 48504458.0
        // totalDose+=temp_rad_profile.at(t).dose_rate/1.0e3*temp_rad_profile.at(t).duration/1.0e3;  // time (seconds) * dose rate (neq/cm2/s)
        totalDose+=double(temp_rad_profile.at(t).dose_rate)*double(temp_rad_profile.at(t).duration);
        time+=temp_rad_profile.at(t).duration;
        // cout << temp_rad_profile.at(t).duration << endl;
        // if(temp_rad_profile.at(t).iLeak_data > 0.){
        if(temp_rad_profile.at(t).iLeak_data > 0. && temp_rad_profile.at(t).duration == 1200 
            && (time < 14722058 || (time > 16450058 && time < 45480458) || time > 48504458)){
          // cout << time << endl;
          time_vector_data.push_back(time/(24.0*3600.0));
          fluence_vector_data.push_back(totalDose/6.);
          iLeak_data.push_back(temp_rad_profile.at(t).iLeak_data/1000.);
          diLeak_data.push_back(temp_rad_profile.at(t).diLeak_data/1000.);
          flux_vector.push_back(temp_rad_profile.at(t).dose_rate);
        }
        time_vector.push_back(time/(24.0*3600.0));  // time vector in days
        fluence_vector.push_back(totalDose);
        // cout << t << ": " << fluence_vector.at(t) << endl;
        // fluence_vector.push_back(totalDose*1.0e6);
        // if(t%(int)(max_steps/100.)==0) cout << (int)((int)t*100./max_steps) << " percent done..." << endl;
    }
    // data output, plotting and visualisation

    cout << "Processing finished, writing data..." << endl;

    // plots as function of time
    plot_vectors(time_vector,Neff_vector,"N_{eff} [1/cm^{3}]","Date [days]",totalDose,"Neff",1,output_file_name);
    plot_vectors(time_vector,Ndonor_vector,"N_{donor} [1/cm^{3}]","Date [days]",totalDose,"Ndonors",2,output_file_name);
    plot_vectors(time_vector,Nacceptor_vector,"N_{acceptor} [1/cm^{3}]","Date [days]",totalDose,"Nacceptors",2,output_file_name);
    plot_vectors(time_vector,V_dep_vector,"U_{dep} [V]","Date [days]",totalDose,"U_dep",66,output_file_name);

    plot_vectors(time_vector,N_benef_anneal_g1_vec,"N_{dep,benef_anneal_g1} [V]","Date [days]",totalDose,"N_benef_anneal_g1",66,output_file_name);
    plot_vectors(time_vector,N_revers_anneal_g1_vec,"N_{dep,revers_anneal_g1} [V]","Date [days]",totalDose,"N_revers_anneal_g1",66,output_file_name);
    plot_vectors(time_vector,N_nadefects_g1_vec,"N_{dep,nadefects_g1} [V]","Date [days]",totalDose,"N_nadefects_g1",66,output_file_name);
    plot_vectors(time_vector,N_constdamage_g1_vec,"N_{dep,constdamage_g1} [V]","Date [days]",totalDose,"N_constdamage_g1",66,output_file_name);
    plot_vectors(time_vector,N_donor_vec,"N_{dep,donor} [V]","Date [days]",totalDose,"N_donor",66,output_file_name);
    // plot_vectors(time_vector,N_donor_stable_donorremoval,"N_{dep,donor_stable_donorremoval} [V]","Date [days]",totalDose,"N_donor_stable_donorremoval",66,output_file_name);

    plot_vectors(time_vector,sensor->get_alpha_vec(),"#alpha [A/cm]","Date [days]",totalDose,"alpha",2,output_file_name);
    plot_vectors(time_vector,sensor->get_leakage_current(),"I_{leak} (@21C) [mA/cm^{2}]","Date [days]",totalDose,"I_leak",2,output_file_name);
    // plot_vectors(time_vector,sensor->get_leakage_current(),"I_{leak} (@ -7C) [mA] per module","Date [days]",totalDose,"I_leak_per_module",888,output_file_name);
    plot_vectors(time_vector,sensor->get_leakage_current(),"I_{leak} (@ -7C) [mA], 1 ROG","Date [days]",totalDose,"I_leak_per_module",888,output_file_name);
    plot_vectors(time_vector,sensor->get_leakage_current(),"I_{leak} (@userTref) [mA/cm^{3}]","Date [days]",totalDose,"I_leak_volume",8,output_file_name);
    plot_vectors(time_vector,sensor->get_G_i(),"G_{i} [A/cm^{3}]","Date [days]",totalDose,"G_i",2,output_file_name);
    plot_vectors(time_vector,sensor->get_powerconsumption(),"P [mW/cm^{2}]","Date [days]",totalDose,"powerconsumption",2,output_file_name);
    plot_vectors(time_vector,Temperature_vector,"Temperature [K]","Date [days]",totalDose,"temperature",2,output_file_name);
    plot_vectors(time_vector,cpTemperature_vector,"Temperature [K]","Date [days]",totalDose,"cp_Temperature",2,output_file_name);
    plot_vectors(time_vector,fluence_vector,"Fluence [n_{eq}/cm^{2}]","Date [days]",totalDose,"fluence",2,output_file_name);
    plot_vectors(time_vector_data,flux_vector,"Fluence [n_{eq}/cm^{2}/s]","Date [days]",totalDose,"flux",2,output_file_name);
    plot_vectors(time_vector_data,iLeak_data,"I_{leak} (@ -7C) [mA], 1 ROG","Date [days]",totalDose,"I_leak_per_module_data",999,output_file_name);
    plot_vectors(time_vector_data,diLeak_data,"dI_{leak} (@ -8.5C) [mA] per module","Date [days]",totalDose,"dI_leak_per_module_data",999,output_file_name);

    // plots as function of dose
    plot_vectors(fluence_vector,Neff_vector,"N_{eff} [1/cm^{3}]","Fluence [n_{eq}/cm^{2}]",totalDose,"Neff_vs_fluence",3,output_file_name);
    plot_vectors(fluence_vector,Ndonor_vector,"N_{donor} [1/cm^{3}]","Fluence [n_{eq}/cm^{2}]",totalDose,"Ndonors_vs_fluence",4,output_file_name);
    plot_vectors(fluence_vector,Nacceptor_vector,"N_{acceptor} [1/cm^{3}]","Fluence [n_{eq}/cm^{2}]",totalDose,"Nacceptors_vs_fluence",4,output_file_name);
    plot_vectors(fluence_vector,V_dep_vector,"U_{dep} [V]","Fluence [n_{eq}/cm^{2}]",totalDose,"U_dep_vs_fluence",4,output_file_name);
    plot_vectors(fluence_vector,sensor->get_alpha_vec(),"#alpha [A/cm]","Fluence [n_{eq}/cm^{2}]",totalDose,"alpha_vs_fluence",4,output_file_name);
    plot_vectors(fluence_vector,sensor->get_leakage_current(),"I_{leak} [mA/cm^{2}]","Fluence [n_{eq}/cm^{2}]",totalDose,"I_leak_vs_fluence",4,output_file_name);
    plot_vectors(fluence_vector,sensor->get_leakage_current(),"I_{leak} (@21C) [mA] per module","Date [days]",totalDose,"I_leak_per_module_vs_fluence",889,output_file_name);
    plot_vectors(fluence_vector,sensor->get_leakage_current(),"I_{leak} (@userTref) [mA/cm^{3}]","Fluence [n_{eq}/cm^{2}]",totalDose,"I_leak_volume_vs_fluence",9,output_file_name);
    plot_vectors(fluence_vector,sensor->get_G_i(),"G_{i} [A/cm^{3}]","Fluence [n_{eq}/cm^{2}]",totalDose,"G_i_vs_fluence",4,output_file_name);
    plot_vectors(fluence_vector,sensor->get_powerconsumption(),"P [mW/cm^{2}]","Fluence [n_{eq}/cm^{2}]",totalDose,"powerconsumption_vs_fluence",4,output_file_name);
    plot_vectors(fluence_vector,Temperature_vector,"Temperature [K]","Fluence [n_{eq}/cm^{2}",totalDose,"temperature_vs_fluence",4,output_file_name);
    plot_vectors(fluence_vector_data,iLeak_data,"I_{leak} (@ -7C) [mA] per module","Date [days]",totalDose,"I_leak_per_module_data_vs_fluence",9999,output_file_name);
    plot_vectors(fluence_vector_data,diLeak_data,"dI_{leak} (@ -8.5C) [mA] per module","Date [days]",totalDose,"dI_leak_per_module_data_vs_fluence",9999,output_file_name);

    return 0;
}
