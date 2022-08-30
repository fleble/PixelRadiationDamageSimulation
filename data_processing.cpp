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
//#include <TH2D.h>           		// root stuff for file reading
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

using namespace std;
using namespace boost::gregorian;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool use_CMS_paper_alpha0 = false; //if false: temperature averaging for leakage current calculation a la Moll (theta function for time scale change), if true: use CMS approach of temperature average (attention: constants can only be changed directly in the code in this case)

bool debug=false;                                    // additional debugging console outputs
bool plot_pdf=true;                                  // plots are by default saved as root file, but can also be directly exported as pdf
bool overwrite_files=true;                           // false = append internal clock identifier to file name, false = overwrite pdfs, append root files // files are by default created with a unique name depending on the internal clock and you can switch it of here such that pdf files are overwritten and the root file is extended (old plots will stay in the root file as well)
double timestep=1;                                   // step size is 1 second, do not change this!
double donorremovalfraction=0.99;                    // fraction of donors which can be removed by donor removal
double userTref=273.15;                              // set a reference temperature for the volume corrected leakage current plot (it will only effect this one plot!) Now: implemented!
double bandGap=1.21;                                 // eV used for scaling temperatures

string output_file_name = "simulation_results";      // set unique file name for each simulation

date d(2015,May,1);                                  // IBL     //set a date for the plots to begin (to be correct it has to be equal to the beginning of your (!) temp/irr profile)
//date d(2011,Feb,11);                               // PIXEL


const double Ndonor_0 = 1.7e12;                      // IBL       // initial donor concentration (right now the code only works for originally n-type bulk material!)
//const double Ndonor_0 = 1.4e12;                    // B-Layer Layer1/2 Disks

double thickness=200;                                // IBL       // sensor thickness
//double thickness=250;                              // B-Layer Layer1/2 Disks

//double global_layer_conversion=0.92e12;            // Layer 1    //conversion factor from luminosity to neq (1fb-1=2.3e12neq/cm2) - is used only for computation of total luminosity, will be wrong if there are multiple conversion factor in the original profiles as for example due to different center of mass energies
//double global_layer_conversion=1.1e12;             // L1 average
//double global_layer_conversion=0.571e12;           // Layer 2
//double global_layer_conversion=0.7e12;             // Layer 2 average
//double global_layer_conversion=0.582e12;           // Disks
double global_layer_conversion=6.262e12;             // IBL
//double global_layer_conversion=2.8e12;//2.929e12;  // B-Layer

int limit_file_size_input = 0;                       // how many lines should be read from the profile, 0 for everything

float DoseRateScaling=1.;                            // this parameter is multiplied to all fluence rates from the input profile

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Do not change things below this line without thinking carefully!

//Definition of Hamburg model depletion voltage functions

TF1 *Naccept_rev_TF1 = new TF1("Naccept_rev_TF1","[0]*[1]*(1-TMath::Exp(-[2]*x))/[2] + [3]*TMath::Exp(-[2]*x)",0,10000000);
TF1 *Naccept_const_TF1 = new TF1("Naccept_const_TF1","[0] * [1]*x",0,10000000);
TF1 *Nneutrals_rev_TF1 = new TF1("Ndonor_rev_TF1","[0]*[1]*(1-TMath::Exp(-[2]*x))/[2] + [3]*TMath::Exp(-[2]*x)",0,10000000);
TF1 *Ndonor_const_TF1 = new TF1("Ndonor_const_TF1","-[0]*(1-TMath::Exp(-[1]*[2]*x))",0,10000000);

TF1 *Naccept_rev_TF1_approx = new TF1("Naccept_rev_TF1","[0]*[1]*x + [3]*TMath::Exp(-[2]*x)",0,10000000);
TF1 *Ndonor_neutrals_rev_TF1_approx = new TF1("Ndonor_rev_TF1","[0]*[1]*x + [3]*TMath::Exp(-[2]*x)",0,10000000);

TF1 *Ndonor_TF1 = new TF1("Ndonor_TF1","[0]*[1]/[2] * ([2]*x+TMath::Exp(-[2]*x)-1) +[3]*(1-TMath::Exp(-[2]*x)) ",0,10000000);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct DataElement
{
    int duration;
    float temperature;
    long int dose_rate;
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//calculate depletion voltage from given effective doping concentration

double NtoV_function(double doping_conc)
{
    return 1.6021766208e-13/(2*11.68*8.854187817)*fabs(doping_conc)*thickness*thickness;                // q/2*epsilon*epsilon_0 * Neff*D^2
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


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
  Nneutral_reversible=0;
  volume= 0.135;
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

void Sensor::irradiate(leakage_current_consts leconsts, Annealing_constants constants,long int phi, float t, long double totalDose)
{
  //t=t*3600;                                                                     //conversion from hours to seconds
  double a=1e-30;
//calculating the effective doping concentration

if(debug) cout << t<<" "<< constants.get_gA(Temperature)<<" "<<phi<<" "<<constants.get_ka(Temperature)<<" "<<Nacceptor<<" "<< constants.get_gC(Temperature)<<" "<< constants.get_gY(Temperature)<<" "<< constants.get_ky1(Temperature) <<" "<<Ndonor <<endl;

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
      tmp2.at(1)=t / 60. * exp(-leconsts.E1_star*(1.0/(double)Temperature-1.0/294.15)/(8.6173303e-5));  //T-ref is hardcoded and fixed, it cannot be changed in a trivial way since the other parameters, especially alpha 0 star were evaluated at this reference temperature!!!    // t/60 to convert time from seconds to minutes since this is the proposed unit here

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
  if(G_i_tmp/(totalDose*1.0e6+(double)phi+0.001)>5e-16)alpha_vec.push_back(1e-17);    //totalDose was devided by 1e6 to fit into a long double
  else alpha_vec.push_back(G_i_tmp/(totalDose*1.0e6+(double)phi+0.001));
  leakage_current.push_back(G_i_tmp*1000.0*thickness*1.e-4);                                      //insert volume here for current per module
  powerconsumption.push_back(leakage_current.back()*NtoV_function(get_Neff()));

  return;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//function to read in the temp/rad profile

vector<DataElement> get_Profile(string filename)                                 // reads in a file of format (timestep/temperature/radiation_dose_rate) and gives back a vector of DataElements containing each the timestep/temperature/dose_rate information of one line of the input file
{
    ifstream input_file(filename.c_str());                                              // open stream to data file
    vector<DataElement> Temperature_Profile_vector;
    DataElement tmp_data;

    if (input_file.is_open())                                                           // check whether file opening was succesfull
    {
        if(debug) cout<< "Reading temperature/radiation file: " << filename << endl;

        int timestep;
        float temp1;
        long int rad1;

        while(true)
        {
            input_file >> timestep >> temp1 >> rad1;

            //if(temp1<294 && temp1>290) //TODO: remove!
            //{
            //  cout << "modifying temperature from: " << temp1 << endl;
            //  //temp1=temp1+5;
            //}
            //else //temp1=temp1+5;
            //if (rad1!=0)temp1=temp1+5;

            tmp_data.duration=timestep;
            tmp_data.temperature=temp1;
            tmp_data.dose_rate=rad1*DoseRateScaling;
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

void plot_vectors(vector<double> vecx, vector<double> vecy, string yName, string xName, long double totalDose, string plotname, int mode)
{
    string output_file_name_pdf=output_file_name;
    output_file_name_pdf.insert(output_file_name_pdf.size(),"_");
    output_file_name_pdf.insert(output_file_name_pdf.size(),plotname);
    output_file_name_pdf.insert(output_file_name_pdf.size(),".pdf");

    string output_file_name_root=output_file_name;
    output_file_name_root.insert(output_file_name_root.size(),".root");


    TFile *file = new TFile(output_file_name_root.c_str(),"UPDATE");
    //file->cd();

    if(mode==1 || mode==2 || mode==66 || mode==8) //in case of plot as a function of time, convert the days from simulation into a more handy date
    {
        for(int k=0; k<vecx.size();k++)
        {
            int day=vecx.at(k);

            date d2 = d + days(day);

            TDatime da (d2.year(), d2.month(), d2.day(),0,0,0);

            vecx.at(k)=da.Convert();
        }
    }
    else  //otherwise the axis is fluence and we need to reconvert it since for safty reasons it was stored devided by 1e6 in order not to overflow the long double
    {
        //for(int k=0; k<vecx.size();k++)
        //{
        //    vecx.at(k)*=1.0e6;
        //}
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

    const TVectorD t_timevector(vecx.size(),&vecx[0]);
    const TVectorD t_vdepvector(vecy.size(),&vecy[0]);
    const TVectorD t_time_error(vecx.size());
    const TVectorD t_vdep_error(vecy.size());

    TCanvas * c1 = new TCanvas("TCanvas_name","Title",0,0,1024,668);
    gStyle->SetOptTitle(0);

    TGraphErrors *gr = new TGraphErrors(t_timevector,t_vdepvector,t_time_error,t_vdep_error);

    if(mode==1 || mode==2 || mode==66 || mode==8) //in case of plot as a function of time, convert the days from simulation into a more handy date
    {
        gr->GetXaxis()->SetTimeDisplay(1);
        gr->GetXaxis()->SetNdivisions(6,2,0);
        gr->GetXaxis()->SetTimeFormat("%d/%m/%Y");
        gr->GetXaxis()->SetTimeOffset(0,"gmt");
    }

    gr->GetXaxis()->SetTitle(xName.c_str());
    gr->GetYaxis()->SetTitle(yName.c_str());

    if(mode!=66)
    {
      gr->SetName(plotname.c_str());
      gr->Draw("AP");
      if(plot_pdf) c1->Print(output_file_name_pdf.c_str());
      gr->Write();
    }

    if(mode==1)
    {
        cout << "Final " << plotname << " is: " << vecy.at(vecy.size()-1)<< yName <<endl;
        cout << "Total collected fluence is: " << totalDose*1.0e6/global_layer_conversion << "fb-1 (only one conversion factor supported - if more than one was used for the profile, this value will be wrong)" << endl;
        cout << "Total collected dose is: " << totalDose*1.0e6 << "neq/cm2" << endl;
    }
    else
    {
      cout << "Final " << plotname << " is: " << vecy.at(vecy.size()-1)<<endl;
    }

    file->Close();
    delete c1;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {

    if(argc!=2)                                                 		// Check for correct number of arguments
    {
        cout << "Wrong number of arguments! There are: " << argc << " arguments." << endl;
        return -1;
    }

    if(!overwrite_files)
    {
        time_t current_time = time(nullptr);
        output_file_name.insert(output_file_name.size(),to_string(current_time));
    }

    string profilinput_file_name = argv[1];

    leakage_current_consts LeCo = { 1.23e-17,   //alpha_1
                                    7.07e-17,   //alpha_0_star
                                    3.3e-18,    //beta
                                    1.2e13,     //k01
                                    1.11,       //E1
                                    1.3};       //E1_star

    Annealing_constants AnCo( 1.0e-2,           //atof(argv[2]),       //gA
                              1.6e-2,           //atof(argv[3]),       //gY
                              1.0e-2,           //atof(argv[4]),       //gC
                              2.4e13,           //ka
                              7.4e14,           //ky
                              1.09,             //Ea
                              1.325,            //Ey
                              6.4118e-14);      //cc

    //DoseRateScaling = atof(argv[5]) ;


    vector<DataElement> temp_rad_profile = get_Profile(profilinput_file_name);    //read in the input data profile, a vector containing the individual data elements (duration, temperature, dose rate) is returned

    int max_steps=temp_rad_profile.size();

    vector<double> V_dep_vector;
    vector<double> Neff_vector;
    vector<double> Ndonor_vector;
    vector<double> Nacceptor_vector;
    vector<double> Temperature_vector;

    vector<double> time_vector;
    vector<double> fluence_vector;

    Sensor *sensor = new Sensor();

    long double time=0.0;
    long double totalDose=0;

    // main loop where all irradiation steps happen

    cout << "Profile succesfully read. Length of the profile is: " << max_steps << endl;

    for (int t=0; t<max_steps; t++)   // iterate through the profile and irradiate the sensor at each step
    {
        sensor->set_temp(temp_rad_profile.at(t).temperature);
        sensor->irradiate(LeCo,AnCo,temp_rad_profile.at(t).dose_rate,temp_rad_profile.at(t).duration,totalDose);

        Temperature_vector.push_back(temp_rad_profile.at(t).temperature);

        Neff_vector.push_back(sensor->get_Neff());
        V_dep_vector.push_back(NtoV_function(sensor->get_Neff()));
        Ndonor_vector.push_back(sensor->get_Ndonor());
        Nacceptor_vector.push_back(sensor->get_Nacceptor());

        totalDose+=temp_rad_profile.at(t).dose_rate/1.0e3*temp_rad_profile.at(t).duration/1.0e3;  // time (seconds) * dose rate (neq/cm2/s)

        time+=temp_rad_profile.at(t).duration;
        time_vector.push_back(time/(24.0*3600.0));  // time vector in days
        fluence_vector.push_back(totalDose*1.0e6);

        if(t%(int)(max_steps/100.)==0) cout << (int)((int)t*100./max_steps) << " percent done..." << endl;
    }

    // data output, plotting and visualisation

    cout << "Processing finished, writing data..." << endl;

    // plots as function of time
    plot_vectors(time_vector,Neff_vector,"N_{eff} [1/cm^{3}]","Date [days]",totalDose,"Neff",1);
    plot_vectors(time_vector,Ndonor_vector,"N_{donor} [1/cm^{3}]","Date [days]",totalDose,"Ndonors",2);
    plot_vectors(time_vector,Nacceptor_vector,"N_{acceptor} [1/cm^{3}]","Date [days]",totalDose,"Nacceptors",2);
    plot_vectors(time_vector,V_dep_vector,"U_{dep} [V]","Date [days]",totalDose,"U_dep",66);
    plot_vectors(time_vector,sensor->get_alpha_vec(),"#alpha [A/cm]","Date [days]",totalDose,"alpha",2);
    plot_vectors(time_vector,sensor->get_leakage_current(),"I_{leak} (@21C) [mA/cm^{2}]","Date [days]",totalDose,"I_leak",2);
    plot_vectors(time_vector,sensor->get_leakage_current(),"I_{leak} (@userTref) [mA/cm^{3}]","Date [days]",totalDose,"I_leak_volume",8);
    plot_vectors(time_vector,sensor->get_G_i(),"G_{i} [A/cm^{3}]","Date [days]",totalDose,"G_i",2);
    plot_vectors(time_vector,sensor->get_powerconsumption(),"P [mW/cm^{2}]","Date [days]",totalDose,"powerconsumption",2);
    plot_vectors(time_vector,Temperature_vector,"Temperature [K]","Date [days]",totalDose,"temperature",2);
    plot_vectors(time_vector,fluence_vector,"Fluence [n_{eq}/cm^{2}]","Date [days]",totalDose,"fluence",2);

    // plots as function of dose
    plot_vectors(fluence_vector,Neff_vector,"N_{eff} [1/cm^{3}]","Fluence [n_{eq}/cm^{2}]",totalDose,"Neff_vs_fluence",3);
    plot_vectors(fluence_vector,Ndonor_vector,"N_{donor} [1/cm^{3}]","Fluence [n_{eq}/cm^{2}]",totalDose,"Ndonors_vs_fluence",4);
    plot_vectors(fluence_vector,Nacceptor_vector,"N_{acceptor} [1/cm^{3}]","Fluence [n_{eq}/cm^{2}]",totalDose,"Nacceptors_vs_fluence",4);
    plot_vectors(fluence_vector,V_dep_vector,"U_{dep} [V]","Fluence [n_{eq}/cm^{2}]",totalDose,"U_dep_vs_fluence",4);
    plot_vectors(fluence_vector,sensor->get_alpha_vec(),"#alpha [A/cm]","Fluence [n_{eq}/cm^{2}]",totalDose,"alpha_vs_fluence",4);
    plot_vectors(fluence_vector,sensor->get_leakage_current(),"I_{leak} [mA/cm^{2}]","Fluence [n_{eq}/cm^{2}]",totalDose,"I_leak_vs_fluence",4);
    plot_vectors(fluence_vector,sensor->get_leakage_current(),"I_{leak} (@userTref) [mA/cm^{3}]","Fluence [n_{eq}/cm^{2}]",totalDose,"I_leak_volume_vs_fluence",9);
    plot_vectors(fluence_vector,sensor->get_G_i(),"G_{i} [A/cm^{3}]","Fluence [n_{eq}/cm^{2}]",totalDose,"G_i_vs_fluence",4);
    plot_vectors(fluence_vector,sensor->get_powerconsumption(),"P [mW/cm^{2}]","Fluence [n_{eq}/cm^{2}]",totalDose,"powerconsumption_vs_fluence",4);
    plot_vectors(fluence_vector,Temperature_vector,"Temperature [K]","Fluence [n_{eq}/cm^{2}",totalDose,"temperature_vs_fluence",4);

    return 0;
}
