// Standard libs and STL
#include <iostream>                     // Basic IO
#include <fstream>                      // Read and write files
#include <cstdlib>                      // convert arguments of function call into int, float,...
#include <string>
#include <vector>
#include <math.h>                       // for fabs()
#include <cmath>
#include <sstream>                      // to get string into stream
#include <ctime>
#include <time.h>
#include <json/value.h>
// ROOT
#include <TFile.h>                      // more root stuff
#include <TH2D.h>                       // root stuff for file reading
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
#include <TDatime.h>
// BOOST
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>


using namespace std;
namespace po = boost::program_options;
// using namespace boost::gregorian;
// using namespace boost::posix_time;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// TODO: replace all this by arguments
// use_CMS_paper_alpha0
// if false: temperature averaging for leakage current calculation a la Moll
// (theta function for time scale change)
// if true: use CMS approach of temperature average
// attention: constants can only be changed directly in the code in this case
bool use_CMS_paper_alpha0 = false;
// begin date: 2017-05-23 14:32:22.210863
bool debug = true;                                    // additional debugging console outputs
bool plot_pdf = false;                                  // plots are by default saved as root file, but can also be directly exported as pdf
bool overwrite_files = true;                           // false = append internal clock identifier to file name, false = overwrite pdfs, append root files // files are by default created with a unique name depending on the internal clock and you can switch it of here such that pdf files are overwritten and the root file is extended (old plots will stay in the root file as well)
double timestep = 1;                                   // step size is 1 second, do not change this!
double donorremovalfraction = 0.99;                    // fraction of donors which can be removed by donor removal
const double userTrefC = 0.;
double userTref = 273.15 + userTrefC;                          // set a reference temperature for the volume corrected leakage current plot (it will only effect this one plot!) Now: implemented!
double bandGap = 1.21;                                 // eV used for scaling temperatures

string startTime = "2017-05-23 14:32:22";
// date d(2017,May,23);                                  // IBL     //set a date for the plots to begin (to be correct it has to be equal to the beginning of your (!) temp/irr profile)
boost::posix_time::ptime t1(boost::posix_time::time_from_string(startTime));
//date d(2011,Feb,11);                               // PIXEL


const double Ndonor_0 = 1.7e12;                      // IBL       // initial donor concentration (right now the code only works for originally n-type bulk material!)
//const double Ndonor_0 = 1.4e12;                    // B-Layer Layer1/2 Disks

// double thickness=200;                                // IBL       // sensor thickness
//double thickness=250;                              // B-Layer Layer1/2 Disks
double thickness = 285;                                // CMS sensor thickness
//double global_layer_conversion=0.92e12;            // Layer 1    //conversion factor from luminosity to neq (1fb-1=2.3e12neq/cm2) - is used only for computation of total luminosity, will be wrong if there are multiple conversion factor in the original profiles as for example due to different center of mass energies
//double global_layer_conversion=1.1e12;             // L1 average
//double global_layer_conversion=0.571e12;           // Layer 2
//double global_layer_conversion=0.7e12;             // Layer 2 average
//double global_layer_conversion=0.582e12;           // Disks
double global_layer_conversion = 6.262e12;             // IBL
//double global_layer_conversion=2.8e12;//2.929e12;  // B-Layer

int limit_file_size_input = 0;                       // how many lines should be read from the profile, 0 for everything

// TODO: Understand this scaling, why is it "temporary"?
float DoseRateScaling = 0.0859955711871/8.9275E-02; //temporary scaling for bpix
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


po::variables_map parseArguments(int argc, char* argv[]) {
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "Produce help message")
        ("input-profile,i", po::value<string>(), "Input profile file name")
        ("output-root-file,o", po::value<string>(), "Output ROOT file name")
    ;
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    

    if (vm.count("help")) {
        cout << desc << endl;
        exit(0);
    }
 
    return vm;
}

bool checkArguments(po::variables_map vm) {
   
    if (not vm.count("input-profile"))
    {
        cout << "Error: --input-profile/-i is a mandatory argument!" << endl;
        return false;
    }
    if (not vm.count("output-root-file"))
    {
        cout << "Error: --output-root-file/-o is a mandatory argument!" << endl;
        return false;
    }
    return true;
}


struct DataElement
{
    int timestamp;
    int duration;
    float temperature;
    long int doseRate;
    float coolPipeTemp;
    float leakageCurrentData;
    float dleakageCurrentData;
};


//calculate depletion voltage from given effective doping concentration

//some problem with memory !!!!!!

double NtoV_function(double doping_conc)
{
    return 1.6021766208e-13/(2*11.68*8.854187817)*fabs(doping_conc)*thickness*thickness;                // q/2*epsilon*epsilon_0 * Neff*D^2
}

class AnnealingConstants {     // Class to store a set of annealing specific constants and compute temperature dependent values
public:
    double gA;
    double gY;
    double gC;
    double ka;
    double ky;
    double Ea;
    double Ey;
    double cc;


    AnnealingConstants(double, double, double, double, double, double, double, double);
    AnnealingConstants();
    void initialize (double, double, double, double, double, double, double, double);
    void initialize (string);
    double get_gA() const;
    double get_gY() const;
    double get_gC() const;
    double get_ka(int) const;
    double get_ky(int) const;
    double get_Ea() const;
    double get_Ey() const;
    double get_cc() const;
};

AnnealingConstants::AnnealingConstants()
{
    gA = 0.;
    gY = 0.;
    gC = 0.;
    ka = 0.;
    ky = 0.;
    Ea = 0.;
    Ey = 0.;
    cc = 0.;
}

AnnealingConstants::AnnealingConstants(double gA_, double gY_, double gC_, double ka_,
                                       double ky_, double Ea_, double Ey_, double cc_)
{
    gA = gA_;
    gY = gY_;
    gC = gC_;
    ka = ka_;
    ky = ky_;
    Ea = Ea_;
    Ey = Ey_;
    cc = cc_;
}

void AnnealingConstants::initialize(double gA_, double gY_, double gC_, double ka_,
                                      double ky_, double Ea_, double Ey_, double cc_)
{
    gA = gA_;
    gY = gY_;
    gC = gC_;
    ka = ka_;
    ky = ky_;
    Ea = Ea_;
    Ey = Ey_;
    cc = cc_;
}

void AnnealingConstants::initialize(string filename)
{
    ifstream input_file(filename.c_str());
    if (input_file.is_open())
    {
        cout << "Reading annealing constant file: " << filename << endl;

        string line;
      
        string quantity;
        string valueStr;
        float value;

        bool endOfFile = false;
        std::size_t found;
        while (not endOfFile)
        {
            if(input_file.eof()) {break;}  // if end of file is reached: brake the loop

            getline (input_file, line);
           
            found =  line.find("#");
            if (found != std::string::npos) {
                line = line.substr(0, found);
            }
            found =  line.find(":");
            if (found == std::string::npos) {continue;}
            quantity = line.substr(0, found);
            valueStr = line.substr(found+1, line.size());
            found =  valueStr.find(",");
            if (found != std::string::npos) {
                valueStr = valueStr.substr(0, found);
            }
            boost::trim(quantity);
            boost::trim(valueStr);
            value = stod(valueStr);

            if (quantity == "\"gA\"" or quantity == "'gA'") {gA = value;}
            if (quantity == "\"gY\"" or quantity == "'gY'") {gY = value;}
            if (quantity == "\"gC\"" or quantity == "'gC'") {gC = value;}
            if (quantity == "\"ka\"" or quantity == "'ka'") {ka = value;}
            if (quantity == "\"ky\"" or quantity == "'ky'") {ky = value;}
            if (quantity == "\"Ea\"" or quantity == "'Ea'") {Ea = value;}
            if (quantity == "\"Ey\"" or quantity == "'Ey'") {Ey = value;}
            if (quantity == "\"cc\"" or quantity == "'cc'") {cc = value;}
        }
    }
    else
    {
        if(debug) cout << "Error opening the temperature/radiation file!"<<endl;
    }

    input_file.close();
}

double AnnealingConstants::get_gA() const
{
    return gA;
}

double AnnealingConstants::get_gY() const
{
    return gY;
}

double AnnealingConstants::get_gC() const
{
    return gC;
}

double AnnealingConstants::get_ka(int T) const
{
    return ka*exp(-Ea/(8.6173303e-5*(double)T));
}

double AnnealingConstants::get_ky(int T) const
{
    return ky*exp(-Ey/(8.6173303e-5*(double)T));
}

double AnnealingConstants::get_Ea() const
{
    return Ea;
}

double AnnealingConstants::get_Ey() const
{
    return Ey;
}

double AnnealingConstants::get_cc() const
{
    return cc;
}


// leakage current constants will be stored in this struct

struct LeakageCurrentConstants
{
  double alpha_1;
  double alpha_0_star;
  double beta;
  double k01;
  double E1;
  double E1_star;
};

void initializeLeakageCurrentConstants(LeakageCurrentConstants& leakageCurrentConstants)
{
    leakageCurrentConstants.alpha_1 = 1.23e-17;       //  A/cm  +/- 0.06
    leakageCurrentConstants.alpha_0_star = 7.07e-17;  //  A/cm
    leakageCurrentConstants.beta = 3.29e-18;          //  A/cm  +/- 0.18
    leakageCurrentConstants.k01 = 1.2e13;             //  /s    +5.3/-1.0   (!!!)
    leakageCurrentConstants.E1 = 1.11;                //  eV    +/- 0.05
    leakageCurrentConstants.E1_star = 1.3;            //  eV    +/- 0.14
}


// all atributes of a sensor (like irradiation, doping concentration, leakage currents, etc) will be stored in this class

class Sensor
{
    private:
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
  
        int temperature;
        vector<vector<double>> leakage_current_alpha;
        vector<vector<double>> time_history;
        vector<double> G_i;
        vector<double> leakage_current;
        double volume;
        vector<double> alpha_vec;
        vector<double> powerconsumption;
        vector<double> inverse_weighted_temp;

    public:
        Sensor();
        double get_Neff() const;
        void set_Nacceptor(double);
        void set_Ndonor(double);
        double get_Nacceptor() const;
        double get_Ndonor() const;

        void set_temperature(int);
        double get_temperature() const;
        void irradiate(LeakageCurrentConstants,AnnealingConstants, long int, float, long double);
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

Sensor::Sensor()
{
  Nacceptor=0;
  Ndonor=Ndonor_0;
  temperature=0;
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

double Sensor::get_temperature() const
{
    return temperature;
}

void Sensor::set_temperature(int value)
{
    temperature = value;
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

void Sensor::set_Nacceptor(double value)
{
    Nacceptor = value;
    return;
}

void Sensor::set_Ndonor(double value)
{
    Ndonor = value;
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

void Sensor::irradiate(
    LeakageCurrentConstants leakageCurrentConstants,
    AnnealingConstants constants,
    long int fluence,
    float duration,
    long double totalDose)
{
    double a = 1e-30;
    //calculating the effective doping concentration
    // debug=1;
    // if(debug) cout << duration <<" "<< constants.get_gA()<<" "<<fluence<<" "<<constants.get_ka(temperature)<<" "<<Nacceptor<<" "<< constants.get_gC()<<" "<< constants.get_gY()<<" "<< constants.get_ky(temperature) <<" "<<Ndonor <<endl;
    // cout << endl << "Beginning irradiation..." << endl;
    // cout << "Duration: " << duration << endl;
    // cout << "Fluence: " << fluence << endl;
    // cout << "gA: " << constants.get_gA() << endl;
    // cout << "ka: " << constants.get_ka(temperature) << endl;
    // cout << "gC: " << constants.get_gC() << endl;
    // cout << "gY: " << constants.get_gY() << endl;
    // cout << "ky: " << constants.get_ky(temperature) << endl;
    // cout << "n_acceptors: " << Nacceptor << endl;
    // cout << "n_donors: " << Ndonor << endl;
    
    Naccept_rev_TF1->SetParameter(0, constants.get_gA());
    Naccept_rev_TF1->SetParameter(1, fluence);
    Naccept_rev_TF1->SetParameter(2, constants.get_ka(temperature));
    Naccept_rev_TF1->SetParameter(3, Nacceptor_reversible);
  
    Naccept_rev_TF1_approx->SetParameter(0, constants.get_gA());
    Naccept_rev_TF1_approx->SetParameter(1, fluence);
    Naccept_rev_TF1_approx->SetParameter(2, constants.get_ka(temperature));
    Naccept_rev_TF1_approx->SetParameter(3, Nacceptor_reversible);
  
    Naccept_const_TF1->SetParameter(0, constants.get_gC());
    Naccept_const_TF1->SetParameter(1, fluence);
  
    Nneutrals_rev_TF1->SetParameter(0, constants.get_gY());
    Nneutrals_rev_TF1->SetParameter(1, fluence);
    Nneutrals_rev_TF1->SetParameter(2, constants.get_ky(temperature));
    Nneutrals_rev_TF1->SetParameter(3, Nneutral_reversible);
  
    Ndonor_neutrals_rev_TF1_approx->SetParameter(0, constants.get_gY());
    Ndonor_neutrals_rev_TF1_approx->SetParameter(1, fluence);
    Ndonor_neutrals_rev_TF1_approx->SetParameter(2, constants.get_ky(temperature));
    Ndonor_neutrals_rev_TF1_approx->SetParameter(3, Nneutral_reversible);
  
    Ndonor_const_TF1->SetParameter(0, Ndonor_stable_donorremoval);
    Ndonor_const_TF1->SetParameter(1, constants.get_cc());
    Ndonor_const_TF1->SetParameter(2, fluence);
  
    Ndonor_TF1->SetParameter(0, constants.get_gY());
    Ndonor_TF1->SetParameter(1, fluence);
    Ndonor_TF1->SetParameter(2, constants.get_ky(temperature));
    Ndonor_TF1->SetParameter(3, Nneutral_reversible);
  
  
    if(constants.get_ka(temperature)>a && constants.get_ky(temperature)>a)
    {
      Nacceptor_reversible             =  Naccept_rev_TF1->Eval(duration);
      Nacceptor_stable_constdamage     +=  Naccept_const_TF1->Eval(duration);
      Nneutral_reversible              =  Nneutrals_rev_TF1->Eval(duration);
      Ndonor_stable_donorremoval       += Ndonor_const_TF1->Eval(duration);
      Nacceptor_stable_reverseannealing   +=  Ndonor_TF1->Eval(duration);
    }
    else
    {
      cout << "Potential numerical problem due to ultra low temp and therby caused very small ky and ka values. Using approach in order to perform calculation. In general, no problem!"<<endl;
  
      Nacceptor_reversible           =  Naccept_rev_TF1_approx->Eval(duration);
      Nacceptor_stable_constdamage   +=  Naccept_const_TF1->Eval(duration);
      Nneutral_reversible            =  Ndonor_neutrals_rev_TF1_approx->Eval(duration);
      Ndonor_stable_donorremoval     += Ndonor_const_TF1->Eval(duration);
      Nacceptor_stable_reverseannealing +=  Ndonor_TF1->Eval(duration);
    }

    Nconstdamage_g1   = Nacceptor_stable_constdamage/constants.get_gC();
    Nbenef_anneal_g1  = Nacceptor_reversible/constants.get_gA();
    Ndonor_g1_TF1->SetParameter(0, fluence);
    Ndonor_g1_TF1->SetParameter(1, constants.get_ky(temperature));
    Ndonor_g1_TF1->SetParameter(2, Nnadefects_g1);
    Nrevers_anneal_g1 += Ndonor_g1_TF1->Eval(duration);
    Nnadefects_g1     = Nneutral_reversible/constants.get_gY();
  
    // cout << "ka = " << constants.get_ka(temperature) << "; ky = " << constants.get_ky(temperature) << ": " << NtoV_function(Nacceptor_stable_reverseannealing) << " = " << NtoV_function(Nacceptor_stable_reverseannealing_g) << "+" << NtoV_function(Nacceptor_stable_reverseannealing_no_g) <<  endl;
    // cout << "n_acceptor_reversible: " << Nacceptor_reversible << endl;
    // cout << "Nacceptor_stable_constdamage: " << Nacceptor_stable_constdamage << endl;
    // cout << "Nacceptor_stable_reverseannealing: " << Nacceptor_stable_reverseannealing << endl;
    Nacceptor =  Nacceptor_reversible + Nacceptor_stable_constdamage + Nacceptor_stable_reverseannealing;
    Ndonor    =  Ndonor_const         + Ndonor_stable_donorremoval;
  
  
    //calculating the leackage current in the following part
  
    vector<double> tmp(3);
    vector<double> tmp2(2);
    vector<double> tmp3(2);
  
    tmp3.at(0) = 0;
    tmp3.at(1) = 0;
  
    //////////////////////////////
    // G_i_tmp needs to be outside of the if clause to be existend also
    // afterwards when it is pushed back in the G_i vector, can be put
    // back to pos1 after finished with the if clauses.
    double G_i_tmp = 0;  
  
    if(use_CMS_paper_alpha0)
      {
        double CMS_alpha0_c1=-8.9e-17;
        double CMS_alpha0_c2=4.6e-14;
        double CMS_beta=2.9e-18;
  
        double inverse_weighted_temp_tmp=0;
        inverse_weighted_temp_tmp = duration / (double)temperature;
        inverse_weighted_temp.push_back(inverse_weighted_temp_tmp);
  
        tmp.at(0) = fluence * duration * leakageCurrentConstants.alpha_1;
        // needs to be updated at every itaration in the while loop since the
        // averaged temperature needs to be recalculated
        tmp.at(1) = fluence * duration;   
        tmp.at(2) = fluence * duration * CMS_beta;
  
        leakage_current_alpha.push_back(tmp);

        // time comes in seconds and k01 is seconds-1 so units are fine here
        tmp2.at(0) = duration *leakageCurrentConstants.k01*exp(-leakageCurrentConstants.E1/(8.6173303e-5*(double)temperature));
        tmp2.at(1) = duration / 60.0;   // convert time from seconds to minutes
  
        if(time_history.size() > 0.1) 
        {
            tmp3 = time_history.at(time_history.size()-1);
        }
  
        tmp2.at(0) += tmp3.at(0);
        tmp2.at(1) += tmp3.at(1);
  
        time_history.push_back(tmp2);
  
  
        int i = 0;
        double temperature_average_i = 0;
  
        while (i < leakage_current_alpha.size())
          {
            temperature_average_i = 0;
            // Calculate the average inverse temperature weighted with the
            // individual time from the moment of irradiation (i) to now
            // (leakage_current_alpha.size())
            for(int j=i; j<leakage_current_alpha.size(); j++)                     
              {
                temperature_average_i += inverse_weighted_temp.at(j);
              }
            //for i=0 the time is not the difference but just the total time
            if (i >= 1)
            {
              temperature_average_i /= (double)(time_history.back().at(1)*60.0-time_history.at(i-1).at(1)*60.0);
              G_i_tmp +=   leakage_current_alpha.at(i).at(0) * exp(-(time_history.back().at(0) - time_history.at(i-1).at(0)))
                         + leakage_current_alpha.at(i).at(1) * (CMS_alpha0_c1 + CMS_alpha0_c2*temperature_average_i)
                         - leakage_current_alpha.at(i).at(2) * log(time_history.back().at(1) - time_history.at(i-1).at(1));
            }
            else
            {
              temperature_average_i /= (double)(time_history.back().at(1)*60.0);
              G_i_tmp +=   leakage_current_alpha.at(i).at(0) * exp(-time_history.back().at(0))
                         + leakage_current_alpha.at(i).at(1) * (CMS_alpha0_c1 + CMS_alpha0_c2*temperature_average_i)
                         - leakage_current_alpha.at(i).at(2) * log(time_history.back().at(1));
            }
  
            i++;
          }
      }
  
    else
    {
        tmp.at(0) = fluence * duration * leakageCurrentConstants.alpha_1;
        // Temperature dependence of alpha0 is in the theta function of the
        // time calculation - high temp = longer times and vice versa
        tmp.at(1) = fluence * duration * leakageCurrentConstants.alpha_0_star;
        tmp.at(2) = fluence * duration * leakageCurrentConstants.beta;
  
        leakage_current_alpha.push_back(tmp);

        // time comes in seconds and k01 is seconds-1 so units are fine here
        tmp2.at(0) = duration * leakageCurrentConstants.k01 * exp(-leakageCurrentConstants.E1 / (8.6173303e-5 * (double)temperature));
        // T-ref is hardcoded and fixed, it cannot be changed in a trivial way
        // since the other parameters, especially alpha 0 star were evaluated
        // at this reference temperature!!!
        float Tref = 293.15;
        tmp2.at(1) = duration / 60. * exp(-leakageCurrentConstants.E1_star * (1.0 / (double)temperature - 1.0/Tref) / (8.6173303e-5));
  
        if(time_history.size() > 0.1) 
        {
            tmp3 = time_history.at(time_history.size()-1);
        }
  
        tmp2.at(0) += tmp3.at(0);
        tmp2.at(1) += tmp3.at(1);
        time_history.push_back(tmp2);
  
        int i = 0;
        while (i < leakage_current_alpha.size())
        {
            if (i >= 1)
            {
                G_i_tmp +=   leakage_current_alpha.at(i).at(0) * exp(-(time_history.back().at(0) - time_history.at(i-1).at(0)))
                           + leakage_current_alpha.at(i).at(1)
                           - leakage_current_alpha.at(i).at(2) * log(time_history.back().at(1) - time_history.at(i-1).at(1));
            }
            else
            {
                G_i_tmp +=   leakage_current_alpha.at(i).at(0) * exp(-time_history.back().at(0))
                           + leakage_current_alpha.at(i).at(1)
                           - leakage_current_alpha.at(i).at(2) * log(time_history.back().at(1));
            }
            i++;
        }
    }
  
    G_i.push_back(G_i_tmp);
    
    // TODO: Why these numbers
    double alpha = G_i_tmp / (totalDose + (double)fluence + 0.001);
    if (alpha > 5e-16)
    {
        alpha_vec.push_back(1e-17);
    }
    else
    {
        alpha_vec.push_back(alpha);
    }

    //insert volume here for current per module
    leakage_current.push_back(G_i_tmp * 1000.0 * thickness * 1.e-4);
    powerconsumption.push_back(leakage_current.back() * NtoV_function(get_Neff()));
  
    return;
}


//function to read in the temp/rad profile

double leakageCurrentScale(double leakageCurrent, double T, double Tref, double bandGap){
    // cout << (Tref-T) << ": " << leakageCurrent << " ~ " << leakageCurrent*(T*T/(Tref*Tref)*exp(-(bandGap/(2*8.6173303e-5))*(1/T - 1/Tref))) << endl;
    return leakageCurrent*(T*T/(Tref*Tref)*exp(-(bandGap/(2*8.6173303e-5))*(1/T - 1/Tref)));
}

vector<DataElement> getProfile(string filename)
{
    vector<DataElement> profile;
    DataElement profileSnapshot;

    ifstream input_file(filename.c_str());
    if (input_file.is_open())
    {
        cout << "Reading temperature/radiation file: " << filename << endl;

        string line;
      
        string fillStr;
        string timestampStr;
        string timestepStr;
        string temperatureStr;
        string doseRateStr;
        string leakageCurrentStr;
      
        int timestamp;
        int timestep;
        float temperature;
        long int doseRate;
        float leakageCurrent;
      
        float dleakageCurrent = 0.;

        // Read out the header
        getline (input_file, line);

        bool endOfFile = false;
        while (not endOfFile)
        {
            input_file >> fillStr >> timestampStr >> timestepStr >> temperatureStr >> doseRateStr >> leakageCurrentStr;
            timestamp = stoi(timestampStr);
            timestep = stoi(timestepStr);
            temperature = stof(temperatureStr);
            doseRate = stol(doseRateStr);
            leakageCurrent = stof(leakageCurrentStr);

            // TODO: Temporary hack to fix issue in profile preparation
            // doseRate = doseRate / timestep;

            // if(debug) cout << timestamp << " " << timestep << " " << temperature << " " << doseRate << " " << leakageCurrent << endl;

            profileSnapshot.timestamp = timestamp;
            profileSnapshot.duration = timestep;
            profileSnapshot.doseRate = doseRate * DoseRateScaling;
            profileSnapshot.temperature = temperature;
            profileSnapshot.leakageCurrentData = leakageCurrentScale(leakageCurrent, userTref, profileSnapshot.temperature, bandGap);
            profileSnapshot.dleakageCurrentData = leakageCurrentScale(dleakageCurrent, userTref, profileSnapshot.temperature, bandGap);

            profile.push_back(profileSnapshot);

            // if(debug) cout << "Line read: " << timestamp << timestep << temp1 << cpTemp << rad1 << leakageCurrent << endl;
            //cout << timestamp << endl;
            //cout << timestep << endl;
            //cout << temp1 << endl;
            //cout << cpTemp << endl;
            //cout << rad1 << endl;
            //cout << leakageCurrent << endl;
            //exit(0);

            if(input_file.eof()) {break;}                               // if end of file is reached: brake the loop

            limit_file_size_input--;
            if(limit_file_size_input==0) {break;}
        }
    }
    else
    {
        cout << "Error opening the temperature/radiation file!"<<endl;
        return profile;
    }

    input_file.close();

    return profile;
}

boost::posix_time::ptime getBeginTime(vector<DataElement> profile)
{
    int timestamp = profile.at(0).timestamp;
    boost::posix_time::ptime beginTime(boost::posix_time::from_time_t(time_t(timestamp)));

    return beginTime;
} 


//function to plot different informations in root and pdf files

void plot_vectors(
    vector<double> vecx,
    vector<double> vecy,
    string yName,
    string xName,
    long double totalDose,
    string plotname,
    int mode,
    boost::posix_time::ptime beginTime,
    string rootOutputFileName
)
{
    string rootOutputFileNamePdf=rootOutputFileName;
    rootOutputFileNamePdf.insert(rootOutputFileNamePdf.size(), "_");
    rootOutputFileNamePdf.insert(rootOutputFileNamePdf.size(), plotname);
    rootOutputFileNamePdf.insert(rootOutputFileNamePdf.size(), ".pdf");

    TFile *file = new TFile(rootOutputFileName.c_str(), "UPDATE");
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
            boost::posix_time::ptime t2 = beginTime + boost::posix_time::hours(hours+day*24)
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
    const TVectorD t_timevector(vecx.size(), &vecx[0]);
    const TVectorD t_vdepvector(vecy.size(), &vecy[0]);
    const TVectorD t_time_error(vecx.size());
    const TVectorD t_vdep_error(vecy.size());

    TCanvas * c1 = new TCanvas("TCanvas_name", "Title", 0, 0, 1024, 668);
    gStyle->SetOptTitle(0);

    TGraphErrors *gr = new TGraphErrors(t_timevector, t_vdepvector, t_time_error, t_vdep_error);

    if(mode==1 || mode==2 || mode==66 || mode==8 || mode==888 || mode==999) //in case of plot as a function of time, convert the days from simulation into a more handy date
    {
        gr->GetXaxis()->SetRangeUser(1488326400., 1577836799.);
        gr->GetXaxis()->SetTimeDisplay(1);
        gr->GetXaxis()->SetNdivisions(6, 2, 0);
        gr->GetXaxis()->SetTimeFormat("%d/%m/%Y");
        gr->GetXaxis()->SetTimeOffset(0, "gmt");
    }

    gr->GetXaxis()->SetTitle(xName.c_str());
    gr->GetYaxis()->SetTitle(yName.c_str());

    // if(mode!=66)
    // {
    //   gr->SetName(plotname.c_str());
    //   gr->Draw("AP");
    //   if(plot_pdf) c1->Print(rootOutputFileNamePdf.c_str());
    //   gr->Write();
    // }

    gr->SetName(plotname.c_str());
    gr->Draw("AP");
    if(plot_pdf) c1->Print(rootOutputFileNamePdf.c_str());
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


int main(int argc, char *argv[]) {

    po::variables_map args = parseArguments(argc, argv);
    bool validArgs = checkArguments(args);
    if (!validArgs) { return 1; }

    string profileInputFileName = args["input-profile"].as<string>();
    string rootOutputFileName = args["output-root-file"].as<string>();

    if(!overwrite_files)
    {
        time_t current_time = time(nullptr);
    }

    // TODO: Read constants values from a config file / understand those values
    LeakageCurrentConstants leakageCurrentConstants;
    initializeLeakageCurrentConstants(leakageCurrentConstants);

    AnnealingConstants annealingConstants;
    annealingConstants.initialize("config/annealing_constants.py");

    vector<DataElement> profile = getProfile(profileInputFileName);

    int max_steps = profile.size();
    cout << "Profile succesfully read. Length of the profile is: " << max_steps << endl;

    boost::posix_time::ptime beginTime = getBeginTime(profile);

    vector<double> V_dep_vector;
    vector<double> Neff_vector;
    vector<double> Ndonor_vector;
    vector<double> Nacceptor_vector;
    vector<double> temperature_vector;
    // vector<double> cpTemperature_vector;

    vector<double> time_vector;
    vector<double> fluence_vector;
    vector<double> flux_vector;

    vector<double> time_vector_data;
    vector<double> fluence_vector_data;
    vector<double> leakageCurrentData;
    vector<double> dleakageCurrentData;

    vector<double> N_benef_anneal_g1_vec;
    vector<double> N_revers_anneal_g1_vec;
    vector<double> N_nadefects_g1_vec;
    vector<double> N_constdamage_g1_vec;
    vector<double> N_donor_vec;
    // vector<double> N_donor_stable_donorremoval;

    Sensor* sensor = new Sensor();

    long double time = 0.0;
    long double totalDose = 0.;

    // main loop where all irradiation steps happen

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
    //
    // if (adjust_gAgC){
    //   for(int i=0; i < gA_Nsteps; i++){
    //     gC_float = gC_float_min;
    //     for(int j=0; j < gC_Nsteps; j++){
    //       AnnealingConstants annealingConstants( gA_float,           //atof(argv[2]),       //gA
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
    //           sensor->set_temperature(profile.at(t).temperature);
    //           sensor->irradiate(leakageCurrentConstants,annealingConstants,profile.at(t).doseRate,profile.at(t).duration,totalDose);
    //
    //           temperature_vector.push_back(profile.at(t).temperature);
    //           cpTemperature_vector.push_back(profile.at(t).coolPipeTemp);
    //           Neff_vector.push_back(sensor->get_Neff());
    //           V_dep_vector.push_back(NtoV_function(sensor->get_Neff()));
    //           // Ndonor_vector.push_back(sensor->get_Ndonor());
    //           // Nacceptor_vector.push_back(sensor->get_Nacceptor());
    //           // totalDose+=profile.at(t).doseRate/1.0e3*profile.at(t).duration/1.0e3;  // time (seconds) * dose rate (neq/cm2/s)
    //           totalDose+=double(profile.at(t).doseRate)*double(profile.at(t).duration);
    //           time+=profile.at(t).duration;
    //           // if(profile.at(t).leakageCurrentData > 0.){
    //           //   time_vector_data.push_back(time/(24.0*3600.0));
    //           //   fluence_vector_data.push_back(totalDose);
    //           //   leakageCurrentData.push_back(profile.at(t).leakageCurrentData/1000.);
    //           //   dleakageCurrentData.push_back(profile.at(t).dleakageCurrentData/1000.);
    //           //   flux_vector.push_back(profile.at(t).doseRate);
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
    //       plot_vectors(time_vector,V_dep_vector,"U_{dep} [V]","Date [days]",totalDose,(TString("U_dep")+gA_str+gC_str).Data(),66,rootOutputFileName);
    //       sensor = new Sensor();
    //       V_dep_vector.clear();
    //       Neff_vector.clear();
    //       Ndonor_vector.clear();
    //       Nacceptor_vector.clear();
    //       temperature_vector.clear();
    //       cpTemperature_vector.clear();
    //       time_vector.clear();
    //       fluence_vector.clear();
    //       flux_vector.clear();
    //       time_vector_data.clear();
    //       fluence_vector_data.clear();
    //       leakageCurrentData.clear();
    //       dleakageCurrentData.clear();
    //       totalDose = 0;
    //       time = 0;
    //       gC_float += gC_step;
    //     }
    //     gA_float += gA_step;
    //   }
    // }

    for (int t=0; t<max_steps; t++)   // iterate through the profile and irradiate the sensor at each step
    {
        sensor->set_temperature(profile.at(t).temperature);
        sensor->irradiate(
            leakageCurrentConstants,
            annealingConstants,
            profile.at(t).doseRate,
            profile.at(t).duration,
            totalDose
        );

        temperature_vector.push_back(profile.at(t).temperature);
        // cpTemperature_vector.push_back(profile.at(t).coolPipeTemp);
        Neff_vector.push_back(sensor->get_Neff());
        V_dep_vector.push_back(NtoV_function(sensor->get_Neff()));
        Ndonor_vector.push_back(sensor->get_Ndonor());
        Nacceptor_vector.push_back(sensor->get_Nacceptor());

        N_benef_anneal_g1_vec.push_back(sensor->get_Nbenef_anneal_g1());
        N_revers_anneal_g1_vec.push_back(sensor->get_Nrevers_anneal_g1());
        N_nadefects_g1_vec.push_back(sensor->get_Nnadefects_g1());
        N_constdamage_g1_vec.push_back(sensor->get_Nconstdamage_g1());
        N_donor_vec.push_back(sensor->get_Ndonor());
        // N_donor_stable_donorremoval.push_back(NtoV_function(sensor->get_Ndonor_stable_donorremoval()));
        totalDose += double(profile.at(t).doseRate) * double(profile.at(t).duration);
        time += profile.at(t).duration;
        // TODO: There used to be the following hard-coded excluded times, why?
        //    && (time < 14722058 || (time > 16450058 && time < 45480458) || time > 48504458)
        if (profile.at(t).leakageCurrentData > 0. && profile.at(t).doseRate > 0.)
        {
            // cout << profile.at(t).leakageCurrentData << endl;
            time_vector_data.push_back(time/(24.0*3600.0));
            fluence_vector_data.push_back(totalDose/6.);
            leakageCurrentData.push_back(profile.at(t).leakageCurrentData/1000.);
            dleakageCurrentData.push_back(profile.at(t).dleakageCurrentData/1000.);
            flux_vector.push_back(profile.at(t).doseRate);
        }
        time_vector.push_back(time/(24.0*3600.0));  // time vector in days
        fluence_vector.push_back(totalDose);

        float current_progress = 100 * t / max_steps;
        float last_progress = 100 * (t-1) / max_steps;
        float update_frequency = 5.;
        bool print_update = (fmod(current_progress, update_frequency) - fmod(last_progress, update_frequency) > 0);
        if (print_update)
        {
            cout << static_cast<int>(t * 100 / max_steps) << "% done..." << endl;
        }
    }

    // data output, plotting and visualisation

    cout << "Processing finished, writing data..." << endl;

    // plots as function of time
    plot_vectors(time_vector, Neff_vector, "N_{eff} [1/cm^{3}]", "Date [days]", totalDose, "Neff", 1, beginTime, rootOutputFileName);
    plot_vectors(time_vector, Ndonor_vector, "N_{donor} [1/cm^{3}]", "Date [days]", totalDose, "Ndonors", 2, beginTime, rootOutputFileName);
    plot_vectors(time_vector, Nacceptor_vector, "N_{acceptor} [1/cm^{3}]", "Date [days]", totalDose, "Nacceptors", 2, beginTime, rootOutputFileName);
    plot_vectors(time_vector, V_dep_vector, "U_{dep} [V]", "Date [days]", totalDose, "U_dep", 66, beginTime, rootOutputFileName);

    plot_vectors(time_vector, N_benef_anneal_g1_vec, "N_{dep, benef_anneal_g1} [V]", "Date [days]", totalDose, "N_benef_anneal_g1", 66, beginTime, rootOutputFileName);
    plot_vectors(time_vector, N_revers_anneal_g1_vec, "N_{dep, revers_anneal_g1} [V]", "Date [days]", totalDose, "N_revers_anneal_g1", 66, beginTime, rootOutputFileName);
    plot_vectors(time_vector, N_nadefects_g1_vec, "N_{dep, nadefects_g1} [V]", "Date [days]", totalDose, "N_nadefects_g1", 66, beginTime, rootOutputFileName);
    plot_vectors(time_vector, N_constdamage_g1_vec, "N_{dep, constdamage_g1} [V]", "Date [days]", totalDose, "N_constdamage_g1", 66, beginTime, rootOutputFileName);
    plot_vectors(time_vector, N_donor_vec, "N_{dep, donor} [V]", "Date [days]", totalDose, "N_donor", 66, beginTime, rootOutputFileName);
    // plot_vectors(time_vector, N_donor_stable_donorremoval, "N_{dep, donor_stable_donorremoval} [V]", "Date [days]", totalDose, "N_donor_stable_donorremoval", 66, beginTime, rootOutputFileName);

    plot_vectors(time_vector, sensor->get_alpha_vec(), "#alpha [A/cm]", "Date [days]", totalDose, "alpha", 2, beginTime, rootOutputFileName);
    plot_vectors(time_vector, sensor->get_leakage_current(), "I_{leak} (@21C) [mA/cm^{2}]", "Date [days]", totalDose, "I_leak", 2, beginTime, rootOutputFileName);
    // plot_vectors(time_vector, sensor->get_leakage_current(), "I_{leak} (@ -7C) [mA] per module", "Date [days]", totalDose, "I_leak_per_module", 888, beginTime, rootOutputFileName);
    plot_vectors(time_vector, sensor->get_leakage_current(), "I_{leak} (@ -7C) [mA],  1 ROG", "Date [days]", totalDose, "I_leak_per_module", 888, beginTime, rootOutputFileName);
    plot_vectors(time_vector, sensor->get_leakage_current(), "I_{leak} (@userTref) [mA/cm^{3}]", "Date [days]", totalDose, "I_leak_volume", 8, beginTime, rootOutputFileName);
    plot_vectors(time_vector, sensor->get_G_i(), "G_{i} [A/cm^{3}]", "Date [days]", totalDose, "G_i", 2, beginTime, rootOutputFileName);
    plot_vectors(time_vector, sensor->get_powerconsumption(), "P [mW/cm^{2}]", "Date [days]", totalDose, "powerconsumption", 2, beginTime, rootOutputFileName);
    plot_vectors(time_vector, temperature_vector, "Temperature [K]", "Date [days]", totalDose, "temperature", 2, beginTime, rootOutputFileName);
    // plot_vectors(time_vector, cpTemperature_vector, "Temperature [K]", "Date [days]", totalDose, "cp_Temperature", 2, beginTime, rootOutputFileName);
    plot_vectors(time_vector, fluence_vector, "Fluence [n_{eq}/cm^{2}]", "Date [days]", totalDose, "fluence", 2, beginTime, rootOutputFileName);
    plot_vectors(time_vector_data, flux_vector, "Fluence [n_{eq}/cm^{2}/s]", "Date [days]", totalDose, "flux", 2, beginTime, rootOutputFileName);
    plot_vectors(time_vector_data, leakageCurrentData, "I_{leak} (@ -7C) [mA],  1 ROG", "Date [days]", totalDose, "I_leak_per_module_data", 999, beginTime, rootOutputFileName);
    plot_vectors(time_vector_data, dleakageCurrentData, "dI_{leak} (@ -8.5C) [mA] per module", "Date [days]", totalDose, "dI_leak_per_module_data", 999, beginTime, rootOutputFileName);

    // plots as function of dose
    plot_vectors(fluence_vector, Neff_vector, "N_{eff} [1/cm^{3}]", "Fluence [n_{eq}/cm^{2}]", totalDose, "Neff_vs_fluence", 3, beginTime, rootOutputFileName);
    plot_vectors(fluence_vector, Ndonor_vector, "N_{donor} [1/cm^{3}]", "Fluence [n_{eq}/cm^{2}]", totalDose, "Ndonors_vs_fluence", 4, beginTime, rootOutputFileName);
    plot_vectors(fluence_vector, Nacceptor_vector, "N_{acceptor} [1/cm^{3}]", "Fluence [n_{eq}/cm^{2}]", totalDose, "Nacceptors_vs_fluence", 4, beginTime, rootOutputFileName);
    plot_vectors(fluence_vector, V_dep_vector, "U_{dep} [V]", "Fluence [n_{eq}/cm^{2}]", totalDose, "U_dep_vs_fluence", 4, beginTime, rootOutputFileName);
    plot_vectors(fluence_vector, sensor->get_alpha_vec(), "#alpha [A/cm]", "Fluence [n_{eq}/cm^{2}]", totalDose, "alpha_vs_fluence", 4, beginTime, rootOutputFileName);
    plot_vectors(fluence_vector, sensor->get_leakage_current(), "I_{leak} [mA/cm^{2}]", "Fluence [n_{eq}/cm^{2}]", totalDose, "I_leak_vs_fluence", 4, beginTime, rootOutputFileName);
    plot_vectors(fluence_vector, sensor->get_leakage_current(), "I_{leak} (@21C) [mA] per module", "Date [days]", totalDose, "I_leak_per_module_vs_fluence", 889, beginTime, rootOutputFileName);
    plot_vectors(fluence_vector, sensor->get_leakage_current(), "I_{leak} (@userTref) [mA/cm^{3}]", "Fluence [n_{eq}/cm^{2}]", totalDose, "I_leak_volume_vs_fluence", 9, beginTime, rootOutputFileName);
    plot_vectors(fluence_vector, sensor->get_G_i(), "G_{i} [A/cm^{3}]", "Fluence [n_{eq}/cm^{2}]", totalDose, "G_i_vs_fluence", 4, beginTime, rootOutputFileName);
    plot_vectors(fluence_vector, sensor->get_powerconsumption(), "P [mW/cm^{2}]", "Fluence [n_{eq}/cm^{2}]", totalDose, "powerconsumption_vs_fluence", 4, beginTime, rootOutputFileName);
    plot_vectors(fluence_vector, temperature_vector, "Temperature [K]", "Fluence [n_{eq}/cm^{2}", totalDose, "temperature_vs_fluence", 4, beginTime, rootOutputFileName);
    plot_vectors(fluence_vector_data, leakageCurrentData, "I_{leak} (@ -7C) [mA] per module", "Date [days]", totalDose, "I_leak_per_module_data_vs_fluence", 9999, beginTime, rootOutputFileName);
    plot_vectors(fluence_vector_data, dleakageCurrentData, "dI_{leak} (@ -8.5C) [mA] per module", "Date [days]", totalDose, "dI_leak_per_module_data_vs_fluence", 9999, beginTime, rootOutputFileName);

    cout << rootOutputFileName << " has been created" << endl;

    return 0;
}
