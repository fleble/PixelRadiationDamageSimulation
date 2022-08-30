/*
Compile Code with: g++ profile_combiner.cpp -Wall -std=c++11 -o converter

Program reads two files (see pattern below) in and orders them and prints out a file with at least 'time_spacing' hours spacing between the data points in the format time temperature(in K) delta_radiation(in neq/cm2/sec)

temperature file:
21-11-2011 15:45:40:838 +01.99
21-11-2011 19:45:40:838 +02.99

radiation file:
2011-11-21 15:45:40.538 00001.234
2011-11-21 22:48:52.330 00010.142

with luminosity in pb-1

both files have to start at a date later than 1.1.2011 01:01:01.001 and should start at the same date!!!

what needs to be changed before starting:
    1. change conversion factor (from pb-1 to neq/cm2) before usage (two conversion factors are supported for 7/8 and 13 TeV)
    2. change output name
    3. call program after compiling with ./converter temperature_file_name irradiation_file_name
*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>			// Basic IO
#include <fstream>			// Read and write files
#include <cstdlib>			// convert arguments of function call into int, float,...
#include <string>
#include <sstream>         		// to get string into stream
#include <vector>
#include <cmath>

using namespace std;

///////////////////////////////////////////////////////

string define_output_file_name = "profile_output.txt";
long int conversion_factor1   = 0.582e9;                    //conversion factor from pb-1 to neq/cm2 for 7TeV
long int conversion_factor2   = 0.839e9;                    //for 13TeV
long int changeinenergy       = 2014;                       //year< since which there is 13TeV instead of 7TeV
int time_spacing              = 4;                          //minimal required spacing between two data points in hours
int replacement88             = 2200;                       //replace no temperature reading from DCS (eg during shutdown) (=-88 C) with this temperature

bool debug            = false;                             //switch debugging information on/off

///////////////////////////////////////////////////////

struct newData
{
  unsigned long long int timestamp;
  int temperature;
  int luminosity;
};

///////////////////////////////////////////////////////

struct data
{
  int year;
  int month;
  int day;
  int hour;
  int minute;
  int second;
  int milisecond;
  unsigned long long int timestamp;
  int payload;
};

///////////////////////////////////////////////////////

unsigned long long int diffintime(int time1, int time2, int offset, bool uberold, unsigned long long int conversion)
{
  long long int diff=0;
  bool uber = false;

  diff = time1 - time2;

  if(diff<0)
  {
    diff+=offset;
    uber = true;
  }

  unsigned long long diff_long=abs(diff);

  if(uberold)diff_long-=1;

  diff_long*=conversion;

  return diff_long;
}

///////////////////////////////////////////////////////

unsigned long long int BuildTimeStamp(data newpoint, data oldpoint)
{
  bool next_uber=false;
  long long int difference_milisecond = diffintime(newpoint.milisecond, oldpoint.milisecond, 1000, false, 1);
  if(newpoint.milisecond-oldpoint.milisecond<0) next_uber = true;

  unsigned long long int difference_second     = diffintime(newpoint.second, oldpoint.second, 60, next_uber, 1000);
  next_uber=false;
  if(newpoint.second-oldpoint.second<0) next_uber = true;

  unsigned long long int difference_minute     = diffintime(newpoint.minute, oldpoint.minute, 60, next_uber, 60*1000);
  next_uber=false;
  if(newpoint.minute-oldpoint.minute<0) next_uber = true;

  unsigned long long int difference_hour       = diffintime(newpoint.hour, oldpoint.hour, 24, next_uber, 60*60*1000);
  next_uber=false;
  if(newpoint.hour-oldpoint.hour<0) next_uber = true;

  unsigned long long int dayoffset=0;
  if(oldpoint.month == 0) dayoffset=0;
  else if(oldpoint.month == 2 )
    {
      dayoffset=28;
      if(oldpoint.year==2008 || oldpoint.year==2012 || oldpoint.year==2016 || oldpoint.year==2020) dayoffset=29;
    }
  else if(oldpoint.month == 4 || oldpoint.month == 6 || oldpoint.month == 9 || oldpoint.month == 11)dayoffset=30;
  else if(oldpoint.month == 1 || oldpoint.month == 3 || oldpoint.month == 5 || oldpoint.month == 7 || oldpoint.month == 8 || oldpoint.month == 10 || oldpoint.month == 12)dayoffset=31;
  else
    {
      cout<<"Strange month, aborting..."<<oldpoint.month<<endl;
      return -1;
    }


  unsigned long long int difference_day        = diffintime(newpoint.day, oldpoint.day, dayoffset, next_uber, 24*60*60*1000);
    //cout<<"difference in days is: "<<difference_day<< " "<<next_uber<<endl;
  next_uber=false;
  if(newpoint.day-oldpoint.day<0) next_uber = true;
  //cout<<difference_day<< "ddd"<<endl;

  unsigned long long int difference_month      = diffintime(newpoint.month, oldpoint.month, 12, next_uber, dayoffset*24*60*60*1000);
  next_uber=false;
  if(newpoint.month-oldpoint.month<0) next_uber = true;

  unsigned long long int difference_year       = diffintime(newpoint.year, oldpoint.year, 100000, next_uber, 365.25*24*60*60*1000);


  unsigned long long int newtimestamp = oldpoint.timestamp + difference_milisecond + difference_second + difference_minute + difference_hour + difference_day + difference_month + difference_year;

  if(debug)cout<<"Timestamp calculation: "<<oldpoint.timestamp<<"  "<<newtimestamp<<"  "<<difference_milisecond<<"  "<<difference_second<<"  "<<difference_minute<<"  "<<difference_hour<<"  "<<difference_day<<"  "<<difference_month<<"  "<<difference_year<<endl;

  return newtimestamp;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {

      string temp_file_name = argv[1];
      string irr_file_name = argv[2];

      ofstream output_file;
      output_file.open(define_output_file_name.c_str());
      ifstream temp_file(temp_file_name.c_str());
      ifstream irr_file(irr_file_name.c_str());


      if (!temp_file.is_open())
      {
        cout<<"Temp file not found"<<endl;
        return -1;
      }
      if (!irr_file.is_open())
      {
        cout<<"Irr file not found"<<endl;
        return -1;
      }


      cout<<"Temperature File: " << temp_file_name << endl;
      cout<<"Irradiation File: " << irr_file_name << endl;


      char column_dump_temp[100];								// max 100 charachters per line!!!
      char column_dump_irr[100];								// max 100 charachters per line!!!

      vector<data> temp_data_vec;
      vector<data> irr_data_vec;

      data temp_data_element;
      data irr_data_element;

      int year_irr, year_temp;
      int month_irr, month_temp;
      int day_irr, day_temp;
      int hour_irr, hour_temp;
      int minute_irr, minute_temp;
      int second_irr, second_temp;
      int miliseconds_irr, miliseconds_temp;
      unsigned long long int timestamp_irr=0;
      unsigned long long int timestamp_temp=0;
      int irr_i, temp_i;

      irr_data_element={2011,1,1,1,1,1,1,0,0};
      temp_data_element={2011,1,1,1,1,1,1,0,0};

      irr_data_vec.push_back(irr_data_element);
      temp_data_vec.push_back(temp_data_element);


//read in the two files (the two while loops) and put them into a vector each, calculate a unique timestamp for each point (in ms)

      cout<<"Start to read in irr file..."<<endl;

      int counter=0;

      while(true)
      {
        irr_file.getline(column_dump_irr,100);

        if(debug) cout << "Originally read line is: " << column_dump_irr << endl;


        year_irr          =   (static_cast<int>(column_dump_irr[0])-48)*1000+(static_cast<int>(column_dump_irr[1])-48)*100+(static_cast<int>(column_dump_irr[2])-48)*10+static_cast<int>(column_dump_irr[3])-48;
        month_irr         =   (static_cast<int>(column_dump_irr[5])-48)*10+static_cast<int>(column_dump_irr[6])-48;
        day_irr           =   (static_cast<int>(column_dump_irr[8])-48)*10+static_cast<int>(column_dump_irr[9])-48;
        hour_irr          =   (static_cast<int>(column_dump_irr[11])-48)*10+static_cast<int>(column_dump_irr[12])-48;
        minute_irr        =   (static_cast<int>(column_dump_irr[14])-48)*10+static_cast<int>(column_dump_irr[15])-48;
        second_irr        =   (static_cast<int>(column_dump_irr[17])-48)*10+static_cast<int>(column_dump_irr[18])-48;
        miliseconds_irr   =   (static_cast<int>(column_dump_irr[20])-48)*100+(static_cast<int>(column_dump_irr[21])-48)*10+(static_cast<int>(column_dump_irr[22])-48);
        irr_i             =   (static_cast<int>(column_dump_irr[24])-48)*10000000+(static_cast<int>(column_dump_irr[25])-48)*1000000+(static_cast<int>(column_dump_irr[26])-48)*100000+(static_cast<int>(column_dump_irr[27])-48)*10000+(static_cast<int>(column_dump_irr[28])-48)*1000+(static_cast<int>(column_dump_irr[30])-48)*100+(static_cast<int>(column_dump_irr[31])-48)*10+(static_cast<int>(column_dump_irr[32])-48);

        if(debug) cout<<"Datapoint irr read: "<<irr_i<<"  "<<miliseconds_irr<<"  "<<second_irr<<"  "<<minute_irr<<"  "<<hour_irr<<"  "<<day_irr<<"  "<<month_irr<<"  "<<year_irr<<endl;

        if(year_irr>changeinenergy)irr_i=(int)round((double)irr_i*(double)conversion_factor2/(double)conversion_factor1);

        irr_data_element={year_irr,month_irr,day_irr,hour_irr,minute_irr,second_irr,miliseconds_irr,timestamp_irr,irr_i};
        timestamp_irr     = BuildTimeStamp(irr_data_element, irr_data_vec.back());
        irr_data_element={year_irr,month_irr,day_irr,hour_irr,minute_irr,second_irr,miliseconds_irr,timestamp_irr,irr_i};

        if(timestamp_irr == -1) return -1;

        irr_data_vec.push_back(irr_data_element);

        if(irr_file.eof()) break;

        counter++;
        if(counter%200 == 0) cout << counter << " lines read. Year: " << irr_data_element.year << " Month:" << irr_data_element.month << " Day:" << irr_data_element.day << endl;
        if(debug && (counter == 1 || counter ==2))system("read");
      }

      cout<<"Irr file completed! Starting with temperature file..."<<endl;

      counter=0;
      system("read");

      int numberofreplacements=0;

      while(true)
      {
        temp_file.getline(column_dump_temp,100);

        if(debug) cout << "Originally read line is: " << column_dump_temp << endl;


        year_temp = (static_cast<int>(column_dump_temp[6])-48)*1000+(static_cast<int>(column_dump_temp[7])-48)*100+(static_cast<int>(column_dump_temp[8])-48)*10+static_cast<int>(column_dump_temp[9])-48;
        month_temp        =   (static_cast<int>(column_dump_temp[3])-48)*10+static_cast<int>(column_dump_temp[4])-48;
        day_temp          =   (static_cast<int>(column_dump_temp[0])-48)*10+static_cast<int>(column_dump_temp[1])-48;
        hour_temp         =   (static_cast<int>(column_dump_temp[11])-48)*10+static_cast<int>(column_dump_temp[12])-48;
        minute_temp       =   (static_cast<int>(column_dump_temp[14])-48)*10+static_cast<int>(column_dump_temp[15])-48;
        second_temp       =   (static_cast<int>(column_dump_temp[17])-48)*10+static_cast<int>(column_dump_temp[18])-48;
        miliseconds_temp  =   (static_cast<int>(column_dump_temp[20])-48)*100+(static_cast<int>(column_dump_temp[21])-48)*10+(static_cast<int>(column_dump_irr[22])-48);
        temp_i            =(static_cast<int>(column_dump_temp[25])-48)*1000+(static_cast<int>(column_dump_temp[26])-48)*100+(static_cast<int>(column_dump_temp[28])-48)*10+(static_cast<int>(column_dump_temp[29])-48);

        if(column_dump_temp[24]=='-')
        {
          temp_i*=-1;
          if(debug)cout<<"negative temp"<<temp_i<<endl;
        }

        if(temp_i+8800<1)
        {
          numberofreplacements++;
          temp_i=replacement88;
          if(debug)cout << "temperature replaced! " << temp_i<< " "<< numberofreplacements<<"/"<<counter<<endl;
        }
        if(debug) cout<<"Datapoint temp read: "<<temp_i<<"  "<<miliseconds_temp<<"  "<<second_temp<<"  "<<minute_temp<<"  "<<hour_temp<<"  "<<day_temp<<"  "<<month_temp<<"  "<<year_temp<<endl;

        temp_data_element={year_temp,month_temp,day_temp,hour_temp,minute_temp,second_temp,miliseconds_temp,timestamp_temp,temp_i};
        timestamp_temp    = BuildTimeStamp(temp_data_element, temp_data_vec.back());
        temp_data_element={year_temp,month_temp,day_temp,hour_temp,minute_temp,second_temp,miliseconds_temp,timestamp_temp,temp_i};

        if(timestamp_temp == -1) return -1;

        temp_data_vec.push_back(temp_data_element);

        if(temp_file.eof()) break;
        counter++;
        if(counter%1000 == 0) cout << counter << " lines read. Year: " << temp_data_element.year << " Month:" << temp_data_element.month << " Day:" << temp_data_element.day << endl;
      }

      cout<<"Both files succesfully read! Beginning with sorting the data..."<<endl;

//create a new vector with the ordered values in the format: timestamp Temperature luminosity

      irr_data_vec.erase(irr_data_vec.begin());                                 // erase first element each because this was only a space holder
      temp_data_vec.erase(temp_data_vec.begin());

      int temp_iterator = 1;
      int irr_iterator  = 1;

      newData element;
      element={0,0,0};

      vector<newData> output_vec;

      while(true)
      {
        if(irr_data_vec.at(irr_iterator).timestamp>temp_data_vec.at(temp_iterator).timestamp)
        {
          element = {temp_data_vec.at(temp_iterator).timestamp, temp_data_vec.at(temp_iterator).payload, irr_data_vec.at(irr_iterator-1).payload};
          temp_iterator++;
        }

        else
        {
          element = {irr_data_vec.at(irr_iterator).timestamp, temp_data_vec.at(temp_iterator-1).payload, irr_data_vec.at(irr_iterator).payload};
          irr_iterator++;
        }

        output_vec.push_back(element);

        if(temp_data_vec.size()<=temp_iterator || irr_data_vec.size()<=irr_iterator) break;
      }


      cout<<"Data is sorted now! Writing file with the data with fixed time spacing..."<<endl;


      if(debug)
      {
        for(int k=0;k<output_vec.size();k++)
        {
          cout<<output_vec.at(k).timestamp<<" "<<output_vec.at(k).temperature<<" "<<output_vec.at(k).luminosity<<endl;
        }
      }


//write out a file with a constant spacing of time_spacing hours, temperature and the converted irradiation is rounded to int

      output_file<<1*60*60<<" "<<(int)round(273.0+(output_vec.at(0).temperature/100.0))<<" "<<(long int)round((conversion_factor1*(output_vec.at(0).luminosity)/1000.0/60.0/60.0))<<endl;

      int old_i = 0;
      int i     = 1;
      int time_difference=0;
      int radiation_difference=0;

      while(i<output_vec.size()-2)
      {
        while(i<output_vec.size()-1 && output_vec.at(i).timestamp-output_vec.at(old_i).timestamp<time_spacing*60*60*1000)               //wait until difference between to points is at least time_spacing hour
        {
          i++;

          if(debug) cout<<"Position in output file:"<<i << " and file size is: " << output_vec.size() <<endl;
        }

        time_difference=(int)round((output_vec.at(i).timestamp-output_vec.at(old_i).timestamp)/1000.0);
        radiation_difference=(long int)round((conversion_factor1*(output_vec.at(i).luminosity-output_vec.at(old_i).luminosity)/1000.0/(float)time_difference));

        if(output_vec.at(i).luminosity-(int)round((double)conversion_factor2/(double)conversion_factor1*output_vec.at(old_i).luminosity)==0)
        {
          radiation_difference=0;
          cout << "corrected for change in conversion factor, should not happen more than once!"<<endl;
        }
        output_file<<time_difference<<" "<<(int)round(273.0+(output_vec.at(i).temperature/100.0))<<" "<<radiation_difference<<endl;

        old_i=i;
        i++;
      }

      output_file.close();

      return 0;
}
