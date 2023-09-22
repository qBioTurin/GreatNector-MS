#include <cstdio>
#include <cstring>
#include <cmath>

// Rate
static double TeDup;
static double TeDup2;
static double TrDup;
static double TrDup2;
static double TeffkillsA;
static double TregKillsTeff;
static double ATZkillrate;
static double A;
static double B;
static double Period;
static map <string, double> RATES_killingTransitions;
static bool populate_data_structures = true;

// Definition of new structures:

void read_double(string fname, double& d)
{
	//cout << "#### " << fname << "####" << endl;
	ifstream f (fname);
	string line;
	if(f.is_open())
	{
		getline(f,line);
		d = stod(line);
		f.close();
		cout << d << endl;
	}
	else
	{
		//std::cerr<<"\nUnable to open " << fname << ": file do not exists\": file do not exists\n";
		exit(EXIT_FAILURE);
	}
}

void init_data_structures( map <string,int>& NumPlaces)
{

	read_double("./TeDup",TeDup);
	read_double("./TeDup2",TeDup2);
	read_double("./TrDup",TrDup);
	read_double("./TrDup2",TrDup2);
	read_double("./TrkTe",TregKillsTeff);
	read_double("./TekA",TeffkillsA);
	read_double("./ATZkill",ATZkillrate);

	read_double("./A",A);
	read_double("./B",B);
	read_double("./Period",Period);
	
	RATES_killingTransitions={{"TeffkillsA_teff_t17", TeffkillsA },
                           {"TeffkillsA_teff_t1", TeffkillsA },
                           {"TregKillsTeff_out_teff_t1", TregKillsTeff},
                           {"TregKillsTeff_out_teff_t17", TregKillsTeff}};
	populate_data_structures = false;
}

// general transitions:

double TeffDup(double *Value,  map <string,int>& NumTrans,  map <string,int>& NumPlaces, const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time)
{
  //cout << "Trans: " <<  NameTrans[T] << endl;
	if( populate_data_structures )
		init_data_structures(NumPlaces);

	double rate = 0.0;
	double A = 0.0;
	double Teff = 0.0 ;

	A = Value[NumPlaces.find("Antigen") -> second];
	Teff = Value[Trans[T].InPlaces[0].Id]; // It has 1 input place 

	rate = Teff  * (TeDup + TeDup2 * A);

	//cout << "Time: " << time << " ; rate: " << rate << "; transition: " << NameTrans[T] << " index " << idxTeff << "\n";
	//cout << "Teff: " << Teff << endl;

	return  rate ; // I have to split since we have two transitions for each color

}

double TregDup(double *Value,  map <string,int>& NumTrans,  map <string,int>& NumPlaces, const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time)
{

	//cout << "Trans: " <<  NameTrans[T] << endl;
	if( populate_data_structures )
		init_data_structures(NumPlaces);

	double rate = 0.0;
	double Teff = 0.0 ;
	double Treg = 0.0 ;

	Treg = Value[NumPlaces.find("Treg_out") -> second] ;
	Teff = Value[NumPlaces.find("Teff_out_t1") -> second] ;
	
	rate = Treg  * (TrDup + TrDup2 * Teff);

	return rate ;

}

double Killing(double *Value, map <string,int>& NumTrans, map <string,int>& NumPlaces,const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time)
{

	//cout << "Trans: " << NameTrans[T] << endl;
	if( populate_data_structures )
		init_data_structures(NumPlaces);

	double rate=0.0;
	double intensity = 1.0;
	
	rate = RATES_killingTransitions.find(NameTrans[T]) -> second ;

	double denom = Value[NumPlaces.find("Antigen") -> second] + Value[NumPlaces.find("Teff_out_t1") -> second] + Value[NumPlaces.find("Teff_out_t17") -> second] + Value[NumPlaces.find("Treg_out") -> second];
	
	for (unsigned int k=0; k<Trans[T].InPlaces.size(); k++)
	{
	  intensity *= pow(Value[Trans[T].InPlaces[k].Id],Trans[T].InPlaces[k].Card);
	}
	
	double frate = rate * intensity;

	if(denom > 0.0)
	  frate = frate/denom;
	
	// cout << "Time: " << time << " ; rate: " << frate << "; transition: " << NameTrans[T] <<"\n"<< endl;
	return frate;
}

double ATZKill(double *Value, map <string,int>& NumTrans, map <string,int>& NumPlaces,const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time)
{

	//cout << "Trans: " <<  NameTrans[T] << endl;
	if( populate_data_structures )
		init_data_structures(NumPlaces);

	double intensity = 1.0;
	double rate = 0.0;
	double Teff1 = Value[NumPlaces.find("Teff_out_t1") -> second];
	double Teff17 = Value[NumPlaces.find("Teff_out_t17") -> second];
	double Treg = Value[NumPlaces.find("Treg_out") -> second];
	double ATZ = Value[NumPlaces.find("ATZ") -> second];

	for (unsigned int k=0; k<Trans[T].InPlaces.size(); k++)
	{
		intensity *= pow(Value[Trans[T].InPlaces[k].Id],Trans[T].InPlaces[k].Card);
	}

	double cells = (Teff1+Teff17+Treg+ATZ);
	if( cells < 1 ) cells = 1;

	rate = ATZkillrate / cells ;

	return  rate * intensity;

	// if(frate <= 0.0 ) frate = 0.0;

	// cout << "Time: " << time << " ; rate: " << frate  << "; transition: " << NameTrans[T] <<"\n"<< endl;
	// cout << "Tcell: " << Tcell  << "; intensity: " << intensity << " ; rate: " << rate  << "\n"<< endl;
	// return frate;
}

double AntigEntry(double *Value, map <string,int>& NumTrans, map <string,int>& NumPlaces,const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time)
{
  //cout << "Trans: " <<  NameTrans[T] << endl;
  if( populate_data_structures )
    init_data_structures(NumPlaces);

  double rate = A * sin(2*M_PI* time / Period) + B;
  return  rate;
  
}
