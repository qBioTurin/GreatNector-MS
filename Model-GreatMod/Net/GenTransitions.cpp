#include <cstdio>
#include <cstring>
#include <cmath>



// general transitions:

double TeffDup(double *Value,  map <string,int>& NumTrans,  map <string,int>& NumPlaces,
               const vector<string> & NameTrans, const struct InfTr* Trans, const int T,
               const double& time, double TeDup, double TeDup2)
{
  //cout << "Trans: " <<  NameTrans[T] << endl;

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

double TregDup(double *Value,  map <string,int>& NumTrans,  map <string,int>& NumPlaces, const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time, double TrDup, double TrDup2)
{


	double rate = 0.0;
	double Teff = 0.0 ;
	double Treg = 0.0 ;

	Treg = Value[NumPlaces.find("Treg_out") -> second] ;
	Teff = Value[NumPlaces.find("Teff_out_t1") -> second] ;
	
	rate = Treg  * (TrDup + TrDup2 * Teff);

	return rate ;

}

double Killing(double *Value, map <string,int>& NumTrans, map <string,int>& NumPlaces,const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time, double rate)
{

	double intensity = 1.0;
	
	//rate = RATES_killingTransitions.find(NameTrans[T]) -> second ;

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

double ATZKill(double *Value, map <string,int>& NumTrans, map <string,int>& NumPlaces,const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time, double rate)
{


	double intensity = 1.0;
	//double rate = 0.0;
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

	rate = rate / cells ;

	return  rate * intensity;

	// if(frate <= 0.0 ) frate = 0.0;

	// cout << "Time: " << time << " ; rate: " << frate  << "; transition: " << NameTrans[T] <<"\n"<< endl;
	// cout << "Tcell: " << Tcell  << "; intensity: " << intensity << " ; rate: " << rate  << "\n"<< endl;
	// return frate;
}

double AntigEntry(double *Value, map <string,int>& NumTrans, map <string,int>& NumPlaces,const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time, double A, double B, double Period)
{


  double rate = A * sin(2*M_PI* time / Period) + B;
  return  rate;
  
}
