// Contains functions for calculating fisher matrix


#include <math.h>
#include <iostream>
#include <fstream>
#include "neutrinos.h"
#include "dcosmology2.h"
#include "dnumrecipes.h"
#include <sstream>

using namespace std;

//References:
//Spitzer:  Spitzer, "Physical Processes in the ISM".


/*********************************************************************
 **************** Constructor/Destructor *****************************
 ********************************************************************/

 Neutrinos::Neutrinos()
{
	//cout<<"Fisher constructor called"<<endl;
}

 Neutrinos::~Neutrinos()
{
	//cout<<"Fisher destructor called"<<endl;

}


///////////////////////////////////////////////////////////////////////
// Member functions
//////////////////////////////////////////////////////////////////////
int Neutrinos::getNoNu()
{
	return number_neutrinos;
}

int Neutrinos::getNoMassStates()
{
	return number_mass_states;
}

///////////////////////////////////////////////////////////////////////
// Initiation routine
//////////////////////////////////////////////////////////////////////
//
// hierarchy_flag: 0 degenerate, 1 normal, 2 inverted 3 one massive
//
void Neutrinos::initNeutrinos(double Mnu, int nNu, int nMass, int hflag)
{
	total_neutrino_mass=Mnu;
	number_neutrinos=3;
	number_mass_states=3;

	hierarchy_flag=hflag;

	if(total_neutrino_mass<minimumNeutrinoMass()){
		total_neutrino_mass=minimumNeutrinoMass();
		cout<<"Error hierarchy forbids total mass <"<<total_neutrino_mass<<endl;
	}
}

///////////////////////////////////////////////////////////////////////
// Mass splittings
//////////////////////////////////////////////////////////////////////
double Neutrinos::lightestNeutrinoMass()
{
	double lightM;

	if(total_neutrino_mass<minimumNeutrinoMass()){
	cout<<"Error Mnu<minimumMass lightest neutrino mass<0"<<endl;
	 return 0.0;
	}

	if(hierarchy_flag==0){
		lightM=total_neutrino_mass/(double)(number_neutrinos);
	}else if(hierarchy_flag==1){
		//normal
	 	lightM=total_neutrino_mass;
		lightM-=sqrt(DELTAMSQR12)+sqrt(DELTAMSQR13);
		lightM/=(double)(number_neutrinos);
	}else if(hierarchy_flag==2){
		//inverted
	 	lightM=total_neutrino_mass;
		lightM-=sqrt(DELTAMSQR12)+2.0*sqrt(DELTAMSQR13);
		lightM/=(double)(number_neutrinos);
	}else{
		lightM=total_neutrino_mass;
	}

	return lightM;
}

double Neutrinos::minimumNeutrinoMass()
{
	double minM;

	if(hierarchy_flag==0){
		minM=0.0;
	}else if(hierarchy_flag==1){
		//normal
		minM=sqrt(DELTAMSQR12)+sqrt(DELTAMSQR13);
	}else if(hierarchy_flag==2){
		//inverted
		minM=sqrt(DELTAMSQR12)+2.0*sqrt(DELTAMSQR13);
	}else{
		minM=0.0;
	}

	return minM;
}

vector<double> Neutrinos::massSplittings()
{
	int i;
	double lightM;
	double frac;
	vector<double> splittings;

	lightM=lightestNeutrinoMass();
	if(hierarchy_flag==0){
		frac=1.0/(double)(number_neutrinos);
		for(i=1;i<=number_neutrinos;i++){
			splittings.push_back(frac);
		}
	}else if(hierarchy_flag==1){
		frac=lightM/total_neutrino_mass;
		splittings.push_back(frac);
		frac=(lightM+sqrt(DELTAMSQR12))/total_neutrino_mass;
		splittings.push_back(frac);
		frac=(lightM+sqrt(DELTAMSQR13))/total_neutrino_mass;
		splittings.push_back(frac);
	}else if(hierarchy_flag==2){
		frac=lightM/total_neutrino_mass;
		splittings.push_back(frac);
		frac=(lightM+sqrt(DELTAMSQR13))/total_neutrino_mass;
		splittings.push_back(frac);
		frac=(lightM+sqrt(DELTAMSQR12)+sqrt(DELTAMSQR13))/total_neutrino_mass;
		splittings.push_back(frac);
	}else{
		splittings.push_back(1.0);
	}

	return splittings;
}

double Neutrinos::omnuhhFromMnu(double mnu)
{
	return mnu/93.14;
}

double Neutrinos::mnuFromOmnuhh(double omnuhh)
{
	return omnuhh*93.14;
}