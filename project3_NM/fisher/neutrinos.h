//neutrinos.h
//
// Massive neutrino information and properties

#ifndef NEUTRINOS_H
#define NEUTRINOS_H

#include <string>
#include <vector>

using namespace std;

const double DELTAMSQR12(7.9e-5);    //mass^2 splitting 1-2  eV^2
const double DELTAMSQR13(2.2e-3);    //mass^2 splitting 1-2  eV^2

class Neutrinos {

public:
	Neutrinos();
	~Neutrinos();

//member functions
int getNoNu();
int getNoMassStates();

//initiation routine
void initNeutrinos(double Mnu, int nNu, int nMass, int hflag);
double lightestNeutrinoMass();
double minimumNeutrinoMass();
vector<double> massSplittings();
double omnuhhFromMnu(double mnu);
double mnuFromOmnuhh(double omnuhh);

protected:

	int hierarchy_flag;
	int number_neutrinos;
	int number_mass_states;
	double total_neutrino_mass;

	vector<double> nu_mass_values;
};

#endif