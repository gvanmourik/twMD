#ifndef ELECTRONS_H
#define ELECTRONS_H

#include "Box.h"
#include "Electron.h"

typedef std::vector<Electron*> ElectronList_t;

class Electrons
{
private:
	ElectronList_t Electrons;

public:
	// Access functions
	ElectronList_t getElectrons() { return Electrons; }


	void addElectron(double Cr, double Ci, double S, double Rho, AtomList_t &Atoms) 
	{ 
		Electron* e;

		//add gaussian to each atom in the box
		for (auto atom : Atoms)
		{
			e->addGaussian(Cr, Ci, S, Rho, atom->pos());
		}

		//add to the electron list
		Electrons.push_back(e);
	}


};



#endif /* ELECTRONS_H */