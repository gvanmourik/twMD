#ifndef ELECTRON_H
#define ELECTRON_H

#include "Gaussian.h"

typedef std::vector<Gaussian*> GaussianList_t;

class Electron
{
private:
	GaussianList_t Gaussians;

public:
	// Access functions
	GaussianList_t getGaussians() { return Gaussians; }


	void addGaussian(double Cr, double Ci, double S, double Rho, Position* Pos) 
	{ 
		Gaussian* g = new Gaussian(Cr, Ci, S, Rho, Pos);
		Gaussians.push_back(g);
	}

};



#endif /* ELECTRON_H */