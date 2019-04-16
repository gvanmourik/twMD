#ifndef GAUSSIAN_H
#define GAUSSIAN_H

#include <boost/serialization/vector.hpp>

#include "Types.h"
#include "ParticleInfo.h"

class Electron;
class Gaussian;
typedef std::vector<Gaussian*> GaussianList_t;

class Gaussian
{
private:
	friend class boost::serialization::access;

	double cr; 		//real-part of the gaussian mag
	double ci; 		//imag-part of the gaussian mag
	double s; 		//gaussian width
	double rho;		//gaussian momentum
	Electron* elec;	//pointer to the associated electron
	Position* Pos;	//pointer to the position of the gaussian

	template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & cr;
        ar & ci;
        ar & s;
        ar & rho;
        ar & elec;
        ar & Pos;
    }

public:
	Gaussian() {}
	Gaussian(Electron* _elec, double _Cr, double _Ci, double _S, double _Rho, Position* _Pos) :
		elec(_elec), cr(_Cr), ci(_Ci),
		s(_S), rho(_Rho), Pos(_Pos) {}
	~Gaussian() {}

	// Access functions
	void setCr(double _Cr) { cr = _Cr; }
	void setCi(double _Ci) { ci = _Ci; }
	void setS(double _S) { s = _S; }
	void setRho(double _Rho) { rho = _Rho; }
	void setPos(Position* _Pos) { Pos = _Pos; }

	double Cr() const { return cr; }
	double Ci() const { return ci; }
	double S() const { return s; }
	double Rho() const { return rho; }
	Electron* e() const { return elec; }
	Position* pos() const { return Pos; }


	void print()
	{
		Pos->print();
	}

};


#endif /* GAUSSIAN_H */