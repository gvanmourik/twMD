#ifndef GAUSSIAN_H
#define GAUSSIAN_H

#include "Types.h"
#include "ParticleInfo.h"

class Gaussian
{
private:
	double Cr; 		//real-part of the gaussian mag
	double Ci; 		//imag-part of the gaussian mag
	double S; 		//gaussian width
	double Rho;		//gaussian momentum
	Position* Pos;	//pointer to the position of the gaussian

public:
	Gaussian() : Cr(NOT_SET_DOUBLE), Ci(NOT_SET_DOUBLE),
				 S(NOT_SET_DOUBLE), Rho(NOT_SET_DOUBLE),
				 Pos(nullptr) {}
	Gaussian(double _Cr, double _Ci, double _S, double _Rho, Position* _Pos) :
				 Cr(_Cr), Ci(_Ci),
				 S(_S), Rho(_Rho),
				 Pos(_Pos) {}
	~Gaussian() {}

	// Access functions
	void setCr(double _Cr) { Cr = _Cr; }
	void setCi(double _Ci) { Ci = _Ci; }
	void setS(double _S) { S = _S; }
	void setRho(double _Rho) { Rho = _Rho; }
	void setPos(Position* _Pos) { Pos = _Pos; }

	bool getCr(double &retCr) 
	{ 
		if (Cr != NOT_SET_DOUBLE) 
		{
			retCr = Cr;
			return true;
		}
		return false;
	}

	bool getCi(double &retCi) 
	{ 
		if (Ci != NOT_SET_DOUBLE) 
		{
			retCi = Ci;
			return true;
		}
		return false;
	}

	bool getS(double &retS) 
	{ 
		if (S != NOT_SET_DOUBLE) 
		{
			retS = S;
			return true;
		}
		return false;
	}

	bool getRho(double &retRho) 
	{ 
		if (Rho != NOT_SET_DOUBLE) 
		{
			retRho = Rho;
			return true;
		}
		return false;
	}

	// pass the pos pointer by reference
	bool getPos(Position* &retPos) 
	{ 
		if (Pos != NOT_SET_POS) 
		{
			retPos = Pos;
			return true;
		}
		return false;
	}

};


#endif /* GAUSSIAN_H */