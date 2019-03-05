#ifndef CONFIG_DATA_H
#define CONFIG_DATA_H

#include "Types.h"
#include "SourceIncludes.h"


class ConfigData
{
private:
	int NDim;
	int NIons;
	int NElec;
	Box_t BoxSize;
	double Radius;
	std::string InitPos; 
	double dT;
	double Cr;  
	double Ci;
	double S;
	double Rho;
	double Z;		//atom charge
	double M; 		//atom mass


public:
	// DEFAULTS??
	ConfigData() {}
	// ConfigData(int _ND, int _NI, int _NE, Box_t _BoxSize, std::string _InitPos, 
	// 		   double _dT, int BoxMag=1) : 
	// 		   		NDim(_ND), NIons(_NI), NElec(_NE), BoxSize(_BoxSize), 
	// 		   		InitPos(_InitPos), dT(_dT) {}
	~ConfigData() {}


	void setNDim(int _NDim) { NDim = _NDim; }
	void setNIons(int _NIons) { NIons = _NIons; }
	void setNElec(int _NElec) { NElec = _NElec; }
	void setBoxSize(Box_t _BoxSize) { BoxSize = _BoxSize; }
	void setCutoffRadius(double _radius) { Radius = _radius; }
	void setInitPos(std::string _InitPos) { InitPos = _InitPos; }
	void setTimeStep(double _dT) { dT = _dT; }
	void setCr(double _Cr) { Cr = _Cr; }
	void setCi(double _Ci) { Ci = _Ci; }
	void setS(double _S) { S = _S; }
	void setRho(double _Rho) { Rho = _Rho; }
	void setZ(double _Z) { Z = _Z; }
	void setM(double _M) { M = _M; }
	
	int getNDim() { return NDim; }
	

	void print()
	{
		std::cout << "# Configuration Parameters:" << std::endl;
		std::cout << "\t# NDim = " << NDim << std::endl;
		std::cout << "\t# NIons = " << NIons << std::endl;
		std::cout << "\t# NElec = " << NElec << std::endl;
		int count = 0;
		for (auto dim : BoxSize)
		{
			std::cout << "\t# Dim[" << count << "] = " << dim << std::endl;
			count++;
		}
		std::cout << "\t# CutoffRadius = " << NElec << std::endl;
		std::cout << "\t# InitPos = " << InitPos << std::endl;
		std::cout << "\t# dT = " << dT << std::endl;
		std::cout << "# Gaussian Parameters:" << std::endl;
		std::cout << "\t# GaussianMagReal = " << Cr << std::endl;
		std::cout << "\t# GaussianMagImag = " << Ci << std::endl;
		std::cout << "\t# GuassianWidth = " << S << std::endl;
		std::cout << "\t# GaussianMomentum = " << Rho << std::endl;
		std::cout << "# Atom Parameters:" << std::endl;
		std::cout << "\t# AtomCharge = " << Z << std::endl;
		std::cout << "\t# AtomMass = " << M << std::endl;
	}

};




#endif /* CONFIG_DATA_H */