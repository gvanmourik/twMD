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
	int NBinsX;
	int NBinsY;
	int NBinsZ;
	BoxSize_t BoxSize;
	double Radius;
	std::string InitPos; 
	double dT;
	double Cr;  
	double Ci;
	double S;
	double Rho;
	double Z;		//atomic charge
	double M; 		//atomic mass


public:
	// ADD errors if parameters are not specified!!!!!!!!!!
	ConfigData() {}
	// ConfigData(int _ND, int _NI, int _NE, BoxSize_t _BoxSize, std::string _InitPos, 
	// 		   double _dT, int BoxMag=1) : 
	// 		   		NDim(_ND), NIons(_NI), NElec(_NE), BoxSize(_BoxSize), 
	// 		   		InitPos(_InitPos), dT(_dT) {}
	~ConfigData() {}


	void setNDim(int _NDim) { NDim = _NDim; }
	void setNIons(int _NIons) { NIons = _NIons; }
	void setNElec(int _NElec) { NElec = _NElec; }
	void setNBinsX(int _NBinsX) { NBinsX = _NBinsX; }
	void setNBinsY(int _NBinsY) { NBinsY = _NBinsY; }
	void setNBinsZ(int _NBinsZ) { NBinsZ = _NBinsZ; }
	void setBoxSize(BoxSize_t _BoxSize) { BoxSize = _BoxSize; }
	void setCutoffRadius(double _radius) { Radius = _radius; }
	void setInitPos(std::string _InitPos) { InitPos = _InitPos; }
	void setTimeStep(double _dT) { dT = _dT; }
	void setCr(double _Cr) { Cr = _Cr; }
	void setCi(double _Ci) { Ci = _Ci; }
	void setS(double _S) { S = _S; }
	void setRho(double _Rho) { Rho = _Rho; }
	void setZ(double _Z) { Z = _Z; }
	void setM(double _M) { M = _M; }
	
	int getNDim() const { return NDim; }
	int getNIons() const { return NIons; }
	int getNElec() const { return NElec; }
	int getNBinsX() const { return NBinsX; }
	int getNBinsY() const { return NBinsY; }
	int getNBinsZ() const { return NBinsZ; }
	BoxSize_t getBoxSize() const { return BoxSize; }
	double getCutoffRadius() const { return Radius; }
	std::string getInitPos() const { return InitPos; }
	double getdT() const { return dT; }
	double getCr() const { return Cr; }
	double getCi() const { return Ci; }
	double getS() const { return S; }
	double getRho() const { return Rho; }
	double getCharge() const { return Z; }
	double getMass() const { return M; }


	void print()
	{
		std::cout << std::endl;
		std::cout << "# Configuration Parameters:" << std::endl;
		std::cout << "\t# NDim = " << NDim << std::endl;
		std::cout << "\t# NIons = " << NIons << std::endl;
		std::cout << "\t# NElec = " << NElec << std::endl;
		std::cout << "\t# NBinsX = " << NBinsX << std::endl;
		std::cout << "\t# NBinsY = " << NBinsY << std::endl;
		std::cout << "\t# NBinsZ = " << NBinsZ << std::endl;
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
		std::cout << "# Atomic Parameters:" << std::endl;
		std::cout << "\t# AtomicCharge = " << Z << std::endl;
		std::cout << "\t# AtomicMass = " << M << std::endl;
		std::cout << std::endl;
	}

};




#endif /* CONFIG_DATA_H */