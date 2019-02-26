#ifndef CONFIG_DATA
#define CONFIG_DATA

#include <string>
#include <fstream>

#include "Types.h"


class ConfigData
{
private:
	int NDim;
	int NIons;
	int NElec;
	Box_t BoxSize;
	std::string InitPos; 
	double dT;


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
	void setInitPos(std::string _InitPos) { InitPos = _InitPos; }
	void setTimeStep(double _dT) { dT = _dT; }
	
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
		std::cout << "\t# InitPos = " << InitPos << std::endl;
		std::cout << "\t# dT = " << dT << std::endl;
	}

};




#endif /* CONFIG_DATA */