#ifndef ELECTRON_H
#define ELECTRON_H

#include "Gaussian.h"

class Electron;
typedef std::vector<Electron*> ElectronList_t;

class Electron
{
private:
	friend class boost::serialization::access;
	GaussianList_t Gaussians;

	template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & Gaussians;
    }

public:
	// Access functions
	GaussianList_t getGaussians() { return Gaussians; }
	

	bool addGaussian(Gaussian* g) 
	{ 
		Gaussians.push_back(g);
		return true;
	}


	void print()
	{
		std::cout << "\tGaussians: " << std::endl;
		for (auto g : Gaussians)
		{
			std::cout << "\t\t";
			g->print();
			std::cout << std::endl;
		}
	}

};



#endif /* ELECTRON_H */