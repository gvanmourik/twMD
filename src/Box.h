#ifndef BOX_H
#define BOX_H

#include <map>
#include <cmath>
#include <algorithm>
#include <sys/time.h>
#include <boost/random.hpp>
#include <boost/functional/hash.hpp>

#include "Atom.h"
#include "Types.h"
#include "BinPos.h"
#include "ConfigData.h"
#include "ParticleInfo.h"
#include "SourceIncludes.h"

enum dim{X,Y,Z};

typedef std::vector<Atom*> AtomList_t;
typedef std::map<std::size_t, AtomList_t> BinMap_t;
// typedef std::unordered_map<std::size_t, BinPos> keyToBin_t;
// typedef std::unordered_map<BinPos, std::size_t> binToKey_t;


class Box
{
private:
	AtomList_t Atoms;
	BinMap_t BinAtomList;
	// keyToBin_t keyToBin;
	// binToKey_t binToKey;

	ConfigData* data;

	double xBoxSize, yBoxSize, zBoxSize;	//box size in x,y,z
	double xBinSize, yBinSize, zBinSize;	//bin size in x,y,z


public:
	Box(ConfigData* CD) : data(CD)
	{
		BoxSize_t boxSize = data->getBoxSize();
		xBoxSize = boxSize[X];
		yBoxSize = boxSize[Y];
		zBoxSize = boxSize[Z];

		initPos();
		initBins();
		updateMapping();
	}

	// access functions
	// AtomList_t getAtoms() { return Atoms; }

	// other functions
	void updateMapping()
	{
		BinAtomList.clear();

		for (auto atom : Atoms)
		{
			auto binKey = getKey( atom->x(), atom->y(), atom->z() );
			BinAtomList[binKey].push_back(atom);
		}
	}

	void addAtom(double x, double y, double z, double mass, double charge)
	{
		Position* pos = new Position(x, y, z);
		Velocity* v = new Velocity(0.0, 0.0, 0.0);
		Atom* a = new Atom(pos, v, mass, charge);

		// BinAtomList[getKey(x,y,z)].push_back(a); 	//update the atom list for each bin
		Atoms.push_back(a);							//update the entire atom list
	}

	//print atom positions
	void printPositions()
	{
		int atomCount = 0;
		for (auto atom : Atoms)
		{
			atomCount++;
			std::cout << "Atom " << atomCount << ": {";
			atom->print();
			std::cout << "}" << std::endl;
		}
	}

	//print atoms in each bin
	void printBins()
	{
		int binCount = 0;
		for (auto bin_pair : BinAtomList)
		{
			binCount++;
			std::cout << "\nBin " << binCount << ":" << std::endl;

			int atomCount = 0;
			auto bin = bin_pair.second;
			for (auto atom : bin)
			{
				atomCount++;
				std::cout << "\tAtom " << atomCount << ": {";
				atom->print();
				std::cout << "}" << std::endl;
			}
		}
	}

private: 
	// set the initial atom positions based on the the InitPos parameter
	//  in the config file
	bool initPos()
	{
		double x,y,z;
		for (int ion=0; ion < data->getNIons(); ++ion)
		{
			//check data->initPos, as it won't always be random
			if ( data->getInitPos() == "random" )
			{
				x = getRandom(0.0, xBoxSize);
				y = getRandom(0.0, yBoxSize);
				z = getRandom(0.0, zBoxSize);
			}

			//add other InitPos modes

			addAtom(x, y, z, data->getMass(), data->getCharge());
		}
		return true;
	}

	bool initBins()
	{
		xBinSize = xBoxSize / data->getNBinsX();
		yBinSize = yBoxSize / data->getNBinsY();
		zBinSize = zBoxSize / data->getNBinsZ();

		return true;
	}

	double getRandom(double min, double max)
	{
	  timeval t;
	  gettimeofday(&t, nullptr);
	  boost::mt19937 seed( (int)t.tv_usec );
	  boost::uniform_real<> dist(min,max);
	  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > random(seed,dist);
	  return random(); 
	}

	// generate a hash value for the bin key;
	std::size_t getKey(double x, double y, double z)
	{
		int xBin, yBin, zBin;
		xBin = floor(x/xBinSize);
		yBin = floor(y/yBinSize);
		zBin = floor(z/zBinSize);
		// BinPos binXYZ;
		// binXYZ.xBin = floor(x/xBinSize);
		// binXYZ.yBin = floor(y/yBinSize);
		// binXYZ.zBin = floor(z/zBinSize);

		std::size_t key = 0;
		boost::hash_combine(key, xBin);
		boost::hash_combine(key, yBin);
		boost::hash_combine(key, zBin);

		// keyToBin[key] = binXYZ;
		// binToKey[binXYZ] = key;

		// std::cout << "getting key = " << key << "\t";
		// std::cout << "{" << xbin << ", ";
		// std::cout << ybin << ", ";
		// std::cout << zbin << "}" << std::endl;

		return key;
	}

};


#endif /* ATOMS_H */