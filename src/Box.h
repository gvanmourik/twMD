#ifndef BOX_H
#define BOX_H

#include <map>
#include <cmath>
#include <algorithm>
#include <sys/time.h>
#include <boost/random.hpp>
#include <boost/functional/hash.hpp>
#include <boost/mpi/collectives.hpp>

#include "Atom.h"
#include "Types.h"
#include "BinPos.h"
#include "ParticleInfo.h"
#include "SourceIncludes.h"

enum dim{X,Y,Z};

// typedef std::vector<Atom*> AtomList_t;
typedef std::map<std::size_t, AtomList_t> BinMap_t;
typedef std::vector<std::size_t> BinIDList_t;
typedef std::vector<BinPos*> BinPosList_t;
// typedef std::unordered_map<std::size_t, BinPos> keyToBin_t;
// typedef std::unordered_map<BinPos, std::size_t> binToKey_t;
// 
void printBinPosList(BinPosList_t &binList)
{
	for (auto binPos : binList)
	{
		binPos->print();
	}
}

void printAtomList(AtomList_t &Atoms)
{
	int atomCount = 0;
	std::cout << "Atoms: " << std::endl;
	for (auto atom : Atoms)
	{
		atomCount++;
		// std::cout << "Bin = " << BinAtomList[getKey(atom->x(), atom->y(), atom->z())]->print() << std::endl;
		std::cout << "Atom " << atomCount << " {" << std::endl; 
		std::cout << "\t(";
		atom->print();
		std::cout << ")" << std::endl;
		std::cout << "\n}" << std::endl;
	}
}


class Box
{
private:
	AtomList_t Atoms;
	BinMap_t BinAtomList;
	BinIDList_t BinIDs;
	BinPosList_t BinPositions;
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

		std::cout << "initializing atom positions..." << std::endl;
		initPos();
		std::cout << "initializing bins..." << std::endl;
		initBins();
		std::cout << "updating the bin mapping..." << std::endl;
		updateMapping();
		std::cout << "building the neighbor list..." << std::endl;
		buildNeighborLists();
	}

	// access functions
	// AtomList_t getAtoms() { return Atoms; }

	// other functions
	void addAtom(double x, double y, double z, double mass, double charge)
	{
		Position* pos = new Position(x, y, z);
		Velocity* v = new Velocity(0.0, 0.0, 0.0);
		Atom* a = new Atom(pos, v, mass, charge);

		// BinAtomList[getKey(x,y,z)].push_back(a); 	//update the atom list for each bin
		Atoms.push_back(a);							//update the entire atom list
	}

	void updateMapping()
	{
		BinAtomList.clear();

		// std::cout << "In updateMapping()..." << std::endl;
		for (auto atom : Atoms)
		{
			auto binKey = getKey( atom->x(), atom->y(), atom->z() );
			BinAtomList[binKey].push_back(atom);
		}
	}

	void buildNeighborLists()
	{
		double cutoffRadius = data->getCutoffRadius();
		
		//debug
		// std::cout << "generating the stencil..." << std::endl;s

		//generate i j k permutation vector
		auto stencil = getStencil(cutoffRadius);

		//debug
		// std::cout << "stencil:" << std::endl;
		// printBinPosList(stencil);
		// std::cout << std::endl;

		int bincount = 0;
	//	int stop;
	//	std::chrono::duration<double> dAlphaTotalTime;
	//	std::chrono::duration<double> dRTotalTime;

		// std::cout << "iterating over the bins..." << std::endl << std::endl;
		for (auto centerBin : BinPositions)
		{
			//from the binPos create a stencil vector of surrounding bins
			bincount++;
			// std::cout << "\ngetting the surroundingBins for bin[" << bincount << "]..." << std::endl;
			auto surroundingBins = getSurroundingBins(centerBin, stencil);
			// std::cout << "getting the atoms within the center bin..." << std::endl;
			auto centerBinAtoms = getBinAtoms(centerBin);

			// std::cout << "\ncenterBin bin = ";
			// centerBin->print();
			// std::cout << "surroundingBins:" << std::endl;
			// printBinPosList(surroundingBins);

			// std::cout << "building the list of atoms in the surrounding bins..." << std::endl;
			//build the surroundingAtomsList
			AtomList_t surroundingAtoms;
			for (auto surroundingBin : surroundingBins)
			{
				auto surBinAtoms = getBinAtoms(surroundingBin);
				if ( surBinAtoms.empty() )
					continue;
				
				// std::cout << "appending to the surrounding atoms list..." << std::endl;
				//append a copy of surBinAtoms into the surroundingAtoms list
				surroundingAtoms.insert( surroundingAtoms.begin(), surBinAtoms.begin(), surBinAtoms.end() );
			}

			// std::cout << "after continue... " << std::endl;
			// if ( centerBinAtoms.empty() )
			// 	std::cout << "no center bin atoms..." << std::endl << std::endl;
			
			// printAtomList(surroundingAtoms);
			// std::cout << "surroundingAtoms.size() = " << surroundingAtoms.size() << std::endl;
			// std::cout << "centerBinAtoms.size() = " << centerBinAtoms.size() << std::endl;
	//		std::cout << "updating the neighbor list for each atom in the center bin..." << std::endl;
			
	//		dAlphaTotalTime = std::chrono::seconds::zero();
	//		dRTotalTime = std::chrono::seconds::zero();


			//mpi test
			int binAtomsCount = centerBinAtoms.size();
			// MPI_Bcast(&binAtomsCount, 1, MPI_INT, 0, MPI_COMM_WORLD); //one-to-all send/recv

			//for each atom in the bin
			#pragma omp parallel for num_threads(4)
			{
				for (int atom=0; atom<binAtomsCount; ++atom)
				{
					centerBinAtoms[atom]->updateNeighborList(cutoffRadius, xBoxSize, yBoxSize, zBoxSize, surroundingAtoms); //dAlphaTotalTime, dRTotalTime);
				}
			}

			// MPI_Reduce(&mypi, &pi, 1, MPI_Type_vector, MPI_SUM, 0, MPI_COMM_WORLD);

				// for (auto atom : centerBinAtoms)
				// {
				// 	// std::cout << "about to update neighbor list..." << std::endl;
				// 	//update neighbor list per atom
				// 	centerBinAtoms[atom]->updateNeighborList(cutoffRadius, xBoxSize, yBoxSize, zBoxSize, surroundingAtoms); //dAlphaTotalTime, dRTotalTime);

				// 	// std::cout << "dAlpha's runtime = " << dAlphaTotalTime.count() << " seconds" << std::endl;
				// 	// std::cout << "dR runtime = " << dRTotalTime.count() << " seconds" << std::endl;
				// 	// int stop;
				// 	// std::cin >> stop;
				// 	// std::cout << "done updating the neighbor list..." << std::endl;
				// }
	//		std::cout << "dAlpha's runtime = " << dAlphaTotalTime.count() << " seconds" << std::endl;
	//		std::cout << "dR runtime = " << dRTotalTime.count() << " seconds" << std::endl;
	//		std::cout << "done updating the neighbor list..." << std::endl;
			// std::cin >> stop;
		}

	}


	//print atom positions
	void printPositions()
	{
		int atomCount = 0;
		std::cout << "Atoms: " << std::endl;
		for (auto atom : Atoms)
		{
			atomCount++;
			// std::cout << "Bin = " << BinAtomList[getKey(atom->x(), atom->y(), atom->z())]->print() << std::endl;
			std::cout << "Atom " << atomCount << " {" << std::endl; 
			std::cout << "\t(";
			atom->print();
			std::cout << ")" << std::endl;
			atom->printNeighbors();
			std::cout << "\n}" << std::endl;
		}
	}

	//print atoms in each bin
	void printBins()
	{
		int binCount = 0;
		for (auto bin_pair : BinAtomList)
		{
			int atomCount = 0;
			auto bin = bin_pair.second;
			if ( !bin.empty() )
			{
				binCount++;
				std::cout << "\nBin " << binCount << ":" << std::endl;
				for (auto atom : bin)
				{
					atomCount++;
					std::cout << "\tAtom " << atomCount << ": {";
					atom->print();
					std::cout << "}" << std::endl;
				}
			}
		}
	}

private: 
	// set the initial atom positions based on the the InitPos parameter
	//  in the config file
	bool initPos()
	{
		double x,y,z;
		for (int ion = 0; ion < data->getNIons(); ++ion)
		{
			//check data->initPos, as it won't always be random
			if ( data->getInitPos() == "random" )
			{
				x = getRandom(0.0, xBoxSize);
				y = getRandom(0.0, yBoxSize);
				z = getRandom(0.0, zBoxSize);
				
				// // test with set positions
				// if (ion == 0)
				// {
				// 	x = 0.1;
				// 	y = 0.1;
				// 	z = 0.1;
				// } 
				// if (ion == 1)
				// {
				// 	x = 0.1;
				// 	y = 1.9;
				// 	z = 0.1;
				// } 
				// if (ion == 2)
				// {
				// 	x = 0.1;
				// 	y = 0.1;
				// 	z = 1.0;
				// } 
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

		// std::cout << "In initBins()..." << std::endl;
		for (int i=0; i < data->getNBinsX(); ++i)
		{
			for (int j=0; j < data->getNBinsY(); ++j)
			{
				for (int k=0; k < data->getNBinsZ(); ++k)
				{
					// BinPos *bin = new BinPos(i,j,k);
					BinPositions.push_back( new BinPos(i,j,k) );
					BinIDs.push_back( getBinID(i,j,k) );
				}
			}
		}

		return true;
	}

	BinPosList_t getStencil(const double &cutoff)
	{
		BinPosList_t perms;

		int dx = ceil(cutoff/xBinSize);
		int dy = ceil(cutoff/yBinSize);
		int dz = ceil(cutoff/zBinSize);

		for (int i=-dx; i <= dx; ++i)
		{
			for (int j=-dy; j <= dy; ++j)
			{
				for (int k=-dz; k <= dz; ++k)
				{
					perms.push_back( new BinPos(i,j,k) );
				}
			}
		}

		return perms;
	}

	BinPosList_t getSurroundingBins(const BinPos *centerBin, const BinPosList_t &Permutations)
	{
		BinPosList_t surroundingBins;
		int xBin, yBin, zBin;
		
		for (auto perm : Permutations)
		{
			xBin = mod( centerBin->xBin + perm->xBin, data->getNBinsX() );
			yBin = mod( centerBin->yBin + perm->yBin, data->getNBinsY() );
			zBin = mod( centerBin->zBin + perm->zBin, data->getNBinsZ() );

			surroundingBins.push_back( new BinPos(xBin, yBin, zBin) );
		}
		return surroundingBins;
	}

	AtomList_t getBinAtoms(const BinPos* bin)
	{
		auto binKey = getBinID( bin->xBin, bin->yBin, bin->zBin);
		return BinAtomList[binKey];
	}

	double getRandom(const double &min, const double &max)
	{
	  	timeval t;
	  	gettimeofday(&t, nullptr);
	  	boost::mt19937 seed( (int)t.tv_usec );
	  	boost::uniform_real<> dist(min,max);
	  	boost::variate_generator<boost::mt19937&, boost::uniform_real<> > random(seed,dist);
	  	return random(); 
	}

	// generate a hash value for the bin key;
	std::size_t getKey(const double &x, const double &y, const double &z)
	{
		int xBin, yBin, zBin;
		xBin = floor(x/xBinSize);
		yBin = floor(y/yBinSize);
		zBin = floor(z/zBinSize);

		// keyToBin[key] = binXYZ;
		// binToKey[binXYZ] = key;

		
		return getBinID(xBin, yBin, zBin);
	}

	std::size_t getBinID(const int &xBin, const int &yBin, const int &zBin)
	{
		std::size_t key = 0;
		boost::hash_combine(key, xBin);
		boost::hash_combine(key, yBin);
		boost::hash_combine(key, zBin);

		// std::cout << "getting key = " << key << "\t";
		// std::cout << "{" << xBin << ", ";
		// std::cout << yBin << ", ";
		// std::cout << zBin << "}" << std::endl;

		return key;
	}

	int mod(const int &a, const int &b)
	{ 
		return (a%b+b)%b; 
	}

};


#endif /* ATOMS_H */