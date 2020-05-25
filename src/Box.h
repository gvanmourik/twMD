#ifndef BOX_H
#define BOX_H

#include <map>
#include <cmath>
#include <algorithm>
#include <sys/time.h>
#include <boost/random.hpp>
#include <boost/functional/hash.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/serialization/map.hpp>

#include "Atom.h"
#include "Types.h"
#include "BinPos.h"
#include "Electron.h"
#include "ParticleInfo.h"
#include "SourceIncludes.h"

enum dim{X,Y,Z};

typedef std::map<std::size_t, AtomList_t> BinMap_t;
typedef std::vector<std::size_t> BinIDList_t;
typedef std::vector<BinPos*> BinPosList_t;
typedef std::map<std::size_t, CD_t> OverlapMap_t;

CD_t j(0.0, 1.0);
double SQRT_SIX = sqrt(6.0);


class Box
{
private:
	AtomList_t Atoms;
	ElectronList_t Electrons;
	BinMap_t BinAtomList;
	BinIDList_t BinIDs;
	BinPosList_t BinPositions;
	OverlapMap_t GaussianOverlapMap;

	//mpi
	boost::mpi::communicator World;

	//config file data
	ConfigData* data;

	//box and bin sizes
	double xBoxSize, yBoxSize, zBoxSize;	//box size in x,y,z
	double xBinSize, yBinSize, zBinSize;	//bin size in x,y,z


public:
	Box(ConfigData* CD, boost::mpi::communicator world) : data(CD), World(world)
	{
		BoxSize_t boxSize = data->getBoxSize();
		xBoxSize = boxSize[X];
		yBoxSize = boxSize[Y];
		zBoxSize = boxSize[Z];

		if (World.rank() == 0)
		{
			std::cout << "initializing atoms..." << std::endl;
			initAtoms();
			std::cout << "initializing electrons..." << std::endl;
			initElectrons();
			std::cout << "initializing bins..." << std::endl;
			initBins();
			std::cout << "updating the bin mapping..." << std::endl;
			updateMapping();
			
			// std::cout << "computing gaussian overlap..." << std::endl;
			// computeGaussianOverlap();

			std::cout << "building the neighbor list..." << std::endl;
		}
		buildNeighborLists();


		if (World.rank() == 0)
		{
			std::cout << "computing gaussian overlap..." << std::endl;
		}
		computeGaussianOverlap();

		V_ii();
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


	//**********************************************************************	
	void buildNeighborLists()
	{
		const int parentProcess = 0;
		const int scatterProcess = parentProcess; //scatter from parent
		double cutoffRadius;
		BinPosList_t stencil;
		BinPosList_t recvBins;
		std::vector<BinPosList_t> binSections( World.size() );

		
		if (World.rank() == 0) 
		{
			cutoffRadius = data->getCutoffRadius();
			stencil = getStencil(cutoffRadius);

			int start, end;
			int numBins = data->getNBins();
			int numProcesses = World.size();
			int myFloor = floor(numBins/numProcesses);
			int stride = ceil(numBins/numProcesses);

			// split bin positions based on number of processors
			for (int process=0; process<numProcesses; ++process)
			{
				//determine end and start indicies of bin sections
				start = process*stride;
				if (myFloor != stride && process == numProcesses-1)
					end = numBins % start;
				else
					end = start + stride;

				//place bins in sections
				BinPosList_t binSection;
				for (int i=start; i<end; ++i)
				{
					binSection.push_back( BinPositions[i] );
				}
				binSections[process] = binSection;
			}
  		}

		//bcast stencil and bin positions
		broadcast(World, stencil, parentProcess);
		broadcast(World, cutoffRadius, parentProcess);
		// broadcast(World, data, parentProcess);
		// broadcast(World, xBoxSize, parentProcess);
		// broadcast(World, yBoxSize, parentProcess);
		// broadcast(World, zBoxSize, parentProcess);
		// broadcast(World, BinAtomList, parentProcess);
		std::cout << "in Process[" << World.rank() << "]..." << std::endl;


		//scatter bin sections to each process
		boost::mpi::scatter(World, binSections, recvBins, scatterProcess);

		std::cout << "recvBins.size() for process[" << World.rank() << "] = " << recvBins.size() << std::endl;

		//********************************************************************************
		//do stuff with recvBins
		for (auto centerBin : recvBins)
		{
			//debug
			// std::cout << "getting surroundingBins process[" << World.rank() << "]..." << std::endl;

			//from the binPos create a stencil vector of surrounding bins
			auto surroundingBins = getSurroundingBins(centerBin, stencil);
			auto centerBinAtoms = getBinAtoms(centerBin);

			//debug
			// std::cout << "building the surroundingAtomsList process[" << World.rank() << "]..." << std::endl;

			//build the surroundingAtomsList
			AtomList_t surroundingAtoms;
			for (auto surroundingBin : surroundingBins)
			{
				auto surBinAtoms = getBinAtoms(surroundingBin);
				if ( surBinAtoms.empty() )
					continue;
				
				surroundingAtoms.insert( surroundingAtoms.begin(), surBinAtoms.begin(), surBinAtoms.end() );
			}


			//debug
			// std::cout << "updating the neighbor list process[" << World.rank() << "]..." << std::endl;

			int binAtomsCount = centerBinAtoms.size();
			int thread_count = data->getThreadCount();
			//for each atom in the bin
			#pragma omp parallel for num_threads(thread_count)
			{
				for (int atom=0; atom<binAtomsCount; ++atom)
				{
					centerBinAtoms[atom]->updateNeighborList(cutoffRadius, xBoxSize, yBoxSize, zBoxSize, 
						surroundingAtoms, thread_count);
				}
			}
		}
		//********************************************************************************

		//-->appears to be no need to gather

	}
	//**********************************************************************


	//**********************************************************************
	//assumes that buildNeighborLists() has been called
	void computeGaussianOverlap()
	{
		const int parentProcess = 0;
		const int scatterProcess = parentProcess; //scatter from parent
		CD_t o;
		std::size_t key, conjKey;
		AtomList_t neighbors;
		GaussianList_t gaussians1, gaussians2;
		AtomList_t recvAtomList;
		std::vector<AtomList_t> AtomSections( World.size() );


		// std::cout << "in computeGaussianOverlap()..." << std::endl;

		//split up Atoms into sections
		if (World.rank() == 0) 
		{
			int start, end;
			int numAtoms = Atoms.size();
			int numProcesses = World.size();
			int myFloor = floor(numAtoms/numProcesses);
			int stride = ceil(numAtoms/numProcesses);

			// split bin positions based on number of processors
			for (int process=0; process<numProcesses; ++process)
			{
				//determine end and start indicies of bin sections
				start = process*stride;
				if (myFloor != stride && process == numProcesses-1)
					end = numAtoms % start;
				else
					end = start + stride;

				//place bins in sections
				AtomList_t atomSection;
				for (int i=start; i<end; ++i)
				{
					atomSection.push_back( Atoms[i] );
				}
				AtomSections[process] = atomSection;
			}
  		}
  		std::cout << "in Process[" << World.rank() << "]..." << std::endl;

  		//scatter atom sections to processes
  		boost::mpi::scatter(World, AtomSections, recvAtomList, scatterProcess);
  		std::cout << "recvAtomList.size() for process[" << World.rank() << "] = " << recvAtomList.size() << std::endl;

  		//perform the computation
		for (auto atom1 : recvAtomList)
		{
			neighbors = atom1->getNeighbors();
			gaussians1 = atom1->getGaussians();
			for (auto atom2 : neighbors)
			{
				gaussians2 = atom2->getGaussians();
				for (auto g1 : gaussians1)
				{
					for (auto g2 : gaussians2)
					{
						//get keys for overlap map
						key = getOverlapKey(atom1, g1, atom2, g2);
						conjKey = getOverlapKey(atom2, g2, atom1, g1);

						// std::cout << "key = " << key << std::endl;

						//check if key exists in map
						if ( !inOverlapMap(key) )
						{
							//if not, add to both
							o = overlap(g1, g2);
							GaussianOverlapMap[key] = o;
							GaussianOverlapMap[conjKey] = std::conj(o);
						}
					}
				}
			}
		}
	}

	//compute the gaussian overlap between two gau
	CD_t overlap(Gaussian* g1, Gaussian* g2)
	{
		// std::cout << "in overlap()..." << std::endl;

		// r12 = 1.0
		// g1 = 1.0
		// g2 = 2.0
		// 1.8 + j0.55
		// 1.89 + j0.52


		// double r12 = 1.0;
		
		// double CR1 = 1.0;
		// double CI1 = 1.0;
		// double S1 = 1.0;
		// double Rho1 = 1.0;

		// double CR2 = 2.0;
		// double CI2 = 2.0;
		// double S2 = 2.0;
		// double Rho2 = 2.0;
		
		// g1->setCr(1.0);
		// g1->setCi(1.0);
		// g1->setS(1.0);
		// g1->setRho(1.0);

		// g2->setCr(1.0);
		// g2->setCi(1.0);
		// g2->setS(1.0);
		// g2->setRho(1.0);

		// std::cout << "radius^2 = " << r12(g1,g2) << std::endl;


		CD_t result = ( 6.0 * SQRT_SIX * (g1->Ci()+j*g1->Cr()) * (j*g2->Ci()+g2->Cr()) * g1->S() * g2->S() ) /
			(
				(exp( (r12(g1->pos(),g2->pos()) * (3.0 + j*2.0*g1->Rho()*g1->S()) * (j*3.0 + 2.0*g2->Rho()*g2->S())) /
				  (4.0 * (j*3.0 * pow(g2->S(),2.0) - 2.0*g1->Rho()*g1->S()*pow(g2->S(),2.0) + 
				   	pow(g1->S(),2.0) * (j*3.0 + 2.0*g2->Rho()*g2->S()))) )
				) *
				(sqrt( j*(-2.0)*g2->Rho()*g1->S() + (3.0*g1->S())/g2->S() + j*2.0*g1->Rho()*g2->S() + 
					(3.0*g2->S())/g1->S() )
				) *
				(j*3.0*pow(g2->S(),2.0) - 2.0*g1->Rho()*g1->S()*pow(g2->S(),2.0) + pow(g1->S(),2.0) * 
					(j*3.0 + 2.0*g2->Rho()*g2->S()) 
				)
			);
		
		// std::cout << "result = " << result << std::endl;
		return result;
	}

	
	//**********************************************************************


	//**********************************************************************
	//returns the total interactive potential across all ions
	double V_ii()
	{
		double r;
		double sum = 0.0;
		double Z = data->getCharge();
		AtomList_t Neighbors;

		for (auto atom : Atoms)
		{
			Neighbors = atom->getNeighbors();
			for (auto neighbor : Neighbors)
			{
				r = r12(atom->pos(), neighbor->pos());
				sum += 1.0/sqrt(r);
			}
		}

		std::cout << "V_ii = " << pow(Z,2.0) / 2.0 * sum << std::endl;
		return pow(Z,2.0) / 2.0 * sum;
	}

	//TODO!!!!!!!!!!
	void V_ie() {}
	void V_ee() {}

	//kinetic energy between two gaussians...currently does not find the total kinetic energy across the electrons
	CD_t T_e(Gaussian* g1, Gaussian* g2) 
	{
		double h = 1.0;
		double m = 1.0;
		CD_t gamma1 = gamma(g1);
		CD_t gamma2 = gamma(g2);
		
		return (-1.0) * gamma2 * pow(h,2.0) * pow(m,-1) * std::conj(gamma1) *
			pow((gamma2+std::conj(gamma1)),-2.0) * 
			((-3.0)*gamma2+((-3.0)+2.0*gamma2*r12(g1->pos(),g2->pos()))*std::conj(gamma1));
	}

	//helper function for T_e()
	CD_t gamma(Gaussian* g)
	{
		return (3.0/(4.0*pow(g->S(),2.0))) - ((j*g->Rho())/(2.0*g->S()));
	}

	//kinetic energy over all ions
	void T_i() 
	{
		double sum = 0.0;
		double m = data->getMass();
		for (auto atom : Atoms)
		{
			// v() returns the velocity squared
			sum += atom->v();
		}
		return m / 2.0 * sum;
	}

	//**********************************************************************


	//returns the distance squared between two positions
	double r12(Position* r1, Position* r2)
	{
		return pow(r1->x() - r2->x(), 2.0) + 
			pow(r1->y() - r2->y(), 2.0) + 
			pow(r1->z() - r2->z(), 2.0);
	}

	//print atom positions
	void printAtoms()
	{
		int atomCount = 0;
		std::cout << "Atoms: " << std::endl;
		for (auto atom : Atoms)
		{
			atomCount++;
			// std::cout << "Bin = " << BinAtomList[getKey(atom->x(), atom->y(), atom->z())]->print() << std::endl;
			std::cout << "Atom " << atomCount << "{" << std::endl; 
			std::cout << "\t(";
			atom->print();
			std::cout << ")" << std::endl;
			atom->printNeighbors();
			std::cout << "\n}" << std::endl;
		}
	}

	//print electron positions
	void printElectrons()
	{
		int electronCount = 0;
		std::cout << "Electrons: " << std::endl;
		for (auto e : Electrons)
		{
			electronCount++;
			// std::cout << "Bin = " << BinAtomList[getKey(atom->x(), atom->y(), atom->z())]->print() << std::endl;
			std::cout << "Electron " << electronCount << "{" << std::endl;
			e->print();
			std::cout << "}" << std::endl;
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

	void printOverlapValues()
	{
		std::cout << "Overlap Values:" << std::endl;
		for (auto overlap_pair : GaussianOverlapMap)
		{
			std::cout << overlap_pair.second << std::endl;
		}
		std::cout << std::endl;
	}

private: 
	// set the initial atom positions based on the the InitPos parameter
	//  in the config file
	bool initAtoms()
	{
		double x,y,z;
		for (int ion = 0; ion < 6/*data->getNIons()*/; ++ion)
		{
			//check data->initPos, as it won't always be random
			if ( data->getInitPos() == "random" )
			{
				x = getRandom(0.0, xBoxSize);
				y = getRandom(0.0, yBoxSize);
				z = getRandom(0.0, zBoxSize);
				
				// test with set positions
				if (ion == 0)
				{
					x = 0.1;
					y = 0.1;
					z = 0.1;
				} 
				if (ion == 1)
				{
					x = 0.1;
					y = 1.9;
					z = 0.1;
				} 
				if (ion == 2)
				{
					x = 0.1;
					y = 0.1;
					z = 1.0;
				}
				if (ion == 3)
				{
					x = 0.1;
					y = 0.1;
					z = 1.1;
				}
				if (ion == 4)
				{
					x = 0.1;
					y = 0.1;
					z = 1.2;
				} 
				if (ion == 5)
				{
					x = 0.1;
					y = 0.1;
					z = 2.2;
				} 
			}

			//add other InitPos modes

			addAtom(x, y, z, data->getMass(), data->getCharge());
		}
		return true;
	}

	//assumes that initAtoms() had already been called
	bool initElectrons()
	{
		int numElec = data->getNElec();
		int thread_count = data->getThreadCount();

		#pragma omp parallel for num_threads(thread_count)
		{
			for (int i=0; i < numElec; ++i)
			{
				Electron* e = new Electron();
				Gaussian* g;

				//add gaussian to each atom in the box
				for (auto atom : Atoms)
				{	
					g = new Gaussian(e,
						data->getCr(),
						data->getCi(),
						data->getS(), 
						data->getRho(),
						atom->pos()
					);
					e->addGaussian(g);
					atom->addGaussian(g); //add to the gaussian list associated with each atom
				}

				//add to the electron list
				#pragma omp critical
				{
					Electrons.push_back(e);
				}
			}
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

	BinPosList_t getSurroundingBins(BinPos *centerBin, const BinPosList_t &Permutations)
	{
		BinPosList_t surroundingBins;
		int xBin, yBin, zBin;
		
		for (auto perm : Permutations)
		{
			xBin = mod( centerBin->xBin() + perm->xBin(), data->getNBinsX() );
			yBin = mod( centerBin->yBin() + perm->yBin(), data->getNBinsY() );
			zBin = mod( centerBin->zBin() + perm->zBin(), data->getNBinsZ() );

			surroundingBins.push_back( new BinPos(xBin, yBin, zBin) );
		}
		return surroundingBins;
	}

	AtomList_t getBinAtoms(BinPos* bin)
	{
		auto binKey = getBinID( bin->xBin(), bin->yBin(), bin->zBin());
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

	std::size_t getOverlapKey(Atom* atom1, Gaussian* g1, Atom* atom2, Gaussian* g2)
	{
		std::size_t key = 0;
		boost::hash_combine(key, atom1);
		boost::hash_combine(key, g1);
		boost::hash_combine(key, atom2);
		boost::hash_combine(key, g2);

		return key;
	}

	bool inOverlapMap(std::size_t key)
	{
		if ( GaussianOverlapMap.find(key) == GaussianOverlapMap.end() ) {
			return false;
		}
		return true;
	}

	int mod(const int &a, const int &b)
	{ 
		return (a%b+b)%b; 
	}

};


#endif /* ATOMS_H */