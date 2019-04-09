#ifndef ATOM_H
#define ATOM_H

// #include <mpi.h>
#include <omp.h>
#include <algorithm>

#include "Types.h"
#include "ConfigData.h"
#include "SourceIncludes.h"
#include "ParticleInfo.h"

#include <boost/serialization/vector.hpp>

class Atom;
typedef std::vector<Atom*> AtomList_t;
// typedef std::set<Atom*> AtomSet_t;

class Atom
{
private:
	friend class boost::serialization::access;
	
	Position* P;			//position
	Velocity* V;			//velocity
	double M; 				//mass
	double Z; 				//charge
	AtomList_t Neighbors; 	//list of neighboring atoms

	template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & P;
        ar & V;
        ar & M;
        ar & Z;
        ar & Neighbors;
    }


public:
	Atom(Position* _P, Velocity* _V, double _M, double _Z) : P(_P), V(_V), M(_M), Z(_Z) {}
	~Atom() {}

	// Access functions
	Position* pos() const { return P; }
	Velocity* vel() const { return V; }
	double x() const { return P->x(); }
	double y() const { return P->y(); }
	double z() const { return P->z(); }
	double mass() const { return M; }
	double charge() const { return Z; }
	double rSquared() const
	{
		return pow(P->x(), 2.0) + pow(P->y(), 2.0) + pow(P->z(), 2.0);
	}

	void printNeighbors()
	{
		if ( !Neighbors.empty() )
		{
			std::cout << "\tNeighbors:" << std::endl;
			for (auto neighbor : Neighbors)
			{
				std::cout << "\t\t";
				neighbor->print();
				std::cout << std::endl;
			}
		}
	}

	void print()
	{
		P->print();
	}


	inline void updateNeighborList(const double &cutoff, const double &Lx, const double &Ly, const double &Lz, const AtomList_t &closeAtoms)
		//std::chrono::duration<double> &t1, std::chrono::duration<double> &t2)
	{

		

		// std::cout << "in updateNeighborList()..." << std::endl;
		Neighbors.clear();
		Atom* atom;
		double dX, dY, dZ;
		double cutoffSquared = pow(cutoff, 2.0);
		// std::cout << "iterating over each of the surrounding atoms..." << std::endl;

		#pragma omp parallel for num_threads(4) shared(Neighbors) private(dX,dY,dZ)
		{
			for (int atomIndex=0; atomIndex < closeAtoms.size(); ++atomIndex)
			{
				atom = closeAtoms[atomIndex];
				//skip the calculations if the atom in question is itself
				if ( this == atom )
					continue;

				// std::cout << "computing dX, dY, dZ, and dR..." << std::endl;
			//	auto start = std::chrono::high_resolution_clock::now();
				dX = dAlpha(x(), atom->x(), Lx);
				dY = dAlpha(y(), atom->y(), Ly);
				dZ = dAlpha(z(), atom->z(), Lz);
			//	auto finish = std::chrono::high_resolution_clock::now();
			//	t1 += finish - start;
		

				// std::cout << "after compute..." << std::endl;
				// std::cout << "dR = " << dR << std::endl;
				// std::cout << "cutoffSquared = " << cutoffSquared << std::endl;
				
			//	start = std::chrono::high_resolution_clock::now();
				if ( dR(dX,dY,dZ) < cutoffSquared )
				{
					#pragma omp critical
					{
						Neighbors.push_back(atom);
					}
				}
			//	finish = std::chrono::high_resolution_clock::now();
			//	t2 += finish - start;
				
			}
		}
		// std::cout << atomIndex << std::endl;
		
		// std::cout << "after iteration..." << std::endl;

		// std::cout << "dAlpha's runtime = " << dAlphaTotalTime.count() << " seconds" << std::endl;
		// std::cout << "dR runtime = " << dRTotalTime.count() << " seconds" << std::endl;
		// int stop;
		// std::cin >> stop;
	}

	inline double dAlpha(const double &x1, const double &x2, const double &L)
	{
		if (x1 < x2)
			return std::min(x2-x1, L-x2-x1);
		else
			return std::min(x1-x2, L-x1-x2);
	}

	inline double dR(const double &dX, const double &dY, const double &dZ)
	{
		return pow(dX, 2) + pow(dY, 2) + pow(dZ, 2);
	}

};


#endif /* ATOM_H */