#ifndef ATOM_H
#define ATOM_H

// #include <mpi.h>
#include <omp.h>
#include <algorithm>
#include <boost/serialization/vector.hpp>

#include "Types.h"
#include "Gaussian.h"
#include "ConfigData.h"
#include "SourceIncludes.h"
#include "ParticleInfo.h"


class Atom;
typedef std::vector<Atom*> AtomList_t;

class Atom
{
private:
	friend class boost::serialization::access;
	
	Position* P;				//position
	Velocity* V;				//velocity
	double M; 					//mass
	double Z; 					//charge
	AtomList_t Neighbors; 		//list of neighboring atoms
	GaussianList_t Gaussians;	//list of gaussians

	template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & P;
        ar & V;
        ar & M;
        ar & Z;
        ar & Neighbors;
        ar & Gaussians;
    }


public:
	Atom() {}
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
	int getGaussianCount() const { return Gaussians.size(); }
	AtomList_t getNeighbors() const { return Neighbors; }
	GaussianList_t getGaussians() const { return Gaussians; }

	

	//other functions
	void addGaussian(Gaussian* g) { Gaussians.push_back(g); }

	inline void updateNeighborList(const double &cutoff, const double &Lx, const double &Ly, const double &Lz, 
		const AtomList_t &closeAtoms, const int thread_count)
	{
		//debug
		// std::cout << "in updateNeighborList process[" << "]..." << std::endl;
		
		Neighbors.clear();
		Atom* atom;
		double dX, dY, dZ;
		double cutoffSquared = pow(cutoff, 2.0);

		#pragma omp parallel for num_threads(thread_count) shared(Neighbors) private(dX,dY,dZ)
		{
			for (int atomIndex=0; atomIndex < closeAtoms.size(); ++atomIndex)
			{
				atom = closeAtoms[atomIndex];
				
				//skip the calculations if the atom in question is itself
				if ( this == atom )
					continue;

				dX = dAlpha(x(), atom->x(), Lx);
				dY = dAlpha(y(), atom->y(), Ly);
				dZ = dAlpha(z(), atom->z(), Lz);

				if ( dR(dX,dY,dZ) < cutoffSquared )
				{
					#pragma omp critical
					{
						Neighbors.push_back(atom);
					}
				}
			}
		}
	}

	inline double dAlpha(const double &x1, const double &x2, const double &L)
	{
		if (x1 < x2)
			return std::min(x2-x1, L-x2+x1);
		else
			return std::min(x1-x2, L-x1+x2);
	}

	inline double dR(const double &dX, const double &dY, const double &dZ)
	{
		return pow(dX, 2) + pow(dY, 2) + pow(dZ, 2);
	}


	//print functions
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
		if ( !Gaussians.empty() )
		{
			std::cout << "\tGaussians:" << std::endl;
			for (auto g : Gaussians)
			{
				std::cout << "\t\t";
				g->print();
				std::cout << std::endl;
			}
		}
	}

	void print()
	{
		P->print();
	}

};


#endif /* ATOM_H */