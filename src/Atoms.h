#ifndef ATOMS_H
#define ATOMS_H

#include <sys/time.h>
#include <boost/random.hpp>

#include "Atom.h"
#include "Types.h"
#include "ConfigData.h"
#include "ParticleInfo.h"
#include "SourceIncludes.h"

typedef std::vector<Atom*> AtomList_t;

class Atoms
{
private:
	AtomList_t Atoms;


public:
	AtomList_t getAtoms() { return Atoms; }

	void addAtom(double x, double y, double z, double mass, double charge)
	{
		Position* pos = new Position(x, y, z);
		Velocity* v = new Velocity(0.0, 0.0, 0.0);
		Atom* a = new Atom(pos, v, mass, charge);
		Atoms.push_back(a);
	}

	bool init(const ConfigData &CD)
	{
		double x,y,z;
		for (int ion=0; ion < CD.getNIons(); ++ion)
		{
			Box_t boxSize = CD.getBoxSize();

			x = getRandom(0.0, boxSize[0]);
			y = getRandom(0.0, boxSize[1]);
			z = getRandom(0.0, boxSize[2]);

			addAtom(x, y, z, CD.getMass(), CD.getCharge());
		}
		return true;
	}

	//print atom positions
	void printPositions()
	{
		for (auto atom : Atoms)
		{
			atom->pos()->print();
		}
	}

private: 
	double getRandom(double min, double max)
	{
	  timeval t;
	  gettimeofday(&t, nullptr);
	  boost::mt19937 seed( (int)t.tv_usec );
	  boost::uniform_real<> dist(min,max);
	  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > random(seed,dist);
	  return random(); 
	}

};


#endif /* ATOMS_H */