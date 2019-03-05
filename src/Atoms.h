#ifndef ATOMS_H
#define ATOMS_H

#include "Atom.h"
#include "Types.h"
#include "ParticleInfo.h"
#include "SourceIncludes.h"

typedef std::vector<Atom*> AtomList_t;

class Atoms
{
private:
	AtomList_t Atoms;


public:
	AtomList_t getAtoms() { return Atoms; }

	void addAtom(double x, double y, double z, double charge)
	{
		Position* pos = new Position(x, y, z);
		Velocity* v = new Velocity(0.0, 0.0, 0.0);
		Atom* a = new Atom(pos, v, charge);
		Atoms.push_back(a);
	}


};


#endif /* ATOMS_H */