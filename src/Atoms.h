#ifndef ATOMS_H
#define ATOMS_H

#include "Atom.h"
#include "Types.h"

typedef std::vector<Atom*> AtomList_t;

class Atoms
{
private:
	AtomList_t Atoms;


public:
	AtomList_t getAtoms() { return Atoms; }

	void addAtom(double x, double y, double z)
	{
		Atom* a = new Atom(x, y, z);
		Atoms.push_back(a);
	}


};


#endif /* ATOMS_H */