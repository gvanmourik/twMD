#ifndef ATOM_H
#define ATOM_H

#include "Types.h"
#include "SourceIncludes.h"
#include "ParticleInfo.h"


class Atom
{
private:
	Position* P;		//position
	Velocity* V;		//velocity
	double M; 			//mass
	double Z; 			//charge


public:
	Atom(Position* _P, Velocity* _V, double _M, double _Z) : P(_P), V(_V), M(_M), Z(_Z) {}
	~Atom() {}

	// Access functions
	Position* pos() { return P; }
	Velocity* vel() { return V; }
	double x() { return P->x(); }
	double y() { return P->y(); }
	double z() { return P->z(); }
	double mass() { return M; }
	double charge() { return Z; }
	double r() 
	{
		return sqrt( pow(P->x(), 2) + pow(P->y(), 2) + pow(P->z(), 2) );
	} 


};


#endif /* ATOM_H */