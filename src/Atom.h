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
	Position* pos() const { return P; }
	Velocity* vel() const { return V; }
	double x() const { return P->x(); }
	double y() const { return P->y(); }
	double z() const { return P->z(); }
	double mass() const { return M; }
	double charge() const { return Z; }
	double r() const
	{
		return sqrt( pow(P->x(), 2) + pow(P->y(), 2) + pow(P->z(), 2) );
	} 


};


#endif /* ATOM_H */