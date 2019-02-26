#ifndef ATOM_H
#define ATOM_H

#include "Types.h"

// modify later for variable number of dims with this link:
// https://stackoverflow.com/questions/3836648/structure-or-class-with-variable-number-of-members
class Atom
{
private:
	double x;
	double y;
	double z;


public:
	Atom() : x(NOT_SET_DOUBLE), y(NOT_SET_DOUBLE), z(NOT_SET_DOUBLE){}
	Atom(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
	~Atom() {}

	// Access functions
	//MODIFY!! to return bool
	double X() { return x; } 
	double Y() { return y; }
	double Z() { return z; }


};


#endif /* ATOM_H */