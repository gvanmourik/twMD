#ifndef PARTICLE_INFO_H
#define PARTICLE_INFO_H

#include "Types.h"

// modify later for variable number of dims with this link:
// https://stackoverflow.com/questions/3836648/structure-or-class-with-variable-number-of-members
class ParticleInfo
{
private:
	double X;
	double Y;
	double Z;


public:
	ParticleInfo(double _X, double _Y, double _Z) : X(_X), Y(_Y), Z(_Z) {}
	~ParticleInfo() {}

	// Access functions
	double x() { return X; } 
	double y() { return Y; }
	double z() { return Z; }
	void set_x(double _X) { X = _X; } 
	void set_y(double _Y) { Y = _Y; }
	void set_z(double _Z) { Z = _Z; }

	friend class Position;
	friend class Velocity;
};


class Position : public ParticleInfo 
{
public:
	Position(double _X, double _Y, double _Z) : ParticleInfo(_X, _Y, _Z) {}
	~Position() {}

	virtual void print()
	{
		std::cout << "X = " << X << std::endl;
		std::cout << "Y = " << Y << std::endl;
		std::cout << "Z = " << Z << std::endl;
	}
};

class Velocity : public ParticleInfo 
{
public:
	Velocity(double _X, double _Y, double _Z) : ParticleInfo(_X, _Y, _Z) {}
	~Velocity() {}

	void print()
	{
		std::cout << "Vx = " << X << std::endl;
		std::cout << "Vy = " << Y << std::endl;
		std::cout << "Vz = " << Z << std::endl;
	}
};


#endif /* PARTICLE_INFO_H */