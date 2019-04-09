#ifndef PARTICLE_INFO_H
#define PARTICLE_INFO_H

#include "Types.h"

// modify later for variable number of dims with this link:
// https://stackoverflow.com/questions/3836648/structure-or-class-with-variable-number-of-members
class ParticleInfo
{
private:
	friend class boost::serialization::access;

	double X;
	double Y;
	double Z;

	template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & X;
        ar & Y;
        ar & Z;
    }


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

	// // generate a hash value for the particle;
	// std::size_t getKey()
	// {
	// 	std::size_t seed = 0;
	// 	boost::hash_combine(seed, X);
	// 	boost::hash_combine(seed, Y);
	// 	boost::hash_combine(seed, Z);

	// 	return seed;
	// }

	friend class Position;
	friend class Velocity;
};


class Position : public ParticleInfo 
{
public:
	Position() : ParticleInfo(0.0, 0.0, 0.0) {}
	Position(double _X, double _Y, double _Z) : ParticleInfo(_X, _Y, _Z) {}
	~Position() {}

	virtual void print()
	{
		std::cout << "X = " << X << ", ";
		std::cout << "Y = " << Y << ", ";
		std::cout << "Z = " << Z;
	}
};

class Velocity : public ParticleInfo 
{
public:
	Velocity() : ParticleInfo(0.0, 0.0, 0.0) {}
	Velocity(double _X, double _Y, double _Z) : ParticleInfo(_X, _Y, _Z) {}
	~Velocity() {}

	void print()
	{
		std::cout << "\tVx = " << X << std::endl;
		std::cout << "\tVy = " << Y << std::endl;
		std::cout << "\tVz = " << Z << std::endl;
	}
};


#endif /* PARTICLE_INFO_H */