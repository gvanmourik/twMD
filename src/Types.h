#ifndef TYPES_H
#define TYPES_H

#include "SourceIncludes.h"


#define NOT_SET_DOUBLE -1.0
#define NOT_SET_POS nullptr

// Typedefs
typedef std::vector<double> BoxSize_t;
typedef std::unordered_map<std::string, bool> CheckSet_t;


// X-macros
#define INIT_POS \
    X(random) //\
    // X(<option2>) \
    // X(<option3>)
 
#define X(a) C##a,
	enum Init_Pos { INIT_POS };
#undef X

// to determine valid init positions
#define X(a) #a,
	static std::vector<std::string> InitPosTypes = { INIT_POS };
#undef X

bool isInitPosType(std::string type)
{
	if( !(std::find(InitPosTypes.begin(), InitPosTypes.end(), type) != InitPosTypes.end()) ) 
	{
		std::cout << " ERROR: InitPos in config file is not one of the accepted types!" << std::endl;
	    return false;
	}
	return true;
}

// #define CONFIG_FILE \
//   CONFIG_FILE_MEMBER(NDim, int) \
//   CONFIG_FILE_MEMBER(NIons, int) \
//   CONFIG_FILE_MEMBER(NElec, int) \
//   CONFIG_FILE_MEMBER(BoxSize, int[]) \
//   CONFIG_FILE_MEMBER(radius, double)







#endif /* TYPES_H */