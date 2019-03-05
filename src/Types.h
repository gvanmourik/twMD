#ifndef TYPES_H
#define TYPES_H

#include "SourceIncludes.h"


#define NOT_SET_DOUBLE -1.0
#define NOT_SET_POS nullptr

// Typedefs
typedef std::vector<double> Box_t;
typedef std::unordered_map<std::string, bool> CheckSet_t;


// X-macros
#define INIT_POS \
    X(random) //\
    // X(<option2>) \
    // X(<option3>)
 
#define X(a) C##a,
	enum Init_Pos { INIT_POS };
#undef X

// #define CONFIG_FILE \
//   CONFIG_FILE_MEMBER(NDim, int) \
//   CONFIG_FILE_MEMBER(NIons, int) \
//   CONFIG_FILE_MEMBER(NElec, int) \
//   CONFIG_FILE_MEMBER(BoxSize, int[]) \
//   CONFIG_FILE_MEMBER(radius, double)







#endif /* TYPES_H */