#ifndef TYPES
#define TYPES

#include <vector>
#include <string>
#include <unordered_map>

// Typedefs
typedef std::vector<int> Box_t;


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







#endif /* TYPES */