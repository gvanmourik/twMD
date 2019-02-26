#ifndef FILE_IO_H
#define FILE_IO_H

#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <algorithm>
#include <boost/any.hpp>

#include "Types.h"
#include "ConfigData.h"

// to determine valid init positions
#define X(a) #a,
	static std::vector<std::string> InitPosStrings = { INIT_POS };
#undef X


class FileIO
{
private:
	std::string FilePath;


public:
	FileIO() {}
	// FileIO(std::string _FilePath) : FilePath(_FilePath) {}
	~FileIO(){}

	std::string getFilePath() { return FilePath; }
	void setFilePath(std::string _FilePath) { FilePath = _FilePath; }


	bool readConfigFile(std::string _FilePath, ConfigData &CD) 
	{
		if (_FilePath.find(".conf") == std::string::npos)
		{
			std::cout << " ERROR: Cannot read config file. Wrong file type!" << std::endl;
			return false;
		}

		// open file
		std::string currentLine;
		std::ifstream configFile(_FilePath);
		if ( configFile.is_open() )
		{
			std::getline(configFile, currentLine);
			if ( currentLine != "Start MD Configuration File" )
			{
				std::cout << " ERROR: Config file not properly initialized!" << std::endl;
				return false;
			}

			// boost::any value;
			// std::map<std::string, boost::any> Values;
			std::string descriptor, value;
			std::getline(configFile, currentLine);
			while ( currentLine != "End MD Configuration File" )
			{
				std::istringstream ss(currentLine);
				ss >> descriptor >> std::skipws >> value;

				if ( descriptor == "NDim" )
				{
					CD.setNDim( std::atoi(value.c_str()) );
				}
				if ( descriptor == "NIons" )
				{
					CD.setNIons( std::atoi(value.c_str()) );
				}
				if ( descriptor == "NElec" )
				{
					CD.setNElec( std::atoi(value.c_str()) );
				}
				if ( descriptor == "BoxSize" )
				{
					Box_t box;
					box.push_back( std::atoi(value.c_str()) );
					for (int i = 0; i < CD.getNDim()-1; ++i)
					{
						ss >> std::skipws >> value;
						box.push_back( std::atoi(value.c_str()) );
					}
					CD.setBoxSize(box);
				}
				if ( descriptor == "InitPos" )
				{
					if( !(std::find(InitPosStrings.begin(), InitPosStrings.end(), value) != InitPosStrings.end()) ) 
					{
						std::cout << " ERROR: InitPos in config file is not one of the accepted types!" << std::endl;
					    return false;
					}
					CD.setInitPos(value);
				}
				if ( descriptor == "TimeStep" )
				{
					CD.setTimeStep( std::atof(value.c_str()) );
				}

				std::getline(configFile, currentLine);
			}
			configFile.close();
		}
		else
		{
			std::cout << " ERROR: Config file did not open properly!" << std::endl;
			return false;
		}
		
		return true;
	}


};




#endif /* FILE_IO_H */