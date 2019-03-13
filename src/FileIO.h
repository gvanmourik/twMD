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

// // to determine valid init positions
// #define X(a) #a,
// 	static std::vector<std::string> InitPosTypes = { INIT_POS };
// #undef X


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
		// check if the file has the correct extension
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
				if ( descriptor == "NBinsX" )
				{
					CD.setNBinsX( std::atoi(value.c_str()) );
				}
				if ( descriptor == "NBinsY" )
				{
					CD.setNBinsY( std::atoi(value.c_str()) );
				}
				if ( descriptor == "NBinsZ" )
				{
					CD.setNBinsZ( std::atoi(value.c_str()) );
				}
				if ( descriptor == "BoxSize" )
				{
					BoxSize_t box;
					box.push_back( std::atof(value.c_str()) );
					for (int i = 0; i < CD.getNDim()-1; ++i)
					{
						ss >> std::skipws >> value;
						box.push_back( std::atof(value.c_str()) );
					}
					CD.setBoxSize(box);
				}
				if ( descriptor == "CutoffRadius" )
				{
					CD.setCutoffRadius( std::atof(value.c_str()) );
				}
				if ( descriptor == "InitPos" )
				{
					if ( !isInitPosType(value) )
						return false;
					CD.setInitPos(value);
				}
				if ( descriptor == "TimeStep" )
				{
					CD.setTimeStep( std::atof(value.c_str()) );
				}

				// gaussian descriptors
				if ( descriptor == "GaussianMagReal" )
				{
					CD.setCr( std::atof(value.c_str()) );
				}
				if ( descriptor == "GaussianMagImag" )
				{
					CD.setCi( std::atof(value.c_str()) );
				}
				if ( descriptor == "GuassianWidth" )
				{
					CD.setS( std::atof(value.c_str()) );
				}
				if ( descriptor == "GaussianMomentum" )
				{
					CD.setRho( std::atof(value.c_str()) );
				}

				// atom descriptors
				if ( descriptor == "AtomicCharge" )
				{
					CD.setZ( std::atof(value.c_str()) );
				}
				if ( descriptor == "AtomicMass" )
				{
					CD.setM( std::atof(value.c_str()) );
				}


				std::getline(configFile, currentLine);
			}
			configFile.close();
		}
		else
		{
			std::cout << " ERROR: Config file did not open properly! Check to make sure it exists." << std::endl;
			return false;
		}
		
		return true;
	}


};




#endif /* FILE_IO_H */