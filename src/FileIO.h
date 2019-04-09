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
			BoxSize_t box;
			std::vector<int> BinNums;
			double minBoxDim;
			int minNumBins;
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
					BinNums.push_back( std::atoi(value.c_str()) );
					CD.setNBinsX( std::atoi(value.c_str()) );
				}
				if ( descriptor == "NBinsY" )
				{
					BinNums.push_back( std::atoi(value.c_str()) );
					CD.setNBinsY( std::atoi(value.c_str()) );
				}
				if ( descriptor == "NBinsZ" )
				{
					BinNums.push_back( std::atoi(value.c_str()) );
					CD.setNBinsZ( std::atoi(value.c_str()) );
				}
				if ( descriptor == "BoxSize" )
				{
					minBoxDim = std::atof(value.c_str());
					box.push_back( std::atof(value.c_str()) );
					for (int i = 0; i < CD.getNDim()-1; ++i)
					{
						ss >> std::skipws >> value;
						checkMinValue<double>(value, minBoxDim);
						box.push_back( std::atof(value.c_str()) );
					}
					CD.setBoxSize(box);

					// std::cout << "minBoxSize = " << minBoxSize << std::endl;
				}
				if ( descriptor == "CutoffRadius" )
				{
					double cutoff = std::atof(value.c_str());
					// double maxbinSize = 
					if ( !checkCutoffRadius(cutoff, minBoxDim, getMaxBinSize(box, BinNums)) )
					{
						std::cout << "ERROR: The provided cutoff radius is to large for the box!" << std::endl;
						return false;
					}
					CD.setCutoffRadius(cutoff);
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

	template <class T>
	void checkMinValue(std::string valueStr, T &minValue)
	{
		T value = std::atof(valueStr.c_str());
		if (value < minValue)
		{
			minValue = value;
		}
	}

	bool checkCutoffRadius(double cutoff, double minBoxDim, double maxBinSize)
	{
		// std::cout << "max bin size = " << maxBinSize << std::endl;
		// std::cout << "compare value = " << (minBoxDim-maxBinSize)/2.0 << std::endl;
		if (cutoff > (minBoxDim-maxBinSize)/2.0)
			return false;
		return true;
	}

	double getMaxBinSize(BoxSize_t box, std::vector<int> binNums)
	{
		// auto BoxSize = CD.getBoxSize();
		double binSize;
		double maxBinSize = 0;

		//find maxBinSize
		for (int dim=0; dim < box.size(); ++dim)
		{
			// std::cout << dim << std::endl;
			binSize = box[dim] / (double)binNums[dim];
			// std::cout << "binSize = " << binSize << std::endl;
			if (binSize > maxBinSize)
			{
				maxBinSize = binSize;
				// std::cout << "maxbinsize = " << maxBinSize << std::endl;
			}
		}
		return maxBinSize;
	}


};




#endif /* FILE_IO_H */