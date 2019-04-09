#include <string>
#include <iostream>
#include <boost/mpi.hpp>

#include "FileIO.h"
#include "Box.h"
#include "Electrons.h"


bool checkArgCount(int argc);

int main(int argc, char** argv)
{
	if ( !checkArgCount(argc) )
		return -1;
	std::string configFilePath = argv[1];


	//Read-in file
	FileIO file;
	ConfigData configData;
	if ( !file.readConfigFile(configFilePath, configData) )
		return -1;
	
	//display config data
	configData.print();

	//Setup MPI
	boost::mpi::environment env(argc, argv);
	boost::mpi::communicator world;

	//init atom positions
	Box* box = new Box(&configData);
	// box->initPos();
	//debug with print
	
	// box->printPositions();
	// box->printBins();


	return 0;
}

bool checkArgCount(int argc)
{
	if (argc > 2)
	{
		std::cout << " ERROR: Too many arguments provided!" << std::endl;
		return false;
	}
	if (argc != 2)
	{
		std::cout << " ERROR: No config file path provided!" << std::endl;
		return false;
	}
	return true;
}
