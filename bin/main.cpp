#include <string>
#include <iostream>
#include <boost/mpi.hpp>

#include "Box.h"
#include "FileIO.h"


bool checkArgCount(int argc);

int main(int argc, char** argv)
{
	//Setup MPI
	boost::mpi::environment env(argc, argv);
	boost::mpi::communicator world;

	if ( !checkArgCount(argc) )
		return -1;
	std::string configFilePath = argv[1];

	//Read-in file
	FileIO file;
	ConfigData configData;
	if ( !file.readConfigFile(configFilePath, configData) )
		return -1;
	
	//display config data
	if (world.rank() == 0)
		configData.print();

	//init atom positions
	Box* box = new Box(&configData, world);
	
	if (world.rank() == 0)
	{
		box->printAtoms();
		box->printElectrons();
		box->printOverlapValues();
		// box->printBins();
	}
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