#include <string>
#include <iostream>

#include "FileIO.h"
#include "Atoms.h"
#include "Electrons.h"


bool checkArgCount(int argc);

int main(int argc, char** argv)
{
	if ( !checkArgCount(argc) )
		return -1;
	std::string configFilePath = argv[1];


	FileIO file;
	ConfigData configData;
	file.readConfigFile(configFilePath, configData);
	configData.print();


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
