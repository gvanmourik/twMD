#include <string>
#include <iostream>

#include "FileIO.h"

int main(int argc, char** argv)
{
	// std::string configFilePath = argv[1];

	if (argc > 2)
	{
		std::cout << " Too many arguments provided..." << std::endl;
	}
	// if (argc != 2)
	// {
	// 	std::cout << " Please input a file name..." << std::endl;
	// }
	// else if (configFilePath.find(".conf") == std::string::npos)
	// {
	// 	std::cout << " ERROR: Config file is the wrong file type!" << std::endl;
	// 	exit(EXIT_FAILURE);
	// }

	std::cout << "Hello world!" << std::endl;

	return 0;
}
