# Project 1

### Prerequisites:
* [CMake](https://cmake.org/download/) (Version >= 3.0)
* cmake (configured as an available bash command, to do so follow the instructions below)

To add cmake as a macOS bash command:
1) Open the paths file with:
```
sudo <text editor> /etc/paths
```
2) Add the following line to the file:
```
/Applications/CMake.app/Contents/bin
```
3) Save and close the file.
4) Close and open the terminal window again.
5) Verify that the cmake command is working with:
```
cmake --version
```


### Building:
```
<run instructions>
```

### Running:
Navigate to the build directory.

To run the test config file:
```
./bin/sim01 ../src/testConfig.conf
```

To run other config files:
```
./bin/sim01 <file_path>
```

### Notes:


