/*
#include <morph/tools.h>
#include <morph/HexGrid.h>
#include <morph/HdfData.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <math.h>
#include <random>
#include <algorithm>
#include <hdf5.h>
#include <unistd.h>
#include <bits/stdc++.h>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#include "analysis.h"
#include "region.h"
*/
#include "hexTest.h"
using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::runtime_error;
using morph::HexGrid;
using morph::HdfData;
using morph::Tools;
using namespace std;


int main (int argc, char **argv)
{
    if (argc < 9) {
      std::cout << "not enough arguments" << argc << endl;
      return -1;
    }
    //const char* logdir = "cd logs";
   // system (logdir);
    string logpath = argv[1];
    string commandStem;
    bool overwrite_logs = true;
    const char* command;
    double dt = stod(argv[2]); //timesetp passed to M.step
    double Dn = stod(argv[3]); //Dn diffusion passed to M.step
    double Dchi = stod(argv[4]); //Dchi chemotaxis passed to M.step
    double Dc = stod(argv[5]);
    int numsteps = atoi(argv[6]); //length of integration 
    int numprint = atoi(argv[7]); //frequency of printing
    int Lcontinue = atoi(argv[8]); //logical to determine if coldstart
      /*
     * Now create a log directory if necessary, and exit on any
     * failures.
     */
    if (morph::Tools::dirExists (logpath) == false) {
        morph::Tools::createDir (logpath);
        if (morph::Tools::dirExists (logpath) == false) {
            cerr << "Failed to create the logpath directory "
                 << logpath << " which does not exist."<< endl;
            return 1;
        }
    } else {
        // Directory DOES exist. See if it contains a previous run and
        // exit without overwriting to avoid confusion.
        if (overwrite_logs == false
            && (morph::Tools::fileExists (logpath + "/params.json") == true
                || morph::Tools::fileExists (logpath + "/guidance.h5") == true
                || morph::Tools::fileExists (logpath + "/positions.h5") == true)) {
            cerr << "Seems like a previous simulation was logged in " << logpath << ".\n"
                 << "Please clean it out manually, choose another directory or set\n"
                 << "overwrite_logs to true in your parameters config JSON file." << endl;
            return 1;
        }
    }
// initialise Region class setting scale
    Region M(6,logpath);
	M.listHex();

    return 0;
};
