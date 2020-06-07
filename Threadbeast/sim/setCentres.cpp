#include <morph/tools.h>
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
#include <string>
#include <cctype>
// #include <boost/math/special_functions/bessel.hpp>
#define PI 3.1415926535897932
#define NUMPOINTS 41 //just the A-E rows.
//#define NUMPOINTS 79 //just the A-E rows.

using std::vector;
using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::runtime_error;
using morph::HdfData;
using morph::Tools;
using namespace std;




int main (int argc, char **argv)
{
    pair <float, float> centres[NUMPOINTS]; //seed points for regions
    if (argc < 5) {
      std::cout << "not enough arguments " << argc << endl;
      return -1;
    }
    double maxX = stod(argv[1]); 
    double minX = stod(argv[2]); 
    double maxY = stod(argv[3]);
	double minY = stod(argv[4]);

    ofstream afile ( "./centres.h");

    cout << "numpoints" << NUMPOINTS << " maxX " << maxX << " minX " << minX << " maxY " << maxY << " minY " << minY << endl;


   for (int i=0;i<NUMPOINTS;i++) {
       double choice = morph::Tools::randDouble();
       if ((0 < choice) && (choice <= 0.25)) {
           centres[i].first = (morph::Tools::randDouble()) * maxX;
           centres[i].second = (morph::Tools::randDouble()) * maxY  ;
	   }	
       else if ((0.25 < choice) && (choice <= 0.5))
	   {
           centres[i].first = (morph::Tools::randDouble()) * maxX;
           centres[i].second = (morph::Tools::randDouble()) * minY  ;
	   }
	   else if ((0.5 < choice) && (choice  <= 0.75)) {
           centres[i].first = (morph::Tools::randDouble()) * minX;
           centres[i].second = (morph::Tools::randDouble()) * maxY  ;
	   }	
       else 
	   {
           centres[i].first = (morph::Tools::randDouble()) * minX;
           centres[i].second = (morph::Tools::randDouble()) * minY;
	   }
    } //end of setting of random value

    for (int i=0; i < NUMPOINTS; i++) {
	  string sindex = to_string(i);
	  string centres1 = " centres[" + sindex + "].first = ";
	  string  centres2 = " centres[" + sindex + "].second = ";
      afile << centres1 << centres[i].first << " ; " << centres2 << centres[i].second << " ; " << endl;
    }

    return 0;
};
