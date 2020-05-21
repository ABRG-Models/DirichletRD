#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;
#define NUMPOINTS 41 //just the A-E rows.

int  main(void) {
    ofstream afile ("./centresResized.h");
    pair <double, double> centres[NUMPOINTS]; //seed points for regions
    

    #include "centres.h"
	double sum1, sum2 = 0;
    for (int i=0;i<NUMPOINTS;i++) {
        sum1 += centres[i].first;
	    sum2 += centres[i].second;
	}
    sum1 = sum1 / (1.0*NUMPOINTS);
	sum2 = sum2 / (1.0*NUMPOINTS);
    for (int i=0;i<NUMPOINTS;i++) {
		std::pair<double,double> tempPair((centres[i].first - sum1),(centres[i].second - sum2));
		centres[i].first = tempPair.first/4.5 + 0.10;
		centres[i].second = tempPair.second/4.5 - 0.15;
	}
   for (int count = 0; count <NUMPOINTS;count++) {
	    std::string str = to_string(count);
	    afile << "centres[" + str + "].first= "<<    centres[count].first << " ; centres[" + str + "].second= " <<  centres[count].second << " ; " << endl;
   }
}
