#include "rd_2d_erm.h"

#include "morph/display.h"
#include <iostream>
#include <vector>
#include <string>

//using namespace std;

int main (int argc, char **argv)
{
    if (argc < 3) {
        cerr << "\nUsage: " << argv[0] << " w0 Dn\n\n";
        cerr << "Be sure to run from the base source directory.\n";
        return -1;
    }
    vector<morph::Gdisplay> displays;
    vector<double> fix(3, 0.0);
    vector<double> eye(3, 0.0);
    vector<double> rot(3, 0.0);

    double rhoInit = 1.4;
    string worldName(argv[1]);
    string winTitle = worldName + ": n";
    displays.push_back (morph::Gdisplay (500, 500, 100, 0, winTitle.c_str(), rhoInit, 0.0, 0.0));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    rhoInit = 1.4;
    winTitle = worldName + ": c";
    displays.push_back (morph::Gdisplay (500, 500, 100, 0, winTitle.c_str(), rhoInit, 0.0, 0.0));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    // How long to run the sim:
    unsigned int maxSteps = 10000000;
    cout << "just before creation of hexgrid" << endl;
    // Instantiate the model object
    RD_2D_Erm M;

    // Modify any parameters before calling M.init()
    M.setLogpath (string("./logs/") + worldName);
    M.Dn = stod(argv[2]);
    M.chi = M.Dn;
    M.Dc = 0.233*M.Dn;
    M.N = 1; // For three chemo attractant molecules (i.e. not using
             // this in the sense of the original Ermentrout system)

    try {
        M.init (displays);
    } catch (const exception& e) {
        cerr << "Exception initialising RD_2D_Karb object: " << e.what() << endl;
    }

    // Start the loop
    bool doing = true;
    while (doing) {

        try {
            // Step the model:
            M.step();
            // Plot every 100 steps:
            if (M.stepCount % 100 == 0) {
                displays[0].resetDisplay (fix, eye, rot);
                M.plot (displays);
            }
            // After a while, stop:
            if (M.stepCount > maxSteps) {
                doing = false;
            }

        } catch (const exception& e) {
            cerr << "Caught exception: " << e.what() << endl;
            doing = false;
        }
    }

    // Before exit, save data
    M.saveState();
    // int numSectors = 12;
    // M.sectorize_radius(numSectors);
    // M.sectorize_angle(numSectors);

     ofstream bfile ("./logs/w0/sectorRadius.out");
     ofstream afile ("./logs/w0/sectorAngle.out");
     int vdegree = 0;
     int mdegree = 0;
     vector <double> rayv;
     vector <double> ringv;

     rayv.resize(0); 
     ringv.resize(0);
     int numSectors = 12;	    
     rayv = M.sectorize_radius(numSectors);
     for (int i=0;i < (int) rayv.size();i++) {
	bfile << " " << rayv[i] << endl;
     }
   
     ringv = M.sectorize_angle(numSectors);
     for (int i=0;i < (int) ringv.size();i++) {
	afile << " " << ringv[i] << endl;
      }
      vdegree = M.find_maxR(rayv);
      mdegree = M.find_maxTh(ringv);
      std::cout<<"  vdegree= "<<vdegree<<endl<<flush;
      std::cout<<"  mdegree= "<<mdegree<<endl<<flush;

    return 0;
}
