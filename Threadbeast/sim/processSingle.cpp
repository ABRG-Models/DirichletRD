#include "region.h"
#include "analysis.h"
#include "ksSolver.h"
#include <morph/display.h>
#include <morph/Config.h>
#include <cctype>
#include <locale>
#include <algorithm>
#include <string>
using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::runtime_error;
using morph::Gdisplay;
using namespace morph;
using namespace std;

int main (int argc, char **argv)
{
    if (argc < 2) {
      std::cout << "not enough arguments" << argc << endl;
      return -1;
    }
    //string jsonfile = argv[1];
    string logpath = argv[1];
    string reg = argv[2];
    double dt = stod(argv[3]); //timesetp passed to M.step
    double Dn = stod(argv[4]); //Dn diffusion passed to M.step
    double Dchi = stod(argv[5]); //Dchi chemotaxis passed to M.step
    double Dc = stod(argv[6]);
	int scale = stoi(argv[7]);
	double xspan = stod(argv[8]);
    int numsteps = atoi(argv[9]); //length of integration
    int numprint = atoi(argv[10]); //frequency of printing
    bool Lcontinue = atoi(argv[11]); //logical to determine if coldstart
    bool LfixedSeed = atoi(argv[12]); //are we using a fixed seed?
    bool Lgraphics = atoi(argv[13]);
    bool LDn = atoi(argv[14]);
	double aNoiseGain = 0.1;
	double boundaryFalloffDist = 0.0078;
    /*
     * open the confgig file and read in the parameters
    morph::Config conf(jsonfile);
    if (!conf.ready) {
        cerr << "Error setting up JSON config: " << conf.emsg << endl;
    }
     */
/*
    double dt = conf.getDouble("dt",0.0001);
    double Dn = conf.getDouble("Dn",5.0);
    double Dchi = conf.getDouble("DChi",5.0);
    double Dc = conf.getDouble("Dc",1.5);
    int scale = conf.getInt("scale",8);
    double xspan = conf.getDouble("xspan",5.0);
    int numsteps = conf.getInt("numsteps",100);
    int numAdjust = conf.getInt("numAdjust",1000000);
    int numprint = conf.getInt("numprint",95);
    string logpath = conf.getString("logpath","../logs");
    //string svgpath = conf.getString("svgpath","./rat.svg");
    double boundaryFalloffDist = conf.getDouble("boundaryFalloffDist",0.0078);
    double aNoiseGain = conf.getDouble("aNoiseGain",0.1);
    int numSectors = conf.getInt("numsectors",12);
    bool Lcontinue = conf.getBool("Lcontinue",false);
*/
    double nnInitialOffset = 1.0;
    double ccInitialOffset = 2.5;
    bool overwrite_logs = true;
    bool skipMorph  = false;
    cout << " Lcontinue " << Lcontinue << " skipMorph " << skipMorph << " reg " <<  reg << endl;
    ofstream afile (logpath + "/centroids.out",ios::app);
    string regname = logpath + "/" + reg + ".h5";
    cout << regname << endl;
// include the analysis methods
    Analysis L;

    unsigned int seed;
    if (LfixedSeed) {
        seed = 1;
    }
    else {
        seed = time(NULL);
    }

    // A ra2yyndo2yym uniform generator returning real/floating point types
    morph::RandUniform<double> ruf(seed);

    /*
     * NOW create a log directory if necessary, and exit on any
     * failures.
     */
    if (morph::Tools::dirExists (logpath) == false) {
        morph::Tools::createDir (logpath);
        if (morph::Tools::dirExists (logpath) == false) {
            cerr << "Failed to create the logpath directory "
                 << logpath << " which does not exist."<< endl;
            return 1;
        }
	}

    ksSolver* S;

    /*
     * create a ksSolver form a HexGrid which is
     * read in from a file
     */
    morph::HexGrid* hexGrid;
    hexGrid = new HexGrid;
    hexGrid->load(regname);
    cout << "after loading " << regname << endl;
    S = new ksSolver(hexGrid, logpath);
// now draw the intial tesselation
    int internalCount = 0;
    double rhoInit = 2.5;
    vector<double> fix(3.0, 0.0);
    vector<double> eye(3.0, 0.0);
    vector<double> rot(3.0, 0.0);
    array<float,3> cl_a = morph::Tools::getJetColorF (0.78);
    array<float,3> cl_c = morph::Tools::getJetColorF (0.28);
    array<float,3> cl_b = morph::Tools::getJetColorF (0.58);
    array<float,3> cl_d = morph::Tools::getJetColorF (0.00);
    array<float,3> offset = {{0, 0, 0}};
    float hexWidth = S->Hgrid->hexen.begin()->d/2.0;
    cerr << "d/2: " << hexWidth << endl;// section for solving on the curved boundaries
    int boundaryCount = 0;
    if (Lgraphics) {
        morph::Gdisplay idisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
        idisp.resetDisplay (fix, eye, rot);
        // plot stuff here.
        cout << "size of Hexgrid  " << " is " << S->Hgrid->num() << endl;
            for (auto h : S->Hgrid->hexen) {
                if (h.boundaryHex()) {
                    cout << "h.boundaryHex in region " <<  " h.vi " << h.vi << endl;
                    idisp.drawHex (h.position(), (h.d/2.0f), cl_a);
                    boundaryCount++;
                }
                else {
                    cout << "h.internal hex in region " <<  " h.vi " << h.vi << endl;
                //idisp.drawHex (h.position(), offset, (h.d/2.0f), cl_b);
                    internalCount++;
                }
            } // end of loop over hexen
        usleep (1000000);
        cout << "before redrawDisplay 2 " << endl;
        idisp.redrawDisplay();
        cout << "after redrawDisplay 2" << endl;
        usleep (100000); // one hundred seconds
        idisp.saveImage(logpath + "/Tesselation1" + reg + ".png");
        idisp.closeDisplay();
    } //end of graphics section
// initialise the fields
    string fname = logpath + "/first" + reg + ".h5";
    cout<< "just before first data read"<< " lcontinue " << Lcontinue << fname << endl;
// initialise with random field
    if (Lcontinue) {
        morph::HdfData ginput(fname,1);
        cout << "just after trying to open ../logs/first.h5" << endl;
		std::string ccstr = "c" ;
        char * ccst = new char[ccstr.length()+1];
        std::strcpy(ccst,ccstr.c_str());
        std::string nstr = "n" ;
        char * nst = new char[nstr.length()+1];
        std::strcpy(nst,nstr.c_str());
        cout << "labels "<< nst <<" , " << nstr <<","<< ccst<< "," << ccstr <<endl;
        ginput.read_contained_vals(ccst,S->CC);
	    ginput.read_contained_vals(nst,S->NN);
		}
    else {
	    for (auto h : S->Hgrid->hexen) {
        // boundarySigmoid. Jumps sharply (100, larger is sharper) over length
        // scale 0.05 to 1. So if distance from boundary > 0.05, noise has
        // normal value. Close to boundary, noise is less.
            double choice = ruf.get();
            if (choice > 0.5) {
                S->NN[h.vi] = - ruf.get() * aNoiseGain +nnInitialOffset;
                S->CC[h.vi] = - ruf.get() * aNoiseGain + ccInitialOffset;
            }
            else {
                S->NN[h.vi] = ruf.get() * aNoiseGain +nnInitialOffset;
                S->CC[h.vi] = ruf.get() * aNoiseGain + ccInitialOffset;
            }
            if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                double bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary- boundaryFalloffDist)) );
                S->NN[h.vi] = (S->NN[h.vi] - nnInitialOffset) * bSig + nnInitialOffset;
                S->CC[h.vi] = (S->CC[h.vi] - ccInitialOffset) * bSig + ccInitialOffset;
            } //end of if on boundary distance
        }//end of loop over region
    } //end of else on Lcontinue
     cout <<  "just after field creation first morph" << endl;

    // begin first time stepping loop
    for (int i=0;i<numsteps;i++) {
	    int countHex = 0;
     	S->step(dt, Dn, Dchi, Dc);
	    if (i%numprint == 0 && Lgraphics) {
          //set up display
            morph::Gdisplay iidisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
            iidisp.resetDisplay (fix, eye, rot);
	 	    cout << "in print routine"<<endl;
            unsigned int regsize = S->Hgrid->num();
		    cout << "in loop over regions " << " size is " << regsize << endl;
			countHex += regsize;
		    vector<double> regionNN;
			vector<double> tempNN;
			for (unsigned int k=0;k<regsize;k++) {
                tempNN.push_back(S->NN[k]);
                cout << "tempNN " << S->NN[k] << " k " << k  << endl;
            }
//normalise over the region then write normalised values to normalised NN over hexGrid
            regionNN = L.normalise(tempNN);
            cout << "total number of hexes counted " << countHex << endl;
            int idx = 0;
            for (auto &h : S->Hgrid->hexen) {
                array<float,3> colour = morph::Tools::getJetColorF(regionNN[idx]);
                if (!h.boundaryHex()) {
                    iidisp.drawHex(h.position(),(h.d/2.0f),colour);
                    //    cout << "just after drawHex 2 colour "<< " value " << regionNN[idx] << endl;
                }
                idx++;
            }
            cout << "just before redraw display 1" << endl;
            iidisp.redrawDisplay();
		    cout << "just after redraw display 1" << endl;
		    cout << "just after to_string"<<endl;
		    iidisp.saveImage(logpath + "/nn" + reg + ".png");
            usleep (1000000); // one hundred seconds
		    cout << "just after saveImage 1" << endl;
		    iidisp.closeDisplay();
		    cout << "just after close display 1 i " << i << endl;
        }//end of loop on numprint drawing fields
		 //cout << "iteration " << i << " of morph 1 time step" << endl;
    } // end of second time stepping

    //code run at end of timestepping
    //first save the  ofstream outFile;
    morph::HdfData fdata(fname);
	std::string nstr = "n";
	char * nst = new char[nstr.length()+1];
	std::strcpy(nst,nstr.c_str());
    std::string ccstr = "c";
	char * ccst = new char[ccstr.length()+1];
	std::strcpy(ccst,ccstr.c_str());
	cout << "labels "<< nst <<" , " << nstr <<","<< ccst<< "," << ccstr <<endl;
    fdata.add_contained_vals(ccst,S->CC);
    cout << "after add_contained vals CC" <<endl;
    fdata.add_contained_vals(nst,S->NN);
    cout << "after add_contained vals NN`" <<endl;
    //fdata.add_val ("/Dchi", Dchi);
    //fdata.add_val ("/Dn", Dn);
    //fdata.add_val ("/Dc",Dc);

    return 0;
};
