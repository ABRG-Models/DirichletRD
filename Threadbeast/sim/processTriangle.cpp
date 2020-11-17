/*
using morph::HexGrid;
using morph::HdfData;:
using morph::Tools;
using morph::Display
*/
#include "region.h"
#include "analysis.h"
#include "ksSolver.h"
#include <morph/display.h>
#include <morph/Config.h>
#include <cctype>
#include <locale>
#include <algorithm>
#include <string>
#include <chrono>
#define numpoints  6
using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::runtime_error;
using morph::Gdisplay;
using namespace morph;
using namespace std;
using namespace std::chrono;

int main (int argc, char **argv)
{
    if (argc < 2) {
      std::cout << "not enough arguments" << argc << endl;
      return -1;
    }
    //string jsonfile = argv[1];
    string logpath = argv[1];
    double dt = stod(argv[2]); //timesetp passed to M.step
    double Dn = stod(argv[3]); //Dn diffusion passed to M.step
    double Dchi = stod(argv[4]); //Dchi chemotaxis passed to M.step
    double Dc = stod(argv[5]);
	int scale = stoi(argv[6]);
	double xspan = stod(argv[7]);
    int numsteps = atoi(argv[8]); //length of integration
    int numprint = atoi(argv[9]); //frequency of printing
    bool Lcontinue = atoi(argv[10]); //logical to determine if coldstart
    unsigned int  LfixedSeed = atoi(argv[11]); //are we using a fixed seed?
    bool Lgraphics = atoi(argv[12]);
    bool LDn = atoi(argv[13]);
	int numSectors = 12;
	double aNoiseGain = 0.1;
	double boundaryFalloffDist = 0.0078;
    vector<morph::BezCurvePath<float>> triangleBound;
    triangleBound.resize(numpoints);
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
    cout << " Lcontinue " << Lcontinue << " skipMorph " << skipMorph << endl;
    ofstream afile (logpath + "/centroids.out",ios::app);

    /*
     * decide if we are using a fixed or
     * run-time dependent seed. Use milliseconds to
     * allow parallel launching

    unsigned int seed;
    milliseconds ms1 = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
    if (LfixedSeed)
        seed = 1;
    else
        seed = static_cast<unsigned int> (ms1.count());
    */
    morph::RandUniform<double> ruf(LfixedSeed);

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
     else {
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
    ofstream gfile ( logpath + "/edges.out");
    ofstream jfile ( logpath + "/results.txt");
	ofstream degfile1 (logpath + "/degree1.data");
	ofstream degfile2 (logpath + "/degree2.data");
	ofstream degfile3 (logpath + "/degree3.data");
    ifstream vfile ("./triangle.inp");
    /*
     * analysis tools
     */
    Analysis L;

    /*
     * now we integrated over a polygonal tesselation but this time with the equations solved
     * in each region with a separate ksSolver
     */
    vector<ksSolver> S;
    /*
     * Set boundaries based on polygons derived from the Dregion tesselation.
     * stored in curvedBoundary
     */
    for (unsigned int i=0; i<numpoints; i++) {
        vector<pair<float,float>> p;
        p.resize(3);
        for (unsigned int j=0; j<3; j++) {
            float x,y;
            vfile >>x ;
            vfile >>y;
            cerr << " x " << x << " y " << y << endl;
            p[j].first = x;
            p[j].second = y;
        }
        cerr << p[0].first << " " << p[0].second << " " << p[1].first << " " << p[1].second << " " << p[2].first << " " << p[2].second << endl;
        morph::BezCurve<float> c1(p[0],p[1]);
        morph::BezCurve<float> c2(p[1],p[2]);
        morph::BezCurve<float> c3(p[2],p[0]);
        triangleBound[i].addCurve(c1);
        triangleBound[i].addCurve(c2);
        triangleBound[i].addCurve(c3);
    }



    for (int j = 0;j<numpoints;j++) {
        pair<double,double> centroid = make_pair(0,0);
        S.push_back(ksSolver(scale, xspan, logpath, triangleBound[j], centroid));
        cout << "in the loop populating the ksVector"<< j <<endl;
    }
    afile << "first setting of centroids" << endl;
    afile << "centroids for zero morphin " << endl;
	cout << "just after populating the regions from the ksSolver vector"<<endl;
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
    if (Lgraphics) {
        morph::Gdisplay idisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
        idisp.resetDisplay (fix, eye, rot);
        for (int j=0;j<numpoints;j++) {
        // plot stuff here.
            cout << "size of Hexgrid j "<< j << " is " << S[j].Hgrid->num() << endl;
            for (auto h : S[j].Hgrid->hexen) {
                if (h.boundaryHex()) {
                    cout << "h.boundaryHex in region " << j <<  " h.vi " << h.vi << endl;
// cout << "region " << j <<  " x " << S[j].Hgrid->d_x[h.vi] << " y " << S[j].Hgrid->d_y[h.vi] << endl;
                    idisp.drawHex (h.position(), (h.d/2.0f), cl_a);
                }
                else
                {
                    cout << "h.internal hex in region " << j <<  " h.vi " << h.vi << endl;
//cout << "region " << j <<  " x " << S[j].Hgrid->d_x[h.vi] << " y " << S[j].Hgrid->d_y[h.vi] << endl;
                //idisp.drawHex (h.position(), offset, (h.d/2.0f), cl_b);
                    internalCount++;
                }
            }

        } // end of loop over regions
        idisp.redrawDisplay();
        cout << "after redrawDisplay 2" << endl;
        usleep (100000); // one hundred seconds
        idisp.saveImage(logpath + "/Tesselation1.png");
        idisp.closeDisplay();
    } //end of graphics section
// initialise the fields
    string fname = logpath + "/first.h5";
    cout<< "just before first data read"<< " lcontinue " << Lcontinue <<endl;
// initialise with random field
    if (Lcontinue)
	{
        morph::HdfData ginput(fname,1);
        cout << "just after trying to open ../logs/first.h5" << endl;
        for (unsigned int j=0;j<numpoints;j++)
        {
		    std::string ccstr = "c" + to_string(j);
		    cout << " j string " << to_string(j) << " length" << ccstr.length()<< endl;
			char * ccst = new char[ccstr.length()+1];
			std::strcpy(ccst,ccstr.c_str());
		    std::string nstr = "n" + to_string(j);
			char * nst = new char[nstr.length()+1];
			std::strcpy(nst,nstr.c_str());
			cout << "labels "<< nst <<" , " << nstr <<","<< ccst<< "," << ccstr <<endl;
	        ginput.read_contained_vals(ccst,S[j].CC);
	        ginput.read_contained_vals(nst,S[j].NN);
		}
	  }
      else
	  {
        for (unsigned int j=0;j<numpoints;j++)
		{
	        for (auto h : S[j].Hgrid->hexen) {
            // boundarySigmoid. Jumps sharply (100, larger is sharper) over length
            // scale 0.05 to 1. So if distance from boundary > 0.05, noise has
            // normal value. Close to boundary, noise is less.
		        double choice = ruf.get();
		        if (choice > 0.5)
				{
                    S[j].NN[h.vi] = - ruf.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = - ruf.get() * aNoiseGain + ccInitialOffset;
				}
		        else
				{
                    S[j].NN[h.vi] = ruf.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = ruf.get() * aNoiseGain + ccInitialOffset;
				}
            if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                double bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary- boundaryFalloffDist)) );
                S[j].NN[h.vi] = (S[j].NN[h.vi] - nnInitialOffset) * bSig + nnInitialOffset;
                S[j].CC[h.vi] = (S[j].CC[h.vi] - ccInitialOffset) * bSig + ccInitialOffset;
		 } //end of if on boundary distance
	    }//end of loop over region
	   }//end of loop over all regions
      } //end of else on Lcontinue
     cout <<  "just after field creation first morph" << endl;

    // begin first time stepping loop
    for (int i=0;i<numsteps;i++)
    {
	    int countHex = 0;
   	    for (int j = 0;j<numpoints;j++) //loop over regions
        {
     	    S[j].step(dt, Dn, Dchi, Dc);
		    if (i%100 == 0)
		    cout << "just after morph 0 squid i "<< j << " NN[5] " << S[j].NN[5] << endl;
	    }
	    if (i%numprint == 0 && Lgraphics)
	    {
          //set up display
            morph::Gdisplay iidisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
            iidisp.resetDisplay (fix, eye, rot);
   	        for (int j = 0;j<numpoints;j++) //loop over regions
	        {
	 	        cout << "in print routine"<<endl;
                unsigned int regsize = S[j].Hgrid->num();
		        cout << "in loop over regions " << " size is " << regsize << endl;
			    countHex += regsize;
		        vector<double> regionNN;
			    vector<double> tempNN;
		      //for (auto h : S[j].Hgrid->hexen)
			    for (unsigned int k=0;k<regsize;k++)
			    {
                    tempNN.push_back(S[j].NN[k]);
                    cout << "tempNN " << S[j].NN[k] << " k " << k  << endl;
                }
//normalise over the region then write normalised values to normalised NN over hexGrid
                regionNN = L.normalise(tempNN);
                cout << "total number of hexes counted " << countHex << endl;
                int idx = 0;
                for (auto &h : S[j].Hgrid->hexen)
                {
                    array<float,3> colour = morph::Tools::getJetColorF(regionNN[idx]);
                    if (!h.boundaryHex()) {
                        iidisp.drawHex(h.position(),(h.d/2.0f),colour);
                    //    cout << "just after drawHex 2 colour "<< " value " << regionNN[idx] << endl;
                    }
                    idx++;
                }
            }
            cout << "just before redraw display 1" << endl;
            iidisp.redrawDisplay();
		    cout << "just after redraw display 1" << endl;
		    cout << "just after to_string"<<endl;
		    iidisp.saveImage(logpath + "/nnField1.png");
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
	for (unsigned int j=0;j<numpoints;j++)
	{
		std::string nstr = "n" + to_string(j);
	    char * nst = new char[nstr.length()+1];
		//std::copy(nstr.begin(),nstr.end(),nst);
		std::strcpy(nst,nstr.c_str());
    	std::string ccstr = "c" + to_string(j);
	    char * ccst = new char[ccstr.length()+1];
	//	std::copy(ccstr.begin(),ccstr.end(),cst);
		std::strcpy(ccst,ccstr.c_str());
		cout << "labels "<< nst <<" , " << nstr <<","<< ccst<< "," << ccstr <<endl;
        fdata.add_contained_vals(ccst,S[j].CC);
        fdata.add_contained_vals(nst,S[j].NN);
    }
    fdata.add_val ("/Dchi", Dchi);
    fdata.add_val ("/Dn", Dn);
    fdata.add_val ("/Dc",Dc);

    return 0;
};
