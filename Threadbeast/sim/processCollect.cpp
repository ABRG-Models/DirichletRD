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
using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::runtime_error;
using morph::Gdisplay;
using namespace morph;
using namespace std;
#define NUMPOINTS 21

int main (int argc, char **argv)
{
    if (argc < 2) {
      std::cout << "not enough arguments" << argc << endl;
      return -1;
    }
    //string jsonfile = argv[1]
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
    bool LfixedSeed = atoi(argv[11]); //are we using a fixed seed?
    bool Lgraphics = atoi(argv[12]);
    bool LDn = atoi(argv[13]);
	int numSectors = 12;
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
    cout << " Lcontinue " << Lcontinue << " skipMorph " << skipMorph << endl;
    ofstream afile (logpath + "/centroids.out",ios::app);

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


// initialise DRegion class setting scale
    DRegion M(scale,xspan,logpath);
    cout << "before dissect_boundary " << endl;
    vector<std::pair<double,double>> cGravity;
    cGravity = M.dissectBoundary(); //dissect region boundary
    cout << "Edges size = " << M.edges.size() << endl;
	M.setRadialSegments(); //set the radial segments for regions
    cout << "after first setRadialSegments " << endl;
// include the analysis methods
    Analysis L;

    double rhoInit = 2.5;
    vector<double> fix(3.0, 0.0);
    vector<double> eye(3.0, 0.0);
    vector<double> rot(3.0, 0.0);
    array<float,3> cl_a = morph::Tools::getJetColorF (0.78);
    array<float,3> cl_c = morph::Tools::getJetColorF (0.28);
    array<float,3> cl_b = morph::Tools::getJetColorF (0.58);
    array<float,3> cl_d = morph::Tools::getJetColorF (0.00);
    array<float,3> offset = {{0, 0, 0}};
    float hexWidth = M.Hgrid->hexen.begin()->d/2.0;
    cerr << "d/2: " << hexWidth << endl;// section for solving on the curved boundaries
     int boundaryCount = 0;
     if (Lgraphics) {
        morph::Gdisplay isp(900, 900, 0, 0, "A boundary", rhoInit, 0.0, 0.0);
        cout << "after setting display" << endl;
        isp.resetDisplay (fix, eye, rot);
        cout << "after setting display" << endl;
    //plot stuff here.
        for (auto h : M.Hgrid->hexen) {
            if (M.Creg[h.vi] ==  1 ) {
                isp.drawHex (h.position(), (h.d/2.0f), cl_a);
                boundaryCount++;
            }
            else if (M.Creg[h.vi] > 1) {
                isp.drawHex (h.position(), (h.d/2.0f), cl_c);
            }
        }
        for (int regcnt = 0; regcnt < NUMPOINTS;regcnt++) {
            std::array<float,3> pos = {M.centres[regcnt].first, M.centres[regcnt].second, 0};
            cout << " drawing centre of region x " << M.centres[regcnt].first << " y " << M.centres[regcnt].second << endl;
            isp.drawHex (pos,offset,1.0*hexWidth,cl_a);
            for (unsigned int i=0; i< M.vCoords[regcnt].size();i++) {
                std::array<float,3> pos = {M.vCoords[regcnt][i].first, M.vCoords[regcnt][i].second, 0};
                isp.drawHex (pos,offset,4.0*hexWidth,cl_b);
            }
        }
        cout << "boundaryCount "<<boundaryCount<<endl;
        usleep (10000000); // ten seconds
        isp.redrawDisplay();
        usleep (10000000); // ten seconds
        isp.saveImage(logpath + "/Tesselation.png");
        isp.closeDisplay();
     }

    /*
     * now we integrated over a polygonal tesselation but this time with the equations solved
     * in each region with a separate ksSolver
     */
    vector<ksSolver> S;
    /*
     * Set boundaries based on polygons derived from the Dregion tesselation.
     * stored in curvedBoundary
     */
    M.populateBoundPolygon(1);
    /*
     * Read in the hexGrids for each region
     */
    cout << " read in the HexGrids" << endl;
    for (int j=0; j<NUMPOINTS;j++){
        string str = to_string(j);
        string ename = logpath + "/" + str + ".h5";
        cout << "In the loop populating the ksSolver vector" << endl;
        morph::HexGrid* hG = new morph::HexGrid;
        cout << "In the loop populating the ksSolver vector" << endl;
        hG->load(ename);
        cout << "In the loop populating the ksSolver vector" << endl;
        S.push_back(ksSolver(hG,logpath));
        cout << "In the loop populating the ksSolver vector" << endl;
    }

    /*
     * fill the regions with the hexes from the
     * imported Hgrids and update the boundaries
     */
    for (int j=0;j<NUMPOINTS;j++)
    {
        M.renewRegion(j,S[j].Hgrid->hexen);
    }
    for (int j=0;j<NUMPOINTS;j++)
    {
        M.renewBoundary(j,S[j].Hgrid->hexen);
    }
    afile << "centroids for zero morphin " << endl;
    for (int j=0; j<NUMPOINTS;j++){
        afile << " Region size " << M.regionHex[j].size() << endl;
    }
	cout << "just after populating the regions from the ksSolver vector"<<endl;
    /*
     * Clear and recreate the map of boundary edges as vectors of
     * hexes
     */
    M.edges_clear();
    // swap the radialAngles to the mCoords
    M.swapRadialSegments(false);
    // redissect the boundaries
    cout << "just before renewDissect first time" << endl;
    for (int j=0;j<NUMPOINTS;j++)
    {
        M.renewDissect(j,0);
    }
    cout << "Edges size " << M.edges.size() << endl;

    /*
     * now draw the intial tesselation
     */
    int internalCount = 0;
    boundaryCount = 0;
    if (Lgraphics) {
        morph::Gdisplay idisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
        idisp.resetDisplay (fix, eye, rot);
        for (int j=0;j<NUMPOINTS;j++) {
        // plot stuff here.
            cout << "size of Hexgrid j "<< j << " is " << S[j].Hgrid->num() << endl;
            for (auto h : S[j].Hgrid->hexen) {
                h.y = -h.y;
                if (h.boundaryHex()) {
                    cout << "h.boundaryHex in region " << j <<  " h.vi " << h.vi << endl;
// cout << "region " << j <<  " x " << S[j].Hgrid->d_x[h.vi] << " y " << S[j].Hgrid->d_y[h.vi] << endl;
                    idisp.drawHex (h.position(), (h.d/2.0f), cl_a);
                    boundaryCount++;
                }
                else
                {
                    cout << "h.internal hex in region " << j <<  " h.vi " << h.vi << endl;
//cout << "region " << j <<  " x " << S[j].Hgrid->d_x[h.vi] << " y " << S[j].Hgrid->d_y[h.vi] << endl;
                //idisp.drawHex (h.position(), offset, (h.d/2.0f), cl_b);
                    internalCount++;
                }
            }
            for (unsigned int i=0; i< M.vCoords[j].size();i++) {
                std::array<float,3> pos = {M.vCoords[j][i].first, M.vCoords[j][i].second, 0};
                idisp.drawHex (pos,offset,4.0*hexWidth,cl_c);
                pos = {M.mCoords[j][i].first, M.mCoords[j][i].second, 0};
                idisp.drawHex (pos,offset,4.0*hexWidth,cl_b);
            }

        } // end of loop over regions
        for (auto h : M.Hgrid->hexen) {
            if (M.Creg[h.vi] ==  1 ) {
                idisp.drawHex (h.position(), (h.d/2.0f), cl_b);
                boundaryCount++;
            }
            else if (M.Creg[h.vi] > 1) {
                idisp.drawHex (h.position(), (h.d/2.0f), cl_c);
            }
        }
        usleep (1000000);
        cout << "before redrawDisplay 2 " << endl;
        idisp.redrawDisplay();
        cout << "after redrawDisplay 2" << endl;
        usleep (100000); // one hundred seconds
        idisp.saveImage(logpath + "/Tesselation1.png");
        idisp.closeDisplay();
    } //end of graphics section

    /*
     * read in the fields
     */

    cout << "just after trying to open ../logs/first.h5" << endl;
    for (unsigned int j=0;j<NUMPOINTS;j++) {
        string str = to_string(j);
        string fname = logpath + "/first" + str + ".h5";
        morph::HdfData ginput(fname,1);
		std::string ccstr = "c";
		cout << " j string " << to_string(j) << " length " << ccstr.length()<< endl;
        char * ccst = new char[ccstr.length()+1];
        std::strcpy(ccst,ccstr.c_str());
        std::string nstr = "n";
        char * nst = new char[nstr.length()+1];
        std::strcpy(nst,nstr.c_str());
        cout << "labels "<< nst <<" , " << nstr <<","<< ccst<< "," << ccstr <<endl;
        ginput.read_contained_vals(ccst,S[j].CC);
        ginput.read_contained_vals(nst,S[j].NN);
    }
    cout <<  "just after field creation first morph" << endl;

    for (unsigned int i=0; i< numsteps; i++) {
        for (unsigned int j=0; j<NUMPOINTS;j++) {
            S[j].step(dt, Dn, Dchi, Dc);
        }
    }


	int countHex = 0;
    morph::Gdisplay iidisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
    iidisp.resetDisplay (fix, eye, rot);
   	for (int j = 0;j<NUMPOINTS;j++) {
	 	cout << "in print routine"<<endl;
        unsigned int regsize = S[j].Hgrid->num();
		cout << "in loop over regions " << " size is " << regsize << endl;
        countHex += regsize;
        vector<double> regionNN;
        vector<double> tempNN;
        //for (auto h : S[j].Hgrid->hexen)
        for (unsigned int k=0;k<regsize;k++) {
            tempNN.push_back(S[j].NN[k]);
            cout << "tempNN " << S[j].NN[k] << " k " << k  << endl;
        }
//normalise over the region then write normalised values to normalised NN over hexGrid
        regionNN = L.normalise(tempNN);
        cout << "total number of hexes counted " << countHex << endl;
        int idx = 0;
        for (auto &h : S[j].Hgrid->hexen) {
            array<float,3> colour = morph::Tools::getJetColorF(regionNN[idx]);
            if (!h.boundaryHex()) {
                iidisp.drawHex(h.position(),(h.d/2.0f),colour);
                //    cout << "just after drawHex 2 colour "<< " value " << regionNN[idx] << endl;
             }
             idx++;
         }
         for (unsigned int i=0; i< M.vCoords[j].size();i++) {
             std::array<float,3> pos = {M.vCoords[j][i].first, M.vCoords[j][i].second, 0};
             iidisp.drawHex (pos,offset,4.0*hexWidth,cl_c);
         }
    }//end of loop over regions
    cout << "just before redraw display 1" << endl;
    iidisp.redrawDisplay();
    cout << "just after redraw display 1" << endl;
    cout << "just after to_string"<<endl;
    iidisp.saveImage(logpath + "/nnFieldFull.png");
    usleep (1000000); // one hundred seconds
    cout << "just after saveImage 1" << endl;
    iidisp.closeDisplay();
    /*
    //code run at end of timestepping
    //first save the  ofstream outFile;
    string fname = logpath;
    morph::HdfData fdata(fname);
	for (unsigned int j=0;j<NUMPOINTS;j++)
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
        //data.add_contained_vals("X",M.X[0]);
        //data.add_contained_vals("Y",M.X[1]);
    }
    fdata.add_val ("/Dchi", Dchi);
    fdata.add_val ("/Dn", Dn);
    fdata.add_val ("/Dc",Dc);

    return 0;
    */
};
