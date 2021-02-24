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

int main (int argc, char **argv)
{
    if (argc < 2) {
      std::cout << "not enough arguments" << argc << endl;
      return -1;
    }
    //string jsonfile = argv[1];
    string logpath = argv[1];
	int scale = stoi(argv[2]);
	double xspan = stod(argv[3]);
    bool LfixedSeed = atoi(argv[4]); //are we using a fixed seed?
    bool lRandomOrientation = atoi(argv[5]);
    int numSectors = atoi(argv[6]);
    int radialOrder = atoi(argv[7]);
    int angularOrder = atoi(argv[8]);
    double Dn = stod(argv[9]);
    bool Lgraphics = true;
    bool overwrite_logs = true;
    bool skipMorph  = true;
    const int NUMPOINTS = 41;
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
    DRegion M(scale,xspan,logpath,NUMPOINTS);
    M.setCreg();
    cout << "before dissect_boundary " << endl;
//	for (auto h : M.Hgrid->hexen) {
//	   cout << " first h.vi output " << h.vi << endl;
//	   }
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
            pos = {M.vCoords[regcnt][0].first, M.vCoords[regcnt][0].second, 0};
            isp.drawHex (pos,offset,4.0*hexWidth,cl_b);
            /*
            for (unsigned int i=0; i< M.vCoords[regcnt].size();i++) {
                std::array<float,3> pos = {M.vCoords[regcnt][i].first, M.vCoords[regcnt][i].second, 0};
                isp.drawHex (pos,offset,4.0*hexWidth,cl_b);
            }
            */
        }
        cout << "boundaryCount "<<boundaryCount<<endl;
        usleep (10000000); // ten seconds
        isp.redrawDisplay();
        usleep (10000000); // ten seconds
        isp.saveImage(logpath + "/Tesselation.png");
        isp.closeDisplay();
    }
    cerr << "d/2: " << hexWidth << endl;
    if (Lgraphics) {
        morph::Gdisplay iidisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
        iidisp.resetDisplay (fix, eye, rot);
   	    for (int j = 0;j<NUMPOINTS;j++) { //loop over regions
            double phase = M.radialAngles[j][0];
            cout << " Bessel Phase " << phase << endl;
     	    int numHexes = M.regionBessel(j, radialOrder, angularOrder, phase, Dn);
            cout << "after regionBessel numHexes" << numHexes <<  endl;
          //set up display
	 	    cout << "in print routine NN vector length "<< M.NN.size() << endl;
            unsigned int regsize = M.regionHex[j].size();
            unsigned int NNsize = M.NN[j].size();
		    cout << "loop over region " << j << " region size " << regsize << " NN size " << NNsize <<  endl;
		    vector<double> regionNN;
			vector<double> tempNN;
			for (unsigned int k=0;k<regsize;k++)
			{
                tempNN.push_back(M.NN[j][k]);
                cout << "tempNN " << M.NN[j][k] << " k " << k  << endl;
            }
//normalise over the region then write normalised values to normalised NN over hexGrid
            regionNN = L.normalise(tempNN);
            int idx = 0;
            for (auto &h : M.regionHex[j])
            {
                array<float,3> colour = morph::Tools::getJetColorF(regionNN[idx]);
                if (M.Creg[h.vi] == 0) {
                    iidisp.drawHex(h.position(),(h.d/2.0f),colour);
                    //    cout << "just after drawHex 2 colour "<< " value " << regionNN[idx] << endl;
                    }
                idx++;
                if (M.Creg[h.vi] == 1) {
                    iidisp.drawHex (h.position(), (h.d / 2.0 ), cl_d);
                }
            }
            std::array<float,3> pos = {M.vCoords[j][0].first, M.vCoords[j][0].second, 0};
            iidisp.drawHex (pos,offset,4.0*hexWidth,cl_c);
            /*
            for (unsigned int i=0; i< M.vCoords[j].size();i++) {
                std::array<float,3> pos = {M.vCoords[j][i].first, M.vCoords[j][i].second, 0};
                iidisp.drawHex (pos,offset,4.0*hexWidth,cl_c);
            }
            */
        }//end of loop over regions
        cout << "just before redraw display 1" << endl;
        iidisp.redrawDisplay();
		cout << "just after redraw display 1" << endl;
		cout << "just after to_string"<<endl;
		iidisp.saveImage(logpath + "/nnField1.png");
        usleep (1000000); // one hundred seconds
		cout << "just after saveImage 1" << endl;
		iidisp.closeDisplay();
    }

    cout << " just after writing data "  << endl;
    for (int j=0; j<NUMPOINTS; j++) {
        M.renewBoundary(j, M.regionHex[j], true);
    }
    for (int j=0; j<NUMPOINTS; j++) {
        M.renewDissect(j, 0);
    }
    gfile << endl << "analysis on first morphing iteration " << endl;
    vector <int> radiusDVector;
    vector <int> angleDVector;
	vector <double> angleVector;
    vector <double> radiusVector;
    int degreeRadius = 0;
    int degreeAngle = 0;
	double avDegreeAngle = 0.0;
	double avDegreeRadius = 0.0;
	double occupancy = 0.0;
    int countRegions = 0;
    double tempArea = 0.0;
    double tempPerimeter = 0.0;
    double avAbsCorrelation = 0.0;
    int angleOffset = 0;
    int radiusOffset = 0;
    for (int j=0;j<NUMPOINTS;j++) {
        if (M.regArea(j) != 0){
			int regionCount = 0;
            gfile<<"in the degree loop" << endl;
            //angle degree
            tempArea = M.regArea(j);
            tempPerimeter = M.renewRegPerimeter(j);

            // digital version
            angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, M.NN[j]);
            degreeAngle = L.find_zeroDAngle(angleDVector);
            gfile << "region "<< j << " degreeDAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;

            // analogue version
            angleVector = M.sectorize_reg_angle(j,numSectors,radiusOffset, numSectors, M.NN[j]);
            degreeAngle = L.find_zeroAngle(angleVector,3);
            gfile << "region "<< j << " degreeAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;

            //radial degree
            degreeRadius = -100;
            radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors, M.NN[j]);
            degreeRadius = L.find_zeroDRadius(radiusDVector);
            gfile  << "region "<< j << " degreeDRadius "<< degreeRadius << "  " <<endl ;

            //radial degree
            degreeRadius = -100;
			radiusVector = M.sectorize_reg_radius(j,numSectors, angleOffset, angleOffset + numSectors, M.NN[j]);
			degreeRadius = L.find_zeroRadius(radiusVector,3);
            gfile <<  " region "<< j << " degreeRadius  "<< degreeRadius << "  " <<endl << endl;
        } //end of if on non-zero regions
    } //end of loop on NUMPOINT


	      avDegreeAngle = 0;
	      avDegreeRadius = 0;
	      occupancy = 0;
	      countRegions = 0;
		  tempArea = 0;
		  tempPerimeter = 0;
		  avAbsCorrelation = 0;
          int max_comp = 5 * NUMPOINTS;
          M.Srandom_correlate(max_comp, 0, true, true);
          for (int j=0;j<NUMPOINTS;j++) {
	        if (M.regArea(j) != 0)
			{
	          countRegions++;
		      occupancy += M.regNNfrac(j);
              tempArea = M.regArea(j);
              tempPerimeter = M.regPerimeter(j);
		      avAbsCorrelation += M.Srenewcorrelate_edges(j,0);
              angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, M.NN[j]);
              degreeAngle = L.find_zeroDAngle(angleDVector);
		      avDegreeAngle += degreeAngle;
              //radial degree
		      degreeRadius = 0;
              radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors, M.NN[j]);
              degreeRadius = L.find_zeroDRadius(radiusDVector);
              avDegreeRadius += degreeRadius;

              degfile1 << degreeAngle/2 << " " << degreeRadius << " " << M.regNNfrac(j) << " " << tempArea << " "<< tempPerimeter<<endl<<flush;
	       } //end of if on non-zero regions
	    } //end of loop on NUMPOINTs
        avDegreeAngle = avDegreeAngle / (1.0 * countRegions);
        avDegreeRadius = avDegreeRadius / (1.0 * countRegions);
        avAbsCorrelation = avAbsCorrelation / (1.0 * countRegions);
	    occupancy = occupancy / (1.0 * countRegions);
        cout << "Regions counted in first analysis loop " << countRegions << endl;
	    // jfile << avDegreeAngle <<" "<<avDegreeRadius<<endl;
	    jfile << " "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation << endl;
//end of integration after solving on polygonal regions via ksSolver

    return 0;
};
