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

/*
using morph::HexGrid;
using morph::HdfData;:
using morph::Tools;
using morph::Display
*/
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
    double dt = stod(argv[2]); //timesetp passed to M.step
    double Dn = stod(argv[3]); //Dn diffusion passed to M.step
    double Dchi = stod(argv[4]); //Dchi chemotaxis passed to M.step
    double Dc = stod(argv[5]);
	int scale = stoi(argv[6]);
	double xspan = stod(argv[7]);
    int numsteps = atoi(argv[8]); //length of integration
    int numprint = atoi(argv[9]); //frequency of printing
    bool Lcontinue = atoi(argv[10]); //true if readdata false if coldstart
    bool lminRadius = atoi(argv[11]); //true for min circles, false for max

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
    bool skipMorph  = true;
    cout << " Lcontinue " << Lcontinue << " skipMorph " << skipMorph << endl;


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
	for (auto h : M.Hgrid->hexen) {
	   cout << " first h.vi output " << h.vi << endl;
	   }
    vector<std::pair<double,double>> centroids;
    centroids = M.dissectBoundary(); //dissect region boundary
    M.setRadialSegments();
    double hexWidth = M.Hgrid->hexen.begin()->d/2.0;
// include the analysis methods
    Analysis L;

    string fname = logpath + "/first.h5";
    cout<< "just before first data read"<< " Lcontinue " << Lcontinue <<endl;
    unsigned int seed = time(NULL);
    // A rando2yym uniform generator returning real/floating point types
    morph::RandUniform<double> ruf(seed);
// initialise with random field
    if (Lcontinue) {
	    morph::HdfData input (fname,1);
	    cout<< "just after first data read fname "<< fname << endl;
	    input.read_contained_vals("n",M.NN);
	    input.read_contained_vals("c",M.CC);
	    cout<< "just after input of NN and CC1"<< endl;
//	    input.close();
    }
    else {
		for (auto h : M.Hgrid->hexen) {
		    double choice = morph::Tools::randDouble();
            // boundarySigmoid. Jumps sharply (100, larger is sharper) over length
            // scale 0.05 to 1. So if distance from boundary > 0.05, noise has
            // normal value. Close to boundary, noise is less.
		    if (choice > 0.5)
			{
                M.NN[h.vi] = - ruf.get() * aNoiseGain +nnInitialOffset;
                M.CC[h.vi] = - ruf.get() * aNoiseGain + ccInitialOffset;
			}
			else
			{
                M.NN[h.vi] = - ruf.get() * aNoiseGain +nnInitialOffset;
                M.CC[h.vi] = - ruf.get() * aNoiseGain + ccInitialOffset;
			}
            if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                double bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary- boundaryFalloffDist)) );
                M.NN[h.vi] = (M.NN[h.vi] - nnInitialOffset) * bSig + nnInitialOffset;
              //  M.CC[h.vi] = (M.CC[h.vi] - ccInitialOffset) * bSig + ccInitialOffset;
		    } //end of if on boundary distance
	    }//end of loop over region
    } //end of else on Lcontinue
    cout <<  "just after field creation" << endl;

     vector<double> fix(3.0, 0.0);
     vector<double> eye(3.0, 0.0);
     vector<double> rot(3.0, 0.0);
     double rhoInit = 2.5;
     //double rhoInit = 5.0;
     array<float,3> cl_a = morph::Tools::getJetColorF (0.78);
     array<float,3> cl_c = morph::Tools::getJetColorF (0.28);
     array<float,3> cl_b = morph::Tools::getJetColorF (0.58);
     array<float,3> cl_d = morph::Tools::getJetColorF (0.00);
     array<float,3> offset = {{0, 0, 0}};
     morph::Gdisplay isp(900, 900, 0, 0, "A boundary", rhoInit, 0.0, 0.0);
     cout << "after setting display" << endl;
     isp.resetDisplay (fix, eye, rot);
    cout << "after setting display" << endl;
    for (int j=0; j<NUMPOINTS;j++) {
        vector<hexGeometry::lineSegment> regionLineSegments;
        regionLineSegments = M.polygonSides(j, true);
        unsigned int sideSize = regionLineSegments.size();
        cout << "sideSize " << sideSize << " region " << j << endl;
        for (auto& h : M.Hgrid->hexen) {
            for (unsigned int i=0; i<sideSize; i++) {
                hexGeometry::point p;
                //p.first = M.Hgrid->d_x[h.vi];
                p.first = h.position()[0];
                p.second = h.position()[1];
               // p.second = M.Hgrid->d_y[h.vi];
               // if (M.hGeo->hexIntersectLineSegment(regionLineSegments[i], p, hexWidth)) {
               if (M.hGeo->hexIntersectLineSegment(regionLineSegments[i],h)){
                        h.setFlags(HEX_IS_REGION_BOUNDARY);
                }
            }
        }
    }

    //plot stuff here.
     int boundaryCount = 0;
     for (auto h : M.Hgrid->hexen) {
         //if (M.Creg[h.vi] ==  1 ) {
         if (h.testFlags(HEX_IS_REGION_BOUNDARY) == true) {
             isp.drawHex (h.position(), (h.d/2.0f), cl_a);
             boundaryCount++;
         }
         else if (M.Creg[h.vi] > 1) {
              isp.drawHex (h.position(), (h.d/2.0f), cl_c);
         }
         else if (M.Creg[h.vi] == 1) {
              isp.drawHex (h.position(), (h.d/2.0f), cl_b);
         }


       //  else
        // {
          //   isp.drawHex (h.position(), offset, (h.d/2.0f), cl_b);
         //}
     }
     for (int regcnt = 0; regcnt < NUMPOINTS;regcnt++) {
         std::array<float,3> pos = {M.centres[regcnt].first, M.centres[regcnt].second, 0};
         cout << " drawing centre of region x " << M.centres[regcnt].first << " y " << M.centres[regcnt].second << endl;
         isp.drawHex (pos,offset,4.0*hexWidth,cl_a);
     }
     cout << "boundaryCount "<<boundaryCount<<endl;
     usleep (10000000); // ten seconds
     isp.redrawDisplay();
     usleep (10000000); // ten seconds
     isp.saveImage(logpath + "/Tesselation.png");
     isp.closeDisplay();

     cerr << "d/2: " << hexWidth << endl;
     cout << "before time-stepping loop" << endl;
     for (int i=0;i<numsteps;i++) {
//cout << " just before time step " << " i " << i << endl;
         M.step(dt, Dn, Dchi, Dc);
		 if (i%numprint == 0) {
             morph::Gdisplay disp(900, 900, 0, 0, "first run", rhoInit, 0.0, 0.0);
             disp.resetDisplay (fix, eye, rot);
		     cout << "in print routine"<<endl;
		     vector<double> normalNN;
		     normalNN.resize(M.n);
		     int countHex = 0;
		     for (int j=0;j<NUMPOINTS;j++) {
                 unsigned int regsize = M.regionIndex[j].size();
		         cout << "in loop over regions " << j << " size is " << regsize << endl;
			     countHex += regsize;
		         vector<double> regionNN(regsize);
			     vector<double> tempNN(regsize);
			     vector<int> regionIdx(regsize);
			     int k = 0;
		         for (auto h : M.regionIndex[j]){
			         int index = h.vi;
			         tempNN[k] = M.NN[index];
			         regionIdx[k] = index;
			         k++;
			     }
//cout << "normalise over the region then write  to normalised NN over hexGrid" << endl;
                 regionNN = L.normalise(tempNN);
		         for (unsigned int k=0;k<regsize;k++) {
		             normalNN[regionIdx[k]] = regionNN[k];
//  cout << M.NN[regionIdx[k]] << "    " << normalNN[regionIdx[k]] << " " << regionIdx[k] <<endl;
//afile << regionIdx[k] <<endl;
			     }
             std::array<float,3> pos = {M.centres[j].first, M.centres[j].second, 0};
             isp.drawHex (pos,offset,4.0*hexWidth,cl_d);
		     } //end of loop over regions
		     cout << "total number of hexes counted " << countHex << endl;
		     for (auto h : M.Hgrid->hexen) {
		        // cout << "just before drawHex  "<< h.vi << "normalNN " << normalNN[h.vi] << endl;
		         if (M.Creg[h.vi] == 0) {
			          if (h.d/2.0 > 0.05){
			              cerr << "hex too large for " << h.vi <<endl;
					  }
			          array<float,3> colour = morph::Tools::getJetColorF(normalNN[h.vi]);
                  //    cout << "just before  drawHex in Tesellation"<<endl;
	                  disp.drawHex(h.position(),hexWidth,colour);
                  //    cout << "just after drawHex in Tessleation"<<endl;
 			     }
		         else{
			         disp.drawHex(h.position(),hexWidth,cl_c);
			         }
		          }
		          cout << "just before redraw display 0" << endl;
                  usleep (1000000);
                  disp.redrawDisplay();
                  usleep (1000000); // one hundred seconds
		          disp.saveImage(logpath + "/nnField.png");
                  disp.closeDisplay();
             } // end of print on numprint
         } //end of numsteps loop
//cout << " just after time step i = " << i << endl;

//code run at end of timestepping
//first save the  ofstream outFile;
         morph::HdfData data(fname);
         data.add_contained_vals("c",M.CC); //chemoattactant
         data.add_contained_vals("n",M.NN); //neuronal density
         data.add_contained_vals("x",M.Hgrid->d_x); //x coord of hex centres
         data.add_contained_vals("y",M.Hgrid->d_y); //y coord of hex centres
         data.add_val ("/Dchi", Dchi); //Chi parameter
         data.add_val ("/Dn", Dn); //Dn paramater
         data.add_val ("/Dc",Dc); //Dc parameter

         cout << " just after writing data "  << endl;
// post run analysis
// look at correlation between adjacent edges
         cout << "before correlate_edges" << endl;
	     double avAbsCorrelation = 0;
         avAbsCorrelation = M.correlate_edges();
         cout << "after correlate_edges" << endl;
// look at correlation between random edges
         const int max_comp = NUMPOINTS*4;
         cout << "before random_correlate_edges" << endl;
         M.random_correlate(max_comp,0);
         cout << "after random_correlate_edges" << endl;
//  cout<<"after correlate_edges" << endl;
         vector <int> radiusDVector;
         vector <int> angleDVector;
	     vector <double> angleVector;
         vector <double> radiusVector;
         int degreeRadius;
         int degreeAngle;
	     double tempArea = 0;
	     double tempPerimeter = 0;
		 int angleOffset = 0;
		 int radiusOffset = 0;
         for (int j=0;j<NUMPOINTS;j++) {
             if (M.regArea(j) != 0){
			     int regionCount = 0;
                 gfile<<"in the degree loop" << endl;
                 tempArea = M.regArea(j); //area by hexes
                 tempPerimeter = M.regPerimeter(j); //perimeter by hexes

                  // digital version
                 angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors);
                 degreeAngle = L.find_zeroDAngle(angleDVector);
                 gfile << "region "<< j << " degreeDAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;

                  // analogue version
                 angleVector = M.sectorize_reg_angle(j,numSectors,radiusOffset, numSectors);
                 angleVector = L.meanzero_vector(angleVector);
                  //degreeAngle = M.find_max(angleVector,3);
                 degreeAngle = L.find_zeroAngle(angleVector,3);
                 gfile << "region "<< j << " degreeAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;
                  //radial degree
                 radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2);
                 degreeRadius = L.find_zeroDRadius(radiusDVector);
                 gfile  << "region "<< j << " degreeDRadius "<< degreeRadius << "  " <<endl ;

                  ///radial degree
                 degreeRadius = -100;
                 int newdegreeRadius = 0;
                 for (int angleOffset=0; angleOffset<numSectors -1; angleOffset += 3) {
			          radiusVector = M.sectorize_reg_radius(j,numSectors, angleOffset, angleOffset + numSectors/2);
			          newdegreeRadius = L.find_zeroRadius(radiusVector,3);
			          if (newdegreeRadius > degreeRadius)
				          degreeRadius = newdegreeRadius;
		         }


                 gfile <<  " region "<< j << " degreeRadius  "<< degreeRadius << "  " <<endl << endl;


                 regionCount++;

	     } //end of if on non-zero regions
      } //end of loop on NUMPOINTs
//writing out of the image file

	  double avDegreeAngle = 0;
	  double avDegreeRadius = 0;
	  double occupancy = 0;
      tempArea = 0;
	  tempPerimeter = 0;
	  int countRegions = 0;
      for (int j=0;j<NUMPOINTS;j++) {
	      if (M.regArea(j) != 0){
	          countRegions++;
		      occupancy += M.regNNfrac(j);
              tempArea = M.regArea(j);
              tempPerimeter = M.regPerimeter(j);

              angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors);
              degreeAngle = L.find_zeroDAngle(angleDVector);
		      avDegreeAngle += degreeAngle;
              //radial degree
		      degreeRadius = 0;
              radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2);
              degreeRadius = L.find_zeroDRadius(radiusDVector);
              avDegreeRadius += degreeRadius;

               degfile1 << degreeAngle/2 << " " << degreeRadius << " " << M.regNNfrac(j) << " " << tempArea<< " "<< tempPerimeter<< " xval " << M.centres[j].first << " yval " << M.centres[j].second <<endl<<flush;
		  }
	  } //end of loop on NUMPOINTs
      avDegreeAngle = avDegreeAngle / (1.0 * countRegions);
      avDegreeRadius = avDegreeRadius / (1.0 * countRegions);
	  occupancy = occupancy / (1.0 * countRegions);
	  jfile <<Dn<<" "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation<<endl;
//end of integration on the original tesselation
// if (skipMorph) return 0;
    cout << "just before setting curved boundaries" <<endl;
// section for solving on the curved boundaries
    vector<ksSolver> S;
// now set the boundaries for the regions, stored in curvedBoundary
	M.setRadialSegments();
    for (int j=0; j<NUMPOINTS; j++) {
       cout << " after setRadialSegments vCoords " << M.vCoords[j].size() << " mCoords " << M.mCoords[j].size() << endl;
    }
	for (int j = 0;j<NUMPOINTS;j++) {
        std::pair<float, float> centroid =  M.baryCentre(j);
        float radius;
        if (lminRadius) {
            radius = static_cast<float>(M.min_radius(j,true));
        }
        else {
            radius = static_cast<float>(M.max_radius(j,true));
        }

        cout << " radius = " << radius << endl;
        S.push_back(ksSolver(scale, xspan, logpath, radius, centroid));
        cout << "in the loop populating the ksVector"<< j <<endl;

    }
	     cout << "just after populating the ksVector"<<endl;
// test the vertices are ordered correctly by increasing angle from the centroid
    for (int j=0; j<NUMPOINTS;j++) {
        M.testRegionVertices(j);
    }
// now draw the intial tesselation
    morph::Gdisplay mdisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
    mdisp.resetDisplay (fix, eye, rot);
    hexWidth = M.Hgrid->hexen.begin()->d/2.0;
    for (int regcnt = 0; regcnt < NUMPOINTS;regcnt++) {
        std::array<float,3> pos = {M.centres[regcnt].first, M.centres[regcnt].second, 0};
        std::array<float,3> pos1 = {M.baryCentre(regcnt).first,M.baryCentre(regcnt).second,0};
        cout << " drawing centre of region x " << M.centres[regcnt].first << " y " << M.centres[regcnt].second << endl;
        isp.drawHex (pos,offset,4.0*hexWidth,cl_a);
        isp.drawHex (pos1,offset,4.0*hexWidth,cl_b);
    }
    for (int j=0;j<NUMPOINTS;j++) {
        if (M.regArea(j) == 0) break;
// plot stuff here.
        int boundaryCount = 0;
        int internalCount = 0;
        cout << "size of Hexgrid j "<< j << " is " << S[j].Hgrid->num() << endl;
        for (auto h : S[j].Hgrid->hexen) {
            if (h.boundaryHex()) {
                mdisp.drawHex (h.position(), (h.d/2.0f), cl_a);
    			boundaryCount++;
            }
        }
    } //end of for loop on j
    usleep (1000000);
    cout << "before redrawDisplay 2 " << endl;
    mdisp.redrawDisplay();
    cout << "after redrawDisplay 2" << endl;
    usleep (100000); // one hundred seconds
    mdisp.saveImage(logpath + "/Tesselation2.png");
    mdisp.closeDisplay();
// initialise the fields
    string gname = logpath + "/second.h5";
    cout<< "just before second data read"<< " lcontinue " << Lcontinue <<endl;
// initialise with random field
    if (Lcontinue)
	{
        morph::HdfData ginput(gname,1);
        cout << "just after trying to open ../logs/second.h5" << endl;
        for (unsigned int j=0;j<NUMPOINTS;j++)
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
        for (unsigned int j=0;j<NUMPOINTS;j++)
		{
            for (auto h : S[j].Hgrid->hexen) {
            // boundarySigmoid. Jumps sharply (100, larger is sharper) over length
            // scale 0.05 to 1. So if distance from boundary > 0.05, noise has
            // normal value. Close to boundary, noise is less.
                double choice = morph::Tools::randDouble();
		        if (choice > 0.5)
				{
                    S[j].NN[h.vi] = - ruf.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = - ruf.get() * aNoiseGain + ccInitialOffset;
				}
		        else
				{
                    S[j].NN[h.vi] = - ruf.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = - ruf.get() * aNoiseGain + ccInitialOffset;
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
    vector<vector<hexGeometry::point>> tickPoints;
    tickPoints.resize(NUMPOINTS);
    int ticks = 100;
    for (int j=0;j<NUMPOINTS;j++) {
        tickPoints[j] = M.divideRegionBoundary(j, ticks);
    }

    cout <<  "just after tickpoints first morph" << endl;
// begin second time stepping loop
    for (int i=0;i<numsteps;i++)
    {
        for (int j = 0;j<NUMPOINTS;j++) //loop over regions
        {
            S[j].step(dt, Dn, Dchi, Dc);
        }
    }

    //set up display
    morph::Gdisplay mmdisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
    mmdisp.resetDisplay (fix, eye, rot);
    for (int j = 0;j<NUMPOINTS;j++) //loop over regions
        {
            cout << "in print routine"<<endl;
            unsigned int regsize = S[j].Hgrid->num();
            cout << "in loop over regions " << " size is " << regsize << endl;
            vector<double> tempNN;
            vector<double> regionNN;
			for (auto h : S[j].Hgrid->hexen) {
                tempNN.push_back(S[j].NN[h.vi]);
                cout << "tempNN " << S[j].NN[h.vi] << " h.vi " << h.vi  << endl;
			}
//normalise over the region then write normalised values to normalised NN over hexGrid
            regionNN = L.normalise(tempNN);
            int idx = 0;
		    for (auto &h : S[j].Hgrid->hexen)
			{
                std::pair<double, double> inPoint;
                inPoint.first = h.x;
                inPoint.second = h.y;
                cout << " inRegion " << " xval " << inPoint.first << " yval " << inPoint.second <<endl;
                if (lminRadius) {
                    array<float,3> colour = morph::Tools::getJetColorF(regionNN[idx]);
                    mmdisp.drawHex(h.position(),(h.d/2.0f),colour);
                    cout << "just after drawHex 2  value " << regionNN[idx] << endl;
                }
                else {
                    //if (M.inRegion(j, inPoint, tickPoints[j], 0.01)){
                        if (1) {
                        array<float,3> colour = morph::Tools::getJetColorF(regionNN[idx]);
                        mmdisp.drawHex(h.position(),(h.d/2.0f),colour);
                        cout << "just after drawHex 2  value " << regionNN[idx] << endl;
                    }
                }
                idx++;
            }
        }//end of loop over regions
    cout << "just before redraw display 1" << endl;
    mmdisp.redrawDisplay();
    cout << "just after redraw display 1" << endl;
    cout << "just after to_string"<<endl;
    mmdisp.saveImage(logpath + "/nnField2.png");
    usleep (1000000); // one hundred seconds
    cout << "just after saveImage 1" << endl;
    mmdisp.closeDisplay();

    //code run at end of timestepping
    //first save the  ofstream outFile;
    morph::HdfData gdata(gname);
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
        gdata.add_contained_vals(ccst,S[j].CC);
        gdata.add_contained_vals(nst,S[j].NN);
        //data.add_contained_vals("X",M.X[0]);
        //data.add_contained_vals("Y",M.X[1]);
     }
     gdata.add_val ("/Dchi", Dchi);
     gdata.add_val ("/Dn", Dn);
     gdata.add_val ("/Dc",Dc);

     cout << " just after writing data "  << endl;
	  // repopulate the regions
	 for (int j=0;j<NUMPOINTS;j++)
     {
	    M.renewRegion(j,S[j].Hgrid->hexen);
	 }
     cout << " just after repopulating regions "  << endl;
	 for (int j=0;j<NUMPOINTS;j++)
	 {
	    M.renewBoundary(j,S[j].Hgrid->hexen);
	 }
     cout << " just after renewBoundary  "  << endl;
	  // redissect the boundaries
	 for (int j=0;j<NUMPOINTS;j++)
	 {
	    M.renewDissect(j);
     }
     cout << " just after dissectBoundary "  << endl;

      gfile << endl << "analysis on first morphing iteration " << endl;
      for (int j=0;j<NUMPOINTS;j++) {
              if (M.regArea(j) != 0){
			      int regionCount = 0;
                  gfile<<"in the degree loop" << endl;
                  //angle degree
                  tempArea = M.regArea(j);
                  tempPerimeter = M.renewRegPerimeter(j);

                  // digital version
                  angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors);
                  degreeAngle = L.find_zeroDAngle(angleDVector);
                  gfile << "region "<< j << " degreeDAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;

                  // analogue version
                  angleVector = M.sectorize_reg_angle(j,numSectors,radiusOffset, numSectors);
                  angleVector = L.meanzero_vector(angleVector);
                  //degreeAngle = M.find_max(angleVector,3);
                  degreeAngle = L.find_zeroAngle(angleVector,3);
                  gfile << "region "<< j << " degreeAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;
                  //radial degree
                  degreeRadius = -100;
                  radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2);
                  //gfile <<"after sectorize_reg_radius"<<endl;
                  // radiusVector = M.meanzero_vector(radiusVector);

                  degreeRadius = L.find_zeroDRadius(radiusDVector);
                  gfile  << "region "<< j << " degreeDRadius "<< degreeRadius << "  " <<endl ;

                  ///radial degree
                  degreeRadius = -100;
                  int newdegreeRadius = 0;
                  for (int angleOffset=0; angleOffset<numSectors -1; angleOffset += 3)
				  {
			          radiusVector = M.sectorize_reg_radius(j,numSectors, angleOffset, angleOffset + numSectors/2);
			          newdegreeRadius = L.find_zeroRadius(radiusVector,3);
			          if (newdegreeRadius > degreeRadius)
				      degreeRadius = newdegreeRadius;
		          }


                  gfile <<  " region "<< j << " degreeRadius  "<< degreeRadius << "  " <<endl << endl;


                  regionCount++;

	      } //end of if on non-zero regions
} //end of loop on NUMPOINT



	      avDegreeAngle = 0;
	      avDegreeRadius = 0;
	      occupancy = 0;
	      countRegions = 0;
		  tempArea = 0;
		  tempPerimeter = 0;
		  avAbsCorrelation = 0;
		  cout << "just after renewcorrelate_edges morph1 " << endl;
		  M.random_correlate(max_comp,1);
		  cout << "just after randomcorrelate_edges morph1 " << endl;
          for (int j=0;j<NUMPOINTS;j++) {
	          countRegions++;
		      occupancy += M.regNNfrac(j);
              tempArea = M.regArea(j);
              tempPerimeter = M.regPerimeter(j);

              avAbsCorrelation += M.renewcorrelate_edges(j,1);
              angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors);
              degreeAngle = L.find_zeroDAngle(angleDVector);
		      avDegreeAngle += degreeAngle;
              //radial degree
		      degreeRadius = 0;
              radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2);
              degreeRadius = L.find_zeroDRadius(radiusDVector);
              avDegreeRadius += degreeRadius;

              degfile2 << degreeAngle/2 << " " << degreeRadius << " " << M.regNNfrac(j) << " " << tempArea << " "<< tempPerimeter<<endl<<flush;
	    } //end of loop on NUMPOINTs
        avDegreeAngle = avDegreeAngle / (1.0 * countRegions);
        avDegreeRadius = avDegreeRadius / (1.0 * countRegions);
        avAbsCorrelation = avAbsCorrelation / (1.0 * countRegions);
	    occupancy = occupancy / (1.0 * countRegions);
	    // jfile << avDegreeAngle <<" "<<avDegreeRadius<<endl;
	    jfile <<Dn<< " "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation<< " after morphing  " << endl;
//end of integration after first morphing
//       if (skipMorph) return 0;
//begin second morphing
// now set the boundaries for the regions, stored in curvedBoundary
    for (int j=0; j<NUMPOINTS;j++) {
        vector<hexGeometry::lineSegment> regionLineSegments;
        regionLineSegments = M.polygonSides(j, true);
        unsigned int sideSize = regionLineSegments.size();
        for (auto& h : S[j].Hgrid->hexen) {
            for (unsigned int i=0; i<sideSize; i++) {
                if (M.hGeo->hexIntersectLineSegment(regionLineSegments[i], h)) {
                        h.setFlags(HEX_IS_REGION_BOUNDARY);
                }
            }
        }
    }
    M.populateBoundPolygon(1);
    std::vector<std::vector<std::list<morph::Hex>::iterator>> regionpHexList; //vector of vectors of pHexes (iterators)
    regionpHexList.resize(NUMPOINTS);
    cout << "just after setting curved boundaries second iteration " << M.curvedBoundary.size()<<endl;
    std::pair<float,float> regionCentroid[NUMPOINTS];
    for (int j = 0;j<7;j++) {
        morph::BezCurvePath<float> bound = M.curvedBoundary[j];
        cout << "after assigning BexCurvePath bound j " << j << endl;
        //regionpHexList[j] = S[j].Hgrid->getRegion(M.curvedBoundary[j], regionCentroid[j], true);
        S[j].Hgrid->setBoundaryOnly(bound, false);
        cout << " after setting BoundaryOnly j " << j << endl;
    }
    for (int j = 8;j<NUMPOINTS;j++) {
        morph::BezCurvePath<float> bound = M.curvedBoundary[j];
        cout << "after assigning BexCurvePath bound j " << j << endl;
        //regionpHexList[j] = S[j].Hgrid->getRegion(M.curvedBoundary[j], regionCentroid[j], true);
        S[j].Hgrid->setBoundaryOnly(bound, false);
        cout << " after setting BoundaryOnly j " << j << endl;
    }
    cout << "just after setting the polygonal boundary"<<endl;
// now draw the intial tesselation
    morph::Gdisplay ndisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
    ndisp.resetDisplay (fix, eye, rot);
    for (int j=0;j<NUMPOINTS;j++)
    {
        // plot stuff here.
        int boundaryCount = 0;
        int internalCount = 0;
        int totalCount = 0;
        cout << "size of region j "<< j << " is " << M.regionpHexList[j].size() << " size of Hexgrid " << S[j].Hgrid->num() << endl;
        //cout << "size of region j "<< j << " is " << S[j].Hgrid->num() << " size of Hexgrid " << S[j].Hgrid->num() << endl;
        //for (auto h : regionpHexList[j]) {
        for (auto h : S[j].Hgrid->hexen) {
            totalCount++;
            if (h.testFlags(HEX_IS_REGION_BOUNDARY) == true)
        //   if (h.boundaryHex())
            {
                ndisp.drawHex (h.position(), (h.d/2.0f), cl_a);
                boundaryCount++;
            }
            else
            {
                //cout << "h.internal hex in region " << j <<  " h.vi " << h.vi << endl;
    		    //cout << "region " << j <<  " x " << S[j].Hgrid->d_x[h.vi] << " y " << S[j].Hgrid->d_y[h.vi] << endl;
               // ndisp.drawHex (h->position(), offset, (h->d/2.0f), cl_b);
                internalCount++;
            }
            }
    		cout << "boundaryCount 2 "<<boundaryCount<< " internalCount 2 "<< internalCount << " totalCount 2 = " << totalCount << "  region size " << regionpHexList[j].size() << endl;


            // Draw small hex at boundary centroid
            array<float,3> c;
            c[2] = 0;
            c[0] = S[j].Hgrid->boundaryCentroid.first;
            c[1] = S[j].Hgrid->boundaryCentroid.second;
            cout << "d/2: " << S[j].Hgrid->hexen.begin()->d/4.0f << endl;
            //ndisp.drawHex (c, offset, (S[j].Hgrid->hexen.begin()->d/2.0f), cl_a);
            cout << "boundaryCentroid x,y: " << c[0] << "," << c[1] << endl;
            cout << "Region centroid x "<< regionCentroid[j].first << " y " << regionCentroid[j].second <<endl;
    }//end of loop on NUMPOINTS to plot boundaries
    usleep (1000000);
    cout << "before redrawDisplay 2 " << endl;
    ndisp.redrawDisplay();
    cout << "after redrawDisplay 2" << endl;
    usleep (1000000); // one hundred seconds
    ndisp.saveImage(logpath + "/Tesselation3.png");
    ndisp.closeDisplay();

      if (skipMorph) {
          return 0;
      }
     //begin second time stepping loop
	/*
    for (int i=0;i<numsteps;i++)
	{
	  //cout << " head of second time stepping loop i " << i << endl;
   	  for (int j = 0;j<NUMPOINTS;j++) //loop over regions
	  {
     	   S[j].step(dt, Dn, Dchi, Dc);
     	 //  S[j].step(dt, Dn, Dchi, Dc);
	  }
     */
	      int countHex = 0;
          //set up display
          morph::Gdisplay nndisp(900, 900, 0, 0, "morph 2 u run", rhoInit, 0.0, 0.0);
          nndisp.resetDisplay (fix, eye, rot);
   	      for (int j = 0;j<NUMPOINTS;j++) //loop over regions
	      {
	 	    cout << "in print routine"<<endl;
            unsigned int regsize = S[j].Hgrid->num();
		    cout << "in loop over regions " << " size is " << regsize << endl;
			countHex += regsize;
		    vector<double> regionNN(regsize);
			vector<double> tempNN(regsize);
			int k = 0;
		    for (auto h : S[j].Hgrid->hexen)
			{
			  tempNN[k] = S[j].NN[h.vi];
			  k++;
			}
//normalise over the region then write normalised values to normalised NN over hexGrid
            regionNN = L.normalise(tempNN);
		    cout << "total number of hexes counted " << countHex << endl;
		      //cout << "just before drawHex 2 "<< h.vi << "normalNN " << normalNN[h.vi] << endl;
		    for (auto h : S[j].Hgrid->hexen)
			{
		      //if (!h.boundaryHex())
			  //{
			    array<float,3> colour = morph::Tools::getJetColorF(regionNN[h.vi]);
	            nndisp.drawHex(h.position(),(h.d/2.0f),colour);
		        //cout << "just after drawHex 2"<<endl;
			   //}
            }
          }//end of loop over regions
		   cout << "just berore redraw display 2  "  << endl;
           nndisp.redrawDisplay();
		   cout << "just after redraw display 2  "  << endl;
           usleep (1000000); // one hundred seconds
		   nndisp.saveImage(logpath + "/nnField3.png");
		   cout << "just after saveImage "  << endl;
           nndisp.closeDisplay();
		   cout << "just after closeDisplay "  << endl;



	 // repopulate the regions
	 for (int j=0;j<NUMPOINTS;j++)
	 {
	    M.renewRegion(j,S[j].Hgrid->hexen);
     }
	 for (int j=0;j<NUMPOINTS;j++)
	 {
	    M.renewBoundary(j,S[j].Hgrid->hexen);
     }
	 // redissect the boundaries
	 for (int j=0;j<NUMPOINTS;j++)
	 {
	    M.renewDissect(j);
     }
	  /*
	  for (int j=0;j<NUMPOINTS;j++)
	  {
	    double correlation = M.renewcorrelate_edges(j,logpath);
		cout << " correlation region second morph " << j << " = " << correlation << endl;
	  }
	  */

      gfile << endl << "analysis on second morphing iteration " << endl;
      for (int j=0;j<NUMPOINTS;j++)
	  {
	         int regionCount = 0;
             gfile<<"in the degree loop" << endl;
	         //angle degree
             tempArea = M.regArea(j);
             tempPerimeter = M.renewRegPerimeter(j);
             // digital version
             angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors);
             degreeAngle = L.find_zeroDAngle(angleDVector);
             gfile << "region "<< j << " degreeDAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;

           // analogue version
              angleVector = M.sectorize_reg_angle(j,numSectors,radiusOffset, numSectors);
              angleVector = L.meanzero_vector(angleVector);
              //degreeAngle = M.find_max(angleVector,3);
              degreeAngle = L.find_zeroAngle(angleVector,3);
              gfile << "region "<< j << " degreeAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;
              //radial degree
              degreeRadius = -100;
              radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2);
              //gfile <<"after sectorize_reg_radius"<<endl;
              // radiusVector = M.meanzero_vector(radiusVector);
              degreeRadius = L.find_zeroDRadius(radiusDVector);
              gfile  << "region "<< j << " degreeDRadius "<< degreeRadius << "  " <<endl ;

                  ///radial degree
              degreeRadius = -100;
              int newdegreeRadius = 0;
			  for (int angleOffset=0; angleOffset<numSectors -1; angleOffset += 3)
		      {
			      radiusVector = M.sectorize_reg_radius(j,numSectors, angleOffset, angleOffset + numSectors/2);
			      newdegreeRadius = L.find_zeroRadius(radiusVector,3);
			      if (newdegreeRadius > degreeRadius)
				          degreeRadius = newdegreeRadius;
		      }
                  gfile <<  " region "<< j << " degreeRadius  "<< degreeRadius << "  " <<endl << endl;
                  regionCount++;
     } //end of loop on NUMPOINT

//computing average values over all regions
	  avDegreeAngle = 0;
	  avDegreeRadius = 0;
	  occupancy = 0;
      avAbsCorrelation = 0;
	  tempPerimeter = 0;
	  tempArea = 0;
	  countRegions = 0;
	  cout << "just after renewcorrelate_edges morph2 " << endl;
	  M.random_correlate(max_comp, 2);
	  cout << "just after random correlate_edges morph2 " << endl;
      for (int j=0;j<NUMPOINTS;j++)
	  {
	      countRegions++;
          avAbsCorrelation += M.renewcorrelate_edges(j,2);
		  occupancy += M.regNNfrac(j);
		  cout << "after renNNfrac " << endl;
          tempArea = M.regArea(j);
          cout << " after regArea " << endl;
          tempPerimeter = M.regPerimeter(j);
          cout << " after regPerimeter " << endl;
          cout << "before sectorixe Dangle morph2 " << endl;
          angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors);
          cout << "after sectorixe Dangle morph2 " << endl;
          degreeAngle = L.find_zeroDAngle(angleDVector);
		  avDegreeAngle += degreeAngle;
          cout << "after sectorixe Dangle morph2 " << endl;
		  degreeRadius = 0;
          radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2);
          degreeRadius = L.find_zeroDRadius(radiusDVector);
          avDegreeRadius += degreeRadius;
          degfile3 << degreeAngle << " " << degreeRadius << " " << M.regNNfrac(j) << " " << tempArea << " " << tempPerimeter<<endl<<flush;
	    } //end of loop on NUMPOINTs
        avDegreeAngle = avDegreeAngle / (1.0 * countRegions);
        avDegreeRadius = avDegreeRadius / (1.0 * countRegions);
	    occupancy = occupancy / (1.0 * countRegions);
		avAbsCorrelation = avAbsCorrelation/(1.0 * countRegions);
	    // jfile << avDegreeAngle <<" "<<avDegreeRadius<<endl;
	    jfile <<Dn<<"  "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation<< " after morphing  2" <<endl;
    return 0;
}
