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
    if (argc < 5) {
      std::cout << "not enough arguments" << argc << endl;
      return -1;
    }
    string jsonfile = argv[1];
    string iter = argv[2];
    bool Lcontinue = stoi(argv[3]);
    int numsteps = stoi(argv[4]);
    int numprint = stoi(argv[5]);
    string logpath = argv[6];
    // adjust the number of steps according to the Dn number
    numsteps = numsteps * ceil(Dn/36.0);
    numprint  = numprint * ceil(Dn/36.0);
    //  open the confgig file and read in the parameters
    morph::Config conf(jsonfile);
    if (!conf.ready) {
        cerr << "Error setting up JSON config: " << conf.emsg << endl;
    }
#ifdef SINGLE
        float dt = conf.getFloat("dt",0.0001);
        float Dn = conf.getFloat("Dn",1.0);
        float Dchi = conf.getFloat("Dchi",0.0);
        float Dc = conf.getFloat("Dc",0.3);
        float xspan = conf.getFloat("xspan",5.0);
        float boundaryFalloffDist = conf.getFloat("boundaryFalloffDist",0.0078);
        float aNoiseGain = conf.getFloat("aNoiseGain",0.1);
        float nnInitialOffset = conf.getFloat("nnInitialOffet", 1.0);
        float ccInitialOffset = conf.getFloat("ccInitialOffset",2.5);
#else
        double dt = conf.getDouble("dt",0.0001);
        double Dn = conf.getDouble("Dn",1.0);
        double Dchi = conf.getDouble("Dchi",0.0);
        double Dc = conf.getDouble("Dc",0.3);
        double xspan = conf.getDouble("xspan",5.0);
        double boundaryFalloffDist = conf.getDouble("boundaryFalloffDist",0.0078);
        double aNoiseGain = conf.getDouble("aNoiseGain",0.1);
        double nnInitialOffset = conf.getDouble("nnInitialOffet", 1.0);
        double ccInitialOffset = conf.getDouble("ccInitialOffset",2.5);
#endif
    int numSectors = conf.getInt("numsectors",12);
    int scale = conf.getInt("scale",8);
    bool Lgraphics = conf.getBool("Lgraphics",false);
    bool LfixedSeed = conf.getBool("LfixedSeed",0);
    bool LDn = conf.getBool("LDn",0);
    bool overwrite_logs = conf.getBool("overwrite_logs",true);
    bool skipMorph  = conf.getBool("skipMorph",false);
    unsigned int NUMPOINTS = conf.getInt("NUMPOINTS",41);
    cout << " Lcontinue " << Lcontinue << " skipMorph " << skipMorph << endl;
    ofstream afile (logpath + "/centroids.out",ios::app);

    unsigned int seed;
    if (LfixedSeed) {
        seed = 1;
    }
    else {
        seed = time(NULL);
    }

    // A ra2yyndo2yym uniform generator returning real/FLTing point types
    morph::RandUniform<FLT> ruf(seed);

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
    DRegion M(scale,xspan,logpath,NUMPOINTS); //create tessellation
    M.setCreg(); //set counts to identify inner boundaries
    M.setInternalBoundary(); //set internal boundaries
    cout << "before dissect_boundary " << endl;
//	for (auto h : M.Hgrid->hexen) {
//	   cout << " first h.vi output " << h.vi << endl;
//	   }
    vector<std::pair<FLT,FLT>> cGravity;
    cGravity = M.dissectBoundary(); //dissect region boundary
    cout << "Edges size = " << M.edges.size() << endl;
	M.setRadialSegments(); //set the radial segments for regions
    int inReg = 0;
    inReg = M.setInnerRegion(); //set mask array for inner regions
    int rC = 0;
    //check for correct number of inner regions
    for (unsigned int j=0; j<NUMPOINTS;j++) {
        if (!M.innerRegion[j]) rC++;
    }
    if (inReg != rC) {
        cout << "Error: setInnerRegion returns " << inReg << " but outerRegions " << rC << endl;
        return -1;
    }
    else {
        cout << "Success: setInnerRegion no of outer regions " << inReg << endl;
    }
    cout << "after first setRadialSegments " << endl;
// include the analysis methods
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
    M.populateBoundPolygon(1);
    cout << "just after setting polygonal boundaries " << M.curvedBoundary.size()<<endl;
    for (unsigned int j = 0;j<NUMPOINTS;j++) {
        S.push_back(ksSolver(scale, xspan, logpath, M.curvedBoundary[j], M.centroids[j]));
        cout << "in the loop populating the ksVector morph0 "<< j <<endl;
    }
    cout << "first setting of centroids" << endl;
    for (unsigned int j=0; j<NUMPOINTS;j++){
        afile << "centroid region " << j << " is ( " << M.centroids[j].first << " , " << M.centroids[j].second << " )" << endl;
    }
    // repopulate the regions
    cout << "just before populating the inner regions Morph 0" << endl;
    for (unsigned int j=0;j<NUMPOINTS;j++)
    {
        if (M.innerRegion[j]) {
            M.renewRegion(j,S[j].Hgrid->hexen);
        }
    }
    cout << "just before populating the inner boundary Morph 0" << endl;
    for (unsigned int j=0;j<NUMPOINTS;j++)
    {
        if (M.innerRegion[j]) {
            M.renewBoundary(j,S[j].Hgrid->hexen);
        }
        //M.renewCentroids(j);
    }
    cout << "just before calculating regionSize " << endl;
    for (unsigned int j=0; j<NUMPOINTS;j++){
        afile << " Region size " << M.regionHex[j].size() << endl;
    }
	cout << "just after populating the regions from the ksSolver vector"<<endl;
    //clear global edges map
    M.edges_clear();
    // swap the radialAngles to the mCoords
   // M.swapRadialSegments(false);
    // redissect the boundaries
    cout << "just before renewDissect first time" << endl;
    for (unsigned int j=0;j<NUMPOINTS;j++)
    {
        if (M.innerRegion[j]) {
            M.renewDissect(j,0);
        }
    }
    cout << "Edges size " << M.edges.size() << endl;
// now draw the intial tesselation
    int internalCount = 0;
    int boundaryCount = 0;
    FLT hexWidth = M.Hgrid->hexen.begin()->d/2.0;
    cerr << "d/2: " << hexWidth << endl;
#ifdef COMPILE_PLOTTING
    float rhoInit = 2.5;
    vector<double> fix(3.0, 0.0);
    vector<double> eye(3.0, 0.0);
    vector<double> rot(3.0, 0.0);
    array<float,3> cl_a = morph::Gdisplay::getJetColorF (0.78); //orange
    array<float,3> cl_c = morph::Gdisplay::getJetColorF (0.28); //blue
    array<float,3> cl_b = morph::Gdisplay::getJetColorF (0.58); //yellow
    array<float,3> cl_d = morph::Gdisplay::getJetColorF (0.00); //black
    array<float,3> offset = {{0, 0, 0}};
    if (Lgraphics) {
        morph::Gdisplay idisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
        idisp.resetDisplay (fix, eye, rot);
        for (unsigned int j=0;j<NUMPOINTS;j++) {
        // plot stuff here.
            cout << "size of Hexgrid j "<< j << " is " << S[j].Hgrid->num() << endl;
            for (auto h : S[j].Hgrid->hexen) {
                if (h.boundaryHex()) {
                    cout << "h.boundaryHex in region " << j <<  " h.vi " << h.vi << endl;
// cout << "region " << j <<  " x " << S[j].Hgrid->d_x[h.vi] << " y " << S[j].Hgrid->d_y[h.vi] << endl;
                    idisp.drawHex (h.position(), (h.d/2.0f), cl_d);
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
                std::array<FLT,3> pos = {M.vCoords[j][i].first, M.vCoords[j][i].second, 0};
                idisp.drawHex (pos,offset,4.0*hexWidth,cl_d);
                pos = {M.mCoords[j][i].first, M.mCoords[j][i].second, 0};
                idisp.drawHex (pos,offset,4.0*hexWidth,cl_d);
            }

        } // end of loop over regions
        for (auto h : M.Hgrid->hexen) {
            if (M.Creg[h.vi] ==  1 ) {
                idisp.drawHex (h.position(), (h.d/2.0f), cl_d);
                boundaryCount++;
            }
            else if (M.Creg[h.vi] > 1) {
                idisp.drawHex (h.position(), (h.d/2.0f), cl_d);
            }
        }
       //  else
        // {
          //   isp.drawHex (h.position(), offset, (h.d/2.0f), cl_b);
         //}
      //end of for loop on j as region number
        usleep (1000000);
        cout << "before redrawDisplay 2 " << endl;
        idisp.redrawDisplay();
        cout << "after redrawDisplay 2" << endl;
        usleep (100000); // one hundred seconds
        idisp.saveImage(logpath + "/Tesselation1" + iter + ".png");
        idisp.closeDisplay();
    } //end of graphics section
#endif
// initialise the fields
    string fname = logpath + "/first.h5";
    cout<< "just before first data read morph 0 "<< " lcontinue " << Lcontinue <<endl;
// initialise with random field
    if (Lcontinue)
	{
        morph::HdfData ginput(fname,1);
        cout << "just after trying to open ../logs/first.h5" << endl;
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
		        FLT choice = ruf.get();
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
                FLT bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary- boundaryFalloffDist)) );
                S[j].NN[h.vi] = (S[j].NN[h.vi] - nnInitialOffset) * bSig + nnInitialOffset;
                S[j].CC[h.vi] = (S[j].CC[h.vi] - ccInitialOffset) * bSig + ccInitialOffset;
		 } //end of if on boundary distance
	    }//end of loop over region
	   }//end of loop over all regions
      } //end of else on Lcontinue
     cout <<  "just after field creation first morph" << endl;

    // begin second time stepping loop unmorphed grid solved via ksSolver
    for (int i=0;i<numsteps;i++)
    {
	    int countHex = 0;
   	    for (unsigned int j = 0;j<NUMPOINTS;j++) { //loop over all regions, only step internal ones
            if (M.innerRegion[j]) {
                S[j].step(dt, Dn, Dchi, Dc);
            }
        }
   #ifdef COMPILE_PLOTTING
	    if (i%numprint == 0 && Lgraphics) //only display internal regions
	    {
          //set up display
            morph::Gdisplay iidisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
            iidisp.resetDisplay (fix, eye, rot);
   	        for (unsigned int j = 0;j<NUMPOINTS;j++) //loop over regions
	        {
                if (M.innerRegion[j]) {
	 	        cout << "in print routine"<<endl;
                unsigned int regsize = S[j].Hgrid->num();
		        cout << "in loop over regions " << " size is " << regsize << endl;
			    countHex += regsize;
		        vector<FLT> regionNN;
			    vector<FLT> tempNN;
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
                    array<FLT,3> colour = morph::Gdisplay::getJetColorF(regionNN[idx]);
                    if (!h.boundaryHex()) {
                        iidisp.drawHex(h.position(),(h.d/2.0f),colour);
                    //    cout << "just after drawHex 2 colour "<< " value " << regionNN[idx] << endl;
                    }
                    idx++;
                }
                for (unsigned int i=0; i< M.vCoords[j].size();i++) {
                    std::array<FLT,3> pos = {M.vCoords[j][i].first, M.vCoords[j][i].second, 0};
                    iidisp.drawHex (pos,offset,4.0*hexWidth,cl_d);
                }
                }//end of loop on inner regions
            }//end of loop over regions
    for (auto h : M.Hgrid->hexen) {
         if (M.Creg[h.vi] ==  1 ) {
             iidisp.drawHex (h.position(), (h.d/2.0f), cl_d);
             boundaryCount++;
         }
         else if (M.Creg[h.vi] > 1) {
              iidisp.drawHex (h.position(), (h.d/2.0f), cl_d);
         }
    }
            cout << "just before redraw display 1" << endl;
            iidisp.redrawDisplay();
		    cout << "just after redraw display 1" << endl;
		    cout << "just after to_string"<<endl;
		    iidisp.saveImage(logpath + "/nnField1" + iter + ".png");
            usleep (1000000); // one hundred seconds
		    cout << "just after saveImage 1" << endl;
		    iidisp.closeDisplay();
		    cout << "just after close display 1 i " << i << endl;
        }//end of loop on numprint drawing fields
   #endif
		 //cout << "iteration " << i << " of morph 1 time step" << endl;
    } // end of second time stepping

    //code run at end of timestepping
    //first save the  ofstream outFile;
    cout << "just before first data write morph 0" << endl;
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
    cout << " just after writing data "  << endl;
    gfile << endl << "analysis on first morphing iteration " << endl;

    //declaration of variables needed for analysis
    vector <int> radiusDVector;
    vector <int> angleDVector;
	vector <FLT> angleVector;
    vector <FLT> radiusVector;
    int degreeRadius;
    int degreeAngle;
	FLT tempArea = 0.0;
	FLT tempPerimeter = 0.0;
    int angleOffset = 0;
    int radiusOffset = 0;
	FLT avDegreeRadius = 0.0;
	FLT avDegreeAngle = 0;
	FLT occupancy = 0.0;
    int countRegions = 0;
    FLT avAbsCorrelation = 0.0;
    const int max_comp = NUMPOINTS*3;
    for (unsigned int j=0;j<NUMPOINTS;j++) {
        M.NN[j] = S[j].NN; //write the NN fields to the DRegion array
        if (M.innerRegion[j]){
			int regionCount = 0;
            gfile<<"in the degree loop" << endl;
            //angle degree
            tempArea = M.regArea(j);
            tempPerimeter = M.renewRegPerimeter(j);

            // digital version
            angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, S[j].NN);
            degreeAngle = L.find_zeroDAngle(angleDVector);
            gfile << "region "<< j << " degreeDAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;

            // analogue version
            angleVector = M.sectorize_reg_angle(j,numSectors,radiusOffset, numSectors, S[j].NN);
            angleVector = L.meanzero_vector(angleVector);
            //degreeAngle = M.find_max(angleVector,3);
            degreeAngle = L.find_zeroAngle(angleVector,3);
            gfile << "region "<< j << " degreeAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;

            //radial degree
            degreeRadius = -100;
            radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
            //gfile <<"after sectorize_reg_radius"<<endl;
            // radiusVector = M.meanzero_vector(radiusVector);
            degreeRadius = L.find_zeroDRadius(radiusDVector);
            gfile  << "region "<< j << " degreeDRadius "<< degreeRadius << "  " <<endl ;

            //radial degree
            degreeRadius = -100;
            int newdegreeRadius = 0;
            for (int angleOffset=0; angleOffset<numSectors -1; angleOffset += 3) {
			    radiusVector = M.sectorize_reg_radius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
			    newdegreeRadius = L.find_zeroRadius(radiusVector,3);
			    if (newdegreeRadius > degreeRadius) {
                    degreeRadius = newdegreeRadius;
                }
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
          //avAbsCorrelation = M.correlate_edges(0);
		  M.random_correlate(max_comp, 1);
		  cout << "just after randomcorrelate_edges morph1 " << endl;
          for (unsigned int j=0;j<NUMPOINTS;j++) {
	        if (M.innerRegion[j]) {
	          countRegions++;
		      occupancy += M.regNNfrac(j);
              tempArea = M.regArea(j);
              tempPerimeter = M.regPerimeter(j);

              avAbsCorrelation += M.renewcorrelate_edges(j,1);
              angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, S[j].NN);
              degreeAngle = L.find_zeroDAngle(angleDVector);
		      avDegreeAngle += degreeAngle;
              //radial degree
		      degreeRadius = 0;
              radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
              degreeRadius = L.find_zeroDRadius(radiusDVector);
              avDegreeRadius += degreeRadius;

              degfile1 << degreeAngle/2 << " " << degreeRadius << " " << M.regNNfrac(j) << " " << tempArea << " "<< tempPerimeter<<endl<<flush;
	       } //end of if on non-zero regions
	    } //end of loop on NUMPOINTs
      if (countRegions == 0) {
          cout << "Error zero regionss counted in second analysis morph 0" << endl;
          return -1;
      }
        avDegreeAngle = avDegreeAngle / (1.0 * countRegions);
        avDegreeRadius = avDegreeRadius / (1.0 * countRegions);
        avAbsCorrelation = avAbsCorrelation / (1.0 * countRegions);
	    occupancy = occupancy / (1.0 * countRegions);
        cout << "Regions counted in first analysis loop " << countRegions << endl;
	    // jfile << avDegreeAngle <<" "<<avDegreeRadius<<endl;
	    jfile <<Dn<< " "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation << endl;

        // write the edge vectors all interpolated to a uniform size
        std::map<int, vector<FLT>>::iterator ptr;
        std::string sideNN = logpath + "/edgeNN.data";
        int ecount = 0;
        for (ptr = M.edgeNN.begin(); ptr != M.edgeNN.end(); ptr++) {
            vector<FLT> tempVect;
            tempVect = M.normaliseVect(ptr->second);
            M.printFLTVect(sideNN, tempVect);
            cout << "edgeNN key " << ptr->first << " element " << ecount << " vector size " << ptr->second.size() << " edges size " << M.edges[ptr->first].size() <<  endl;
            ecount++;
        }
//end of integration after solving on polygonal regions via ksSolver

    /*
     * now create an array of morphed regions, first morph
     */

    if (skipMorph) return 0;
    cout << "just before setting curved boundaries first morph" <<endl;
// section for solving on the curved boundaries
    S.resize(0);
    //set radius for creating circular regions also calculate the adjusted Dn values
    vector<FLT> DnVal;
    DnVal.resize(NUMPOINTS,0.0);
    vector<FLT> DchiVal;
    DchiVal.resize(NUMPOINTS,0.0);
    vector<FLT> DcVal;
    DcVal.resize(NUMPOINTS,0.0);
    vector<FLT> morph0Area;
    morph0Area.resize(NUMPOINTS,0.0);
	for (unsigned int j = 0;j<NUMPOINTS;j++) {
        morph0Area[j] = M.hexArea*M.regArea(j);
    }
// now set the boundaries for the regions, stored in curvedBoundary
    M.populateBoundCurve(1);
    cout << "just after setting curved boundaries morph 1 " << M.curvedBoundary.size()<<endl;
    for (unsigned int j = 0;j<NUMPOINTS;j++) {
        S.push_back(ksSolver(scale, xspan, logpath, M.curvedBoundary[j], M.centroids[j]));
        cout << "in the loop populating the ksVector morph1   "<< j <<endl;
    }
    cout << "first morph  setting of centroids" << endl;
    for (unsigned int j=0; j<NUMPOINTS;j++){
        afile << "centroid region " << j << " is ( " << M.centroids[j].first << " , " << M.centroids[j].second << " )" << endl;
    }
    // repopulate the regions
    cout << "just before populating the inner regions morph 1" << endl;
    for (unsigned int j=0;j<NUMPOINTS;j++)
    {
        if (M.innerRegion[j]) {
            M.renewRegion(j,S[j].Hgrid->hexen);
        }
    }
    cout << "just before populating the inner boundary morph 1" << endl;
    for (unsigned int j=0;j<NUMPOINTS;j++)
    {
        if (M.innerRegion[j]) {
            M.renewBoundary(j,S[j].Hgrid->hexen);
        }
    }
    cout << "second setting of centroids" << endl;
    for (unsigned int j=0; j<NUMPOINTS;j++){
        afile << "centroid region " << j << " is ( " << M.centroids[j].first << " , " << M.centroids[j].second << " )" << endl;
    }
	cout << "just after populating the ksVector"<<endl;
    for (unsigned int j = 0;j<NUMPOINTS;j++) {
	//	for (auto h : S[j].Hgrid->hexen) {
	//	    cout << "hexNumber in region " << j <<  " h.vi " << h.vi << endl;
    //    }
        FLT area = M.hexArea*M.regArea(j);
        if (LDn) {
            DchiVal[j] = Dchi * morph0Area[j] / area;
            DnVal[j] = Dn *  morph0Area[j] / area;
            DcVal[j] = Dc * morph0Area[j] / area;
        }
        else {
            DchiVal[j] = Dchi;
            DnVal[j] = Dn;
            DcVal[j] = Dc;
        }

        cout << "DnVal region " << j << " = " << DnVal[j] << " Dn = " << Dn << " PI " << PI << endl;
    }
// now draw the intial tesselation
    internalCount = 0;
    boundaryCount = 0;
#ifdef COMPILE_PLOTTING
    if (Lgraphics) {
        morph::Gdisplay mdisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
        mdisp.resetDisplay (fix, eye, rot);
        for (unsigned int j=0;j<NUMPOINTS;j++) {
        // plot stuff here.
            cout << "size of Hexgrid j "<< j << " is " << S[j].Hgrid->num() << endl;
            for (auto h : S[j].Hgrid->hexen) {
                if (h.boundaryHex()) {
                    cout << "h.boundaryHex in region " << j <<  " h.vi " << h.vi << endl;
// cout << "region " << j <<  " x " << S[j].Hgrid->d_x[h.vi] << " y " << S[j].Hgrid->d_y[h.vi] << endl;
                    mdisp.drawHex (h.position(), (h.d/2.0f), cl_d);
                    boundaryCount++;
                }
                else
                {
                    cout << "h.internal hex in region " << j <<  " h.vi " << h.vi << endl;
//cout << "region " << j <<  " x " << S[j].Hgrid->d_x[h.vi] << " y " << S[j].Hgrid->d_y[h.vi] << endl;
                //mdisp.drawHex (h.position(), offset, (h.d/2.0f), cl_b);
                    internalCount++;
                }
            }
            for (unsigned int i=0; i< M.vCoords[j].size();i++) {
                std::array<FLT,3> pos = {M.vCoords[j][i].first, M.vCoords[j][i].second, 0};
                mdisp.drawHex (pos,offset,4.0*hexWidth,cl_d);
                pos = {M.mCoords[j][i].first, M.mCoords[j][i].second, 0};
                mdisp.drawHex (pos,offset,4.0*hexWidth,cl_d);
            }

        } // end of loop over regions
        for (auto h : M.Hgrid->hexen) {
            if (M.Creg[h.vi] ==  1 ) {
                mdisp.drawHex (h.position(), (h.d/2.0f), cl_d);
                boundaryCount++;
            }
            else if (M.Creg[h.vi] > 1) {
                mdisp.drawHex (h.position(), (h.d/2.0f), cl_d);
            }
        }
       //  else
        // {
          //   isp.drawHex (h.position(), offset, (h.d/2.0f), cl_b);
         //}
      //end of for loop on j as region number
        usleep (1000000);
        cout << "before redrawDisplay 2 " << endl;
        mdisp.redrawDisplay();
        cout << "after redrawDisplay 2" << endl;
        usleep (100000); // one hundred seconds
        mdisp.saveImage(logpath + "/Tesselation2" + iter + ".png");
        mdisp.closeDisplay();
    } //end of graphics section
#endif
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
		        FLT choice = ruf.get();
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
                FLT bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary- boundaryFalloffDist)) );
                S[j].NN[h.vi] = (S[j].NN[h.vi] - nnInitialOffset) * bSig + nnInitialOffset;
                S[j].CC[h.vi] = (S[j].CC[h.vi] - ccInitialOffset) * bSig + ccInitialOffset;
		 } //end of if on boundary distance
	    }//end of loop over region
	   }//end of loop over all regions
      } //end of else on Lcontinue
     cout <<  "just after field creation first morph" << endl;

    // begin third time stepping loop after first morph
    for (int i=0;i<numsteps;i++)
    {
	    int countHex = 0;
   	    for (unsigned int j = 0;j<NUMPOINTS;j++){ //loop over regions, only set internal ones
            if (M.innerRegion[j]) {
                S[j].step(dt, DnVal[j], DchiVal[j], DcVal[j]);
            }
        }
 #ifdef COMPILE_PLOTTING
	    if (i%numprint == 0 && Lgraphics) { //only display internal regions
          //set up display
            morph::Gdisplay mmdisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
            mmdisp.resetDisplay (fix, eye, rot);
   	        for (unsigned int j = 0;j<NUMPOINTS;j++) //loop over regions
	        {
                if (M.innerRegion[j]) {
	 	        cout << "in print routine"<<endl;
                unsigned int regsize = S[j].Hgrid->num();
		        cout << "in loop over regions " << " size is " << regsize << endl;
			    countHex += regsize;
		        vector<FLT> regionNN;
			    vector<FLT> tempNN;
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
                    array<FLT,3> colour = morph::Gdisplay::getJetColorF(regionNN[idx]);
                    if (!h.boundaryHex()) {
                        mmdisp.drawHex(h.position(),(h.d/2.0f),colour);
                    //    cout << "just after drawHex 2 colour "<< " value " << regionNN[idx] << endl;
                    }
                    idx++;
                }
                for (unsigned int i=0; i< M.vCoords[j].size();i++) {
                    std::array<FLT,3> pos = {M.vCoords[j][i].first, M.vCoords[j][i].second, 0};
                    mmdisp.drawHex (pos,offset,4.0*hexWidth,cl_d);
                }
                } //end of if on inner regions
            }//end of loop over regions
    for (auto h : M.Hgrid->hexen) {
         if (M.Creg[h.vi] ==  1 ) {
             mmdisp.drawHex (h.position(), (h.d/2.0f), cl_d);
             boundaryCount++;
         }
         else if (M.Creg[h.vi] > 1) {
              mmdisp.drawHex (h.position(), (h.d/2.0f), cl_d);
         }
    }
            cout << "just before redraw display 1" << endl;
            mmdisp.redrawDisplay();
		    cout << "just after redraw display 1" << endl;
		    cout << "just after to_string"<<endl;
		    mmdisp.saveImage(logpath + "/nnField2" + iter + ".png");
            usleep (1000000); // one hundred seconds
		    cout << "just after saveImage 1" << endl;
		    mmdisp.closeDisplay();
		    cout << "just after close display 1 i " << i << endl;
        }//end of loop on numprint drawing fields
    #endif
		 //cout << "iteration " << i << " of morph 1 time step" << endl;
    } // end of third time stepping

    //code run at end of timestepping
    //first save the  ofstream outFile;
    cout << "just before second data read " << endl;
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
     // write the NN and CC vals for each region
      gfile << endl << "analysis on first morphing iteration " << endl;
    radiusDVector.resize(0);
    angleDVector.resize(0);
    angleVector.resize(0);
    radiusVector.resize(0);
      for (unsigned int j=0;j<NUMPOINTS;j++) {
          M.NN[j] = S[j].NN; // first fill the DRegion NN values
              if (M.innerRegion[j]){
			      int regionCount = 0;
                  gfile<<"in the degree loop" << endl;
                  //angle degree
                  tempArea = M.regArea(j);
                  tempPerimeter = M.renewRegPerimeter(j);

                  // digital version
                  angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, S[j].NN);
                  degreeAngle = L.find_zeroDAngle(angleDVector);
                  gfile << "region "<< j << " degreeDAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;

                  // analogue version
                  angleVector = M.sectorize_reg_angle(j,numSectors,radiusOffset, numSectors, S[j].NN);
                  angleVector = L.meanzero_vector(angleVector);
                  //degreeAngle = M.find_max(angleVector,3);
                  degreeAngle = L.find_zeroAngle(angleVector,3);
                  gfile << "region "<< j << " degreeAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;
                  //radial degree
                  degreeRadius = -100;
                  radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
                  //gfile <<"after sectorize_reg_radius"<<endl;
                  // radiusVector = M.meanzero_vector(radiusVector);

                  degreeRadius = L.find_zeroDRadius(radiusDVector);
                  gfile  << "region "<< j << " degreeDRadius "<< degreeRadius << "  " <<endl ;

                  ///radial degree
                  degreeRadius = -100;
                  int newdegreeRadius = 0;
                  for (int angleOffset=0; angleOffset<numSectors -1; angleOffset += 3)
				  {
			          radiusVector = M.sectorize_reg_radius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
			          newdegreeRadius = L.find_zeroRadius(radiusVector,3);
			          if (newdegreeRadius > degreeRadius)
				      degreeRadius = newdegreeRadius;
		          }


                  gfile <<  " region "<< j << " degreeRadius  "<< degreeRadius << "  " <<endl << endl;


                  regionCount++;
        } //end of if on non-zero regions
    } //end of loop on NUMPOINT


    M.edges_clear();
    // swap the radialAngles to the mCoords
    // M.swapRadialSegments(false);
    // redissect the boundaries
    cout << "just before renewDissect second time" << endl;
    for (unsigned int j=0;j<NUMPOINTS;j++)
    {
        if (M.innerRegion[j]) {
            M.renewDissect(j,1);
        }
    }
    cout << "Edges size " << M.edges.size() << endl;


	      avDegreeAngle = 0;
	      avDegreeRadius = 0;
	      occupancy = 0;
	      countRegions = 0;
		  tempArea = 0;
		  tempPerimeter = 0;
		  avAbsCorrelation = 0;
		  M.random_correlate(max_comp,2);
		  cout << "just after randomcorrelate_edges morph1 " << endl;
          for (unsigned int j=0;j<NUMPOINTS;j++) {
	        if (M.innerRegion[j])
			{
	          countRegions++;
		      occupancy += M.regNNfrac(j);
              tempArea = M.regArea(j);
              tempPerimeter = M.regPerimeter(j);

              avAbsCorrelation += M.renewcorrelate_edges(j,2);
              angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, S[j].NN);
              degreeAngle = L.find_zeroDAngle(angleDVector);
		      avDegreeAngle += degreeAngle;
              //radial degree
		      degreeRadius = 0;
              radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
              degreeRadius = L.find_zeroDRadius(radiusDVector);
              avDegreeRadius += degreeRadius;

              degfile2 << degreeAngle/2 << " " << degreeRadius << " " << M.regNNfrac(j) << " " << tempArea << " "<< tempPerimeter<<endl<<flush;
	       } //end of if on non-zero regions
	    } //end of loop on NUMPOINTs
      if (countRegions == 0) {
          cout << "Error zero regionss counted in third analysis morph 1" << endl;
          return -1;
      }
        avDegreeAngle = avDegreeAngle / (1.0 * countRegions);
        avDegreeRadius = avDegreeRadius / (1.0 * countRegions);
        avAbsCorrelation = avAbsCorrelation / (1.0 * countRegions);
	    occupancy = occupancy / (1.0 * countRegions);
	    // jfile << avDegreeAngle <<" "<<avDegreeRadius<<endl;
	    jfile <<Dn<< " "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation << endl;
//end of integration after first morphing
//begin second morphing
    if (skipMorph) return 0;
// now set the boundaries for the regions, stored in curvedBoundary
    M.populateBoundCurve(0);
    S.resize(0);
    cout << "just after setting curved boundaries morph 2 " << M.curvedBoundary.size()<<endl;
    for (int j = 0;j<NUMPOINTS;j++)
    {
        std::pair<FLT, FLT> centroid =  M.baryCentre(j);
        S.push_back(ksSolver(scale, xspan, logpath, M.curvedBoundary[j], M.centroids[j]));
        cout << "in the loop populating the ksVector morph2 "<< j <<endl;
    }
    cout << "just after populating the ksVector"<<endl;
	// repopulate the regions
    for (unsigned int j=0;j<NUMPOINTS;j++)
    {
        if (M.innerRegion[j]) {
            M.renewRegion(j,S[j].Hgrid->hexen);
        }
    }
    cout << "after repopulate regions morph2" << endl;
    for (unsigned int j=0;j<NUMPOINTS;j++)
    {
        if (M.innerRegion[j]) {
            M.renewBoundary(j,S[j].Hgrid->hexen);
        }
        //M.renewCentroids(j);
    }
    cout << "after repopulate boundary morph2" << endl;
//    swap the radialAngles to the mCoords
    M.edges_clear();
    M.swapRadialSegments(false);
    // redissect the boundaries
    for (unsigned int j=0;j<NUMPOINTS;j++)
    {
        if (M.innerRegion[j]) {
            M.renewDissect(j,2);
        }
    }
    cout << "Edges size " << M.edges.size() << endl;
    for (unsigned int j = 0;j<NUMPOINTS;j++) {
        FLT area = M.hexArea*M.regArea(j);
        if (LDn) {
            DchiVal[j] =  Dchi * morph0Area[j] /  area;
            DnVal[j] = Dn * morph0Area[j] / area;
            DcVal[j] = Dc * morph0Area[j] / area;
        }
        else {
            DchiVal[j] =  Dchi;
            DnVal[j] = Dn;
            DcVal[j] = Dc;
        }
        cout << "DnVal region " << j << " = " << DnVal[j] << " Dn = " << Dn << " PI " << PI << endl;
    }
    boundaryCount = 0;
    internalCount = 0;
// now draw the intial tesselation
#ifdef COMPILE_PLOTTING
    if (Lgraphics) {
        morph::Gdisplay ndisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
        ndisp.resetDisplay (fix, eye, rot);
		for (unsigned int j=0;j<NUMPOINTS;j++)
		{
          // plot stuff here.
            cout << "size of Hexgrid j "<< j << " is " << S[j].Hgrid->num() << endl;
            for (auto h : S[j].Hgrid->hexen) {
                if (h.boundaryHex()) {
		      //cout << "h.boundaryHex in region " << j <<  " h.vi " << h.vi << endl;
		     // cout << "region " << j <<  " x " << S[j].Hgrid->d_x[h.vi] << " y " << S[j].Hgrid->d_y[h.vi] << endl;
                ndisp.drawHex (h.position(), (h.d/2.0f), cl_d);
                boundaryCount++;
                }
                else
                {
		      //cout << "h.internal hex in region " << j <<  " h.vi " << h.vi << endl;
		      //cout << "region " << j <<  " x " << S[j].Hgrid->d_x[h.vi] << " y " << S[j].Hgrid->d_y[h.vi] << endl;
//                ndisp.drawHex (h.position(), offset, (h.d/2.0f), cl_b);
                    internalCount++;
                }
            }
            for (unsigned int i=0; i< M.vCoords[j].size();i++) {
                std::array<FLT,3> pos = {M.vCoords[j][i].first, M.vCoords[j][i].second, 0};
                ndisp.drawHex (pos,offset,4.0*hexWidth,cl_d);
                pos = {M.mCoords[j][i].first, M.mCoords[j][i].second, 0};
                ndisp.drawHex (pos,offset,4.0*hexWidth,cl_d);
            }
		//cout << "boundaryCount 2 "<<boundaryCount<< " internalCount 2 "<< internalCount <<endl;


        // Draw small hex at boundary centroid
        array<FLT,3> c;
        c[2] = 0;
        c[0] = S[j].Hgrid->boundaryCentroid.first;
        c[1] = S[j].Hgrid->boundaryCentroid.second;
        cout << "d/2: " << S[j].Hgrid->hexen.begin()->d/4.0f << endl;
        //ndisp.drawHex (c, offset, (S[j].Hgrid->hexen.begin()->d/2.0f), cl_a);
        cout << "boundaryCentroid x,y: " << c[0] << "," << c[1] << endl;

        }//end of loop on NUMPOINTS to plot boundaries
        for (auto h : M.Hgrid->hexen) {
            if (M.Creg[h.vi] ==  1 ) {
                ndisp.drawHex (h.position(), (h.d/2.0f), cl_d);
                boundaryCount++;
            }
            else if (M.Creg[h.vi] > 1) {
                ndisp.drawHex (h.position(), (h.d/2.0f), cl_d);
            }
        }
      // red hex at zero
      //ndisp.drawHex (pos, 0.005, cl_aa);

        usleep (100000);
        cout << "before redrawDisplay 2 " << endl;
        ndisp.redrawDisplay();
        cout << "after redrawDisplay 2" << endl;
        usleep (100000); // one hundred seconds
        ndisp.saveImage(logpath + "/Tesselation3" + iter + ".png");
        ndisp.closeDisplay();
    }
#endif
// initialise the fields
        cout<< "just before third data read"<< " lcontinue " << Lcontinue <<endl;
// initialise with random field
        string hname = logpath + "/third.h5";
        if (Lcontinue) {
             morph::HdfData hinput (hname,1);
             for (unsigned int j=0;j<NUMPOINTS;j++){
                std::string nstr = "n" + to_string(j);
                char * nst = new char[nstr.length()+1];
                std::strcpy(nst,nstr.c_str());
                std::string ccstr = "c" + to_string(j);
                char * ccst = new char[ccstr.length()+1];
                cout << "labels "<< nstr <<" , " << ccstr<<endl;
                std::strcpy(ccst,ccstr.c_str());
                hinput.read_contained_vals(nst,S[j].NN);
                hinput.read_contained_vals(ccst,S[j].CC);
                cout<< "just after input of NN and CC1"<< endl;
       //   input.close();
            }
	 }
     else {
        for (unsigned int j=0;j<NUMPOINTS;j++) {
            for (auto h : S[j].Hgrid->hexen) {
            // boundarySigmoid. Jumps sharply (100, larger is sharper) over length
            // scale 0.05 to 1. So if distance from boundary > 0.05, noise has
            // normal value. Close to boundary, noise is less.
                FLT choice = ruf.get();
                if (choice > 0.5)
                {
                    S[j].NN[h.vi] = - ruf.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = - ruf.get() * aNoiseGain + ccInitialOffset;
                }
                else {
                    S[j].NN[h.vi] = ruf.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = ruf.get() * aNoiseGain + ccInitialOffset;
                }

                if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                    FLT bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary- boundaryFalloffDist)) );
                    S[j].NN[h.vi] = (S[j].NN[h.vi] - nnInitialOffset) * bSig + nnInitialOffset;
                    S[j].CC[h.vi] = (S[j].CC[h.vi] - ccInitialOffset) * bSig + ccInitialOffset;
                } //end of if on boundary distance
            }//end of loop over region
        }//end of loop over all regions
     } //end of else on Lcontinue
     cout <<  "just after field creation" << endl;


     //begin fourth time stepping loop after second morph
    int loopsteps = 0;
    for (int i=0;i<numsteps;i++)
	{
	  //head of fourth time stepping loop only step on internal regions
        for (int j = 0;j<NUMPOINTS;j++) {
            if (M.innerRegion[j]) {
                S[j].step(dt, DnVal[j], DchiVal[j], DcVal[j]);
            }
        }

    #ifdef COMPILE_PLOTTING
	    if ((i%numprint == 0) && Lgraphics) { //only display internal regions
            int countHex = 0;
          //set up display
            morph::Gdisplay nndisp(900, 900, 0, 0, "morph 2 u run", rhoInit, 0.0, 0.0);
            nndisp.resetDisplay (fix, eye, rot);
            for (int j = 0;j<NUMPOINTS;j++) {
                if (M.innerRegion[j]) {
                cout << "in print routine"<<endl;
                unsigned int regsize = S[j].Hgrid->num();
                cout << "in loop over regions " << " size is " << regsize << endl;
                countHex += regsize;
                vector<FLT> regionNN(regsize);
                vector<FLT> tempNN(regsize);
                int k = 0;
                for (auto h : S[j].Hgrid->hexen) {
                    tempNN[k] = S[j].NN[h.vi];
                    k++;
                }
//normalise over the region then write normalised values to normalised NN over hexGrid
                regionNN = L.normalise(tempNN);
                cout << "total number of hexes counted " << countHex << endl;
		      //cout << "just before drawHex 2 "<< h.vi << "normalNN " << normalNN[h.vi] << endl;
                for (auto h : S[j].Hgrid->hexen) {
		      //if (!h.boundaryHex())
			  //{
                    array<FLT,3> colour = morph::Gdisplay::getJetColorF(regionNN[h.vi]);
                    nndisp.drawHex(h.position(),(h.d/2.0f),colour);
		        //cout << "just after drawHex 2"<<endl;
			   //}
			    }
                for (unsigned int i=0; i< M.vCoords[j].size();i++) {
                    std::array<FLT,3> pos = {M.vCoords[j][i].first, M.vCoords[j][i].second, 0};
                    nndisp.drawHex (pos,offset,4.0*hexWidth,cl_d);
                }
                }//end of loop on inner regions
            }//end of loop over regions
            for (auto h : M.Hgrid->hexen) {
                if (M.Creg[h.vi] ==  1 ) {
                    nndisp.drawHex (h.position(), (h.d/2.0f), cl_d);
                    boundaryCount++;
                }
                else if (M.Creg[h.vi] > 1) {
                    nndisp.drawHex (h.position(), (h.d/2.0f), cl_d);
                }
            }
            cout << "just berore redraw display 2 i " << i << endl;
            nndisp.redrawDisplay();
		    cout << "just after redraw display 2 i " <<  i << endl;
            usleep (1000000); // one hundred seconds
		    nndisp.saveImage(logpath + "/nnField3" + iter + ".png");
		    cout << "just after saveImage i " << i  << endl;
            nndisp.closeDisplay();
		    cout << "just after closeDisplay i " << i  << endl;
        }//end of loop on numprint drawing fields
#endif
		cout << "after numprint loop " << i << endl;
        loopsteps = i;
    } // end of fourth time stepping loop
	cout << "after third time stepping loop " << loopsteps << endl;
       //nndisp.closeDisplay();

    //code run at end of timestepping
    cout << "just before the fourth data write" << endl;
    //first save the  ofstream outFile;
    morph::HdfData hdata(hname);
	cout <<"after creating hdata "<< hname << endl;
	for (unsigned int j=0;j<NUMPOINTS;j++) {
		std::string nstr = "n" + to_string(j);
	   	char * nst = new char[nstr.length()+1];
		std::strcpy(nst,nstr.c_str());
		std::string ccstr = "c" + to_string(j);
	    char * ccst = new char[ccstr.length()+1];
		std::strcpy(ccst,ccstr.c_str());
		cout <<"in third hdf5 write "<<nstr <<" , "<<ccstr<<endl;
        hdata.add_contained_vals(ccst,S[j].CC);
        hdata.add_contained_vals(nst,S[j].NN);
	}
    //data.add_contained_vals("X",M.X[0]);
    //data.add_contained_vals("Y",M.X[1]);
    // hdata.add_val ("/Dchi", Dchi);
    // hdata.add_val ("/Dn", Dn);
    // hdata.add_val ("/Dc",Dc);
    cout << " just after writing data "  << endl;

    gfile << endl << "analysis on second morphing iteration " << endl;
    radiusDVector.resize(0);
    angleDVector.resize(0);
    angleVector.resize(0);
    radiusVector.resize(0);
    for (unsigned int j=0;j<NUMPOINTS;j++) {
        M.NN[j] = S[j].NN;
	    if (M.innerRegion[j]) {
	        int regionCount = 0;
            gfile<<"in the degree loop" << endl;
	        //angle degree
            tempArea = M.regArea(j);
            tempPerimeter = M.renewRegPerimeter(j);
            // digital version
            angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, S[j].NN);
            degreeAngle = L.find_zeroDAngle(angleDVector);
            gfile << "region "<< j << " degreeDAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;

           // analogue version
            angleVector = M.sectorize_reg_angle(j,numSectors,radiusOffset, numSectors, S[j].NN);
            angleVector = L.meanzero_vector(angleVector);
            //degreeAngle = M.find_max(angleVector,3);
            degreeAngle = L.find_zeroAngle(angleVector,3);
            gfile << "region "<< j << " degreeAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;
            //radial degree
            degreeRadius = -100;
            radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
            //gfile <<"after sectorize_reg_radius"<<endl;
            degreeRadius = L.find_zeroDRadius(radiusDVector);
            gfile  << "region "<< j << " degreeDRadius "<< degreeRadius << "  " <<endl;
            //radial degree
            degreeRadius = -100;
            int newdegreeRadius = 0;
			for (int angleOffset=0; angleOffset<numSectors -1; angleOffset += 3) {
                radiusVector = M.sectorize_reg_radius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
                newdegreeRadius = L.find_zeroRadius(radiusVector,3);
                if (newdegreeRadius > degreeRadius) {
				    degreeRadius = newdegreeRadius;
                }
		    }
            gfile <<  " region "<< j << " degreeRadius  "<< degreeRadius << "  " <<endl << endl;
            regionCount++;
        } //end of if on non-zero regions
    } //end of loop on NUMPOINT

//computing average values over all regions
    avDegreeAngle = 0;
	avDegreeRadius = 0;
	occupancy = 0;
    avAbsCorrelation = 0;
	tempPerimeter = 0;
	tempArea = 0;
	countRegions = 0;
	M.random_correlate(max_comp, 3);
	cout << "just after random correlate_edges morph2 " << endl;
    radiusDVector.resize(0);
    angleDVector.resize(0);
    angleVector.resize(0);
    radiusVector.resize(0);
    for (unsigned int j=0;j<NUMPOINTS;j++) {
	    if (M.innerRegion[j]){
            countRegions++;
            avAbsCorrelation += M.renewcorrelate_edges(j,3);
            occupancy += M.regNNfrac(j);
            cout << "after renNNfrac " << endl;
            tempArea = M.regArea(j);
            cout << " after regArea " << endl;
            tempPerimeter = M.regPerimeter(j);
            cout << " after regPerimeter " << endl;
            cout << "before sectorixe Dangle morph2 " << endl;
            angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, S[j].NN);
            cout << "after sectorixe Dangle morph2 " << endl;
            degreeAngle = L.find_zeroDAngle(angleDVector);
            avDegreeAngle += degreeAngle;
            cout << "after sectorixe Dangle morph2 " << endl;
            //radial degree
            degreeRadius = 0;
            radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
            degreeRadius = L.find_zeroDRadius(radiusDVector);
            avDegreeRadius += degreeRadius;
            degfile3 << degreeAngle << " " << degreeRadius << " " << M.regNNfrac(j) << " " << tempArea << " " << tempPerimeter<<endl<<flush;
         } //end of if on non-zero regions
    } //end of loop on NUMPOINTs
    if (countRegions == 0) {
        cout << "Error zero regionss counted in fourth analysis morph 2" << endl;
        return -1;
    }
    avDegreeAngle = avDegreeAngle / (1.0 * countRegions);
    avDegreeRadius = avDegreeRadius / (1.0 * countRegions);
	occupancy = occupancy / (1.0 * countRegions);
    avAbsCorrelation = avAbsCorrelation/(1.0 * countRegions);
	jfile <<Dn<<"  "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation  <<endl;
    cout << "end of program reached successfully!" << endl;
    return 0;
};
