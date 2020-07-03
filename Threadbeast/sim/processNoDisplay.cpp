/*
#include <morph/tools.h>
#include <morph/HexGrid.h>
#include <morph/HdfData.h>
#include <morph/Random.h>
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
#include <morph/Random.h>
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
    int Lcontinue = atoi(argv[9]); //logical to determine if coldstart

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

    unsigned int seed = time(NULL);

    // A rando2yym uniform generator returning real/floating point types
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
    ofstream jfile ( logpath + "/results.txt",ios::app);
	ofstream degfile1 (logpath + "/degree1.data",ios::app);
	ofstream degfile2 (logpath + "/degree2.data");
	ofstream degfile3 (logpath + "/degree3.data");


// initialise DRegion class setting scale
    DRegion M(8,5.0,logpath);
    cout << "before dissect_boundary " << endl;
	for (auto h : M.Hgrid->hexen) {
	   cout << " first h.vi output " << h.vi << endl;
	   }
    vector<std::pair<double,double>> centroids;
    centroids = M.dissectBoundary(); //dissect region boundary
	M.setRadialSegments(); //set the radial segments for regions
// include the analysis methods
    Analysis L;

    string fname = logpath + "/first.h5";
    cout<< "just before first data read"<< " Lcontinue " << Lcontinue <<endl;
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
		    double choice = ruf.get();
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
                M.NN[h.vi] = ruf.get() * aNoiseGain +nnInitialOffset;
                M.CC[h.vi] = ruf.get() * aNoiseGain + ccInitialOffset;
			}
            if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                double bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary- boundaryFalloffDist)) );
                M.NN[h.vi] = (M.NN[h.vi] - nnInitialOffset) * bSig + nnInitialOffset;
              //  M.CC[h.vi] = (M.CC[h.vi] - ccInitialOffset) * bSig + ccInitialOffset;
		    } //end of if on boundary distance
	    }//end of loop over region
    } //end of else on Lcontinue
    cout <<  "just after field creation" << endl;

    float hexWidth = M.Hgrid->hexen.begin()->d/2.0;
    cerr << "d/2: " << hexWidth << endl;
    int stepped = 0;
    for (int i=0;i<numsteps;i++)
	{
//cout << " just before time step " << " i " << i << endl;
        M.step(dt, Dn, Dchi, Dc);
        stepped++;
	} //end of numsteps loop
    cout << " just after step loop  iterations  = " << stepped << endl;

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
// look at correlation between random edges
         const int max_comp = NUMPOINTS*5;
         M.random_correlate(max_comp,0);
         cout << "after correlate_edges" << endl;
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

               degfile1 << degreeAngle/2 << " " << degreeRadius << " " << M.regNNfrac(j) << " " << tempArea<< " "<< tempPerimeter<<endl<<flush;
		  }
	    } //end of loop on NUMPOINTs
        avDegreeAngle = avDegreeAngle / (1.0 * countRegions);
        avDegreeRadius = avDegreeRadius / (1.0 * countRegions);
	    occupancy = occupancy / (1.0 * countRegions);
	    jfile <<Dn<<" "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation<<endl;
//end of integration on the original tesselation
         if (skipMorph) return 0;
         cout << "just before setting curved boundaries" <<endl;
// section for solving on the curved boundaries
         vector<ksSolver> S;
		//S.resize(NUMPOINTS);
// now set the boundaries for the regions, stored in curvedBoundary
		M.populateBoundVector(1);
         cout << "just after setting curved boundaries " << M.curvedBoundary.size()<<endl;
		for (int j = 0;j<NUMPOINTS;j++)
		{
		  S.push_back(ksSolver(scale, xspan, logpath, M.curvedBoundary[j], M.centres[j]));
		  cout << "in the loop populating the ksVector"<< j <<endl;
		}
		cout << "just after populating the ksVector"<<endl;
		for (int j = 0;j<NUMPOINTS;j++) {
		   for (auto h : S[j].Hgrid->hexen) {
		      cout << "hexNumber in region " << j <<  " h.vi " << h.vi << endl;
			}
	     }
// now draw the intial tesselation
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
		        double choice = ruf.get ();
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

    // begin second time stepping loop
    for (int i=0;i<numsteps;i++)
    {
   	  for (int j = 0;j<NUMPOINTS;j++) //loop over regions
	  {
		  //cout << "just before morph 1 step i "<< i <<endl;
     	  S[j].step(dt, Dn, Dchi, Dc);
	  }
		 //cout << "iteration " << i << " of morph 1 time step" << endl;
     } // end of second time stepping

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
	 M.setRadialSegments();
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
	        if (M.regArea(j) != 0)
			{
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
	       } //end of if on non-zero regions
	    } //end of loop on NUMPOINTs
        avDegreeAngle = avDegreeAngle / (1.0 * countRegions);
        avDegreeRadius = avDegreeRadius / (1.0 * countRegions);
        avAbsCorrelation = avAbsCorrelation / (1.0 * countRegions);
	    occupancy = occupancy / (1.0 * countRegions);
	    // jfile << avDegreeAngle <<" "<<avDegreeRadius<<endl;
	    jfile <<Dn<< " "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation<< " after morphing  " << endl;
//end of integration after first morphing
//begin second morphing
// now set the boundaries for the regions, stored in curvedBoundary
		M.populateBoundVector(0);
		S.resize(0);
         cout << "just after setting curved boundaries second iteration" << M.curvedBoundary.size()<<endl;
		for (int j = 0;j<NUMPOINTS;j++)
		{
		  S.push_back(ksSolver(scale, xspan, logpath, M.curvedBoundary[j], M.centres[j]));
		  cout << "in the loop populating the ksVector"<< j <<endl;
		}
		cout << "just after populating the ksVector"<<endl;
// now draw the intial tesselation
// initialise the fields
    cout<< "just before first data read"<< " lcontinue " << Lcontinue <<endl;
// initialise with random field
    string hname = logpath + "/third.h5";
    if (Lcontinue)
	{
	   morph::HdfData hinput (hname,1);
       for (int j=0;j<NUMPOINTS;j++)
    	{
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
     else
	 {
       for (int j=0;j<NUMPOINTS;j++)
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
    cout <<  "just after field creation" << endl;


     //begin second time stepping loop
    int loopsteps = 0;
    for (int i=0;i<numsteps;i++)
    {
	  //cout << " head of second time stepping loop i " << i << endl;
   	  for (int j = 0;j<NUMPOINTS;j++) //loop over regions
	  {
     	   S[j].step(dt, Dn, Dchi, Dc);
     	 //  S[j].step(dt, Dn, Dchi, Dc);
	  }

     } // end of second time stepping
	   cout << "after second time stepping loop " << loopsteps << endl;
       //nndisp.closeDisplay();

    //code run at end of timestepping
    //first save the  ofstream outFile;
     morph::HdfData hdata(hname);
	 cout <<"after creating hdata "<< hname << endl;
	 for (int j=0;j<NUMPOINTS;j++)
	 {
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
	 M.setRadialSegments();
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
	      if (M.regArea(j) != 0)
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
	  cout << "just after renewcorrelate_edges morph2 " << endl;
	  M.random_correlate(max_comp, 2);
	  cout << "just after random correlate_edges morph2 " << endl;
      for (int j=0;j<NUMPOINTS;j++)
	  {
	      if (M.regArea(j) != 0)
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
              //radial degree
		      degreeRadius = 0;
              radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2);
              degreeRadius = L.find_zeroDRadius(radiusDVector);
              avDegreeRadius += degreeRadius;

              degfile3 << degreeAngle << " " << degreeRadius << " " << M.regNNfrac(j) << " " << tempArea << " " << tempPerimeter<<endl<<flush;
	      } //end of if on non-zero regions
	    } //end of loop on NUMPOINTs
        avDegreeAngle = avDegreeAngle / (1.0 * countRegions);
        avDegreeRadius = avDegreeRadius / (1.0 * countRegions);
	    occupancy = occupancy / (1.0 * countRegions);
		avAbsCorrelation = avAbsCorrelation/(1.0 * countRegions);
	    // jfile << avDegreeAngle <<" "<<avDegreeRadius<<endl;
	    jfile <<Dn<<"  "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation<< " after morphing  2" <<endl;
    return 0;
};
