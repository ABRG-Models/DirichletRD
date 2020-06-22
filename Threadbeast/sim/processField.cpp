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
    int Lcontinue = atoi(argv[10]); //logical to determine if coldstart
	int numSectors = 12;
	int numAdjust = 50000000;
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
		    double choice = morph::Tools::randDouble();
            // boundarySigmoid. Jumps sharply (100, larger is sharper) over length
            // scale 0.05 to 1. So if distance from boundary > 0.05, noise has
            // normal value. Close to boundary, noise is less.
		    if (choice > 0.5)
			{
                M.NN[h.vi] = -(morph::Tools::randDouble()) * aNoiseGain +nnInitialOffset;
                M.CC[h.vi] = -(morph::Tools::randDouble()) * aNoiseGain + ccInitialOffset;
			}
			else
			{
                M.NN[h.vi] = (morph::Tools::randDouble()) * aNoiseGain +nnInitialOffset;
                M.CC[h.vi] = (morph::Tools::randDouble()) * aNoiseGain + ccInitialOffset;
			}
            if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                double bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary- boundaryFalloffDist)) );
                M.NN[h.vi] = (M.NN[h.vi] - nnInitialOffset) * bSig + nnInitialOffset;
              //  M.CC[h.vi] = (M.CC[h.vi] - ccInitialOffset) * bSig + ccInitialOffset;
		    } //end of if on boundary distance
	    }//end of loop over region
    } //end of else on Lcontinue
    cout <<  "just after field creation" << endl;


    vector<double> fix(15, 0.0);
    vector<double> eye(15, 0.0);
    vector<double> rot(3, 0.0);
    double rhoInit = 5.1;
    morph::Gdisplay isp(900, 900, 0, 0, "A boundary", rhoInit, 0.0, 0.0);
    isp.resetDisplay (fix, eye, rot);
//plot stuff here.
    array<float,3> cl_a = morph::Tools::getJetColorF (0.78);
    array<float,3> cl_c = morph::Tools::getJetColorF (0.28);
    array<float,3> cl_b = morph::Tools::getJetColorF (0.58);
    array<float,3> offset = {{0, 0, 0}};
    int boundaryCount = 0;
    for (auto h : M.Hgrid->hexen) {
        if (M.Creg[h.vi] ==  1 ) {
            isp.drawHex (h.position(), (h.d/2.0f), cl_a);
            boundaryCount++;
		}
		else if (M.Creg[h.vi] > 1) {
            isp.drawHex (h.position(), (h.d/2.0f), cl_c);
        }
		else
		{
            isp.drawHex (h.position(), offset, (h.d/2.0f), cl_b);
        }
    }
	cout << "boundaryCount "<<boundaryCount<<endl;



    usleep (10000000);
    isp.redrawDisplay();
    usleep (10000000); // one hundred seconds
    isp.saveImage(logpath + "/Tesselation.png");
    isp.closeDisplay();
    float hexWidth = M.Hgrid->hexen.begin()->d/2.0;
    cerr << "d/2: " << hexWidth << endl;
    for (int i=0;i<numsteps;i++)
	{
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
		    } //end of loop over regions
		    cout << "total number of hexes counted " << countHex << endl;
		    for (auto h : M.Hgrid->hexen) {
		         //cout << "just before drawHex  "<< h.vi << "normalNN " << normalNN[h.vi] << endl;
		             if (M.Creg[h.vi] == 0) {
			             if (h.d/2.0 > 0.05){
			                 cerr << "hex too large for " << h.vi <<endl;
					     }
			             array<float,3> colour = morph::Tools::getJetColorF(normalNN[h.vi]);
	                     disp.drawHex(h.position(),hexWidth,colour);
//  cout << "just after drawHex"<<endl;
 			         }
		             else {
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
        data.add_contained_vals("c",M.CC);
        data.add_contained_vals("n",M.NN);
        data.add_val ("/Dchi", Dchi);
        data.add_val ("/Dn", Dn);
        data.add_val ("/Dc",Dc);

        cout << " just after writing data "  << endl;
     //post run analysis

        cout << "before correlate_edges" << endl;
	    double avAbsCorrelation = 0;
        avAbsCorrelation = M.correlate_edges();
// look at correlation between random edges
         const int max_comp = NUMPOINTS*6;
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
//writing out of the image files
    int imin = 1000000;
    int imax = -1000000;
    int jmin =  1000000;
    int jmax = -1000000;
    int i,j;
    for (auto h : M.Hgrid->hexen) {
       i = h.ri -h.bi;
       j = h.ri + h.bi;
       if (i<imin) imin = i;
       if (i>imax) imax = i;
       if (j<jmin) jmin = j;
       if (j>jmax) jmax = j;
    }
    int iextent = imax - imin + 1;
    int jextent = jmax - jmin + 1;
    //cout << " imax " << imax << " imin " << imin << " jmax " << jmax << " jmin " << jmin << endl;
    //cout << " iextent " << iextent << " jextent " << jextent << " gridsize " << M.n << endl;
    double NNfield[iextent][jextent];
    for (int i=0;i<iextent;i++){
	    for (int j=0;j<jextent;j++){
		    NNfield[i][j] = 0;
	    }
    }
    //cout << " after filling NNfield " << endl;
    for (auto h : M.Hgrid->hexen) {
       i = h.ri - h.bi - imin;
       j = h.ri + h.bi - jmin;
       //cout << " i " << i << " j " << j << endl;
       NNfield[i][j] = M.NN[h.vi];
    }
    ofstream hfile (logpath + "/NNfield.txt");
    //cout << " after creating NNfield.txt " << endl;

    for (int i=0;i<iextent;i++){
	    for (int j=0;j<jextent;j++){
		    hfile << setprecision(10) << setw(15) << NNfield[i][j];
	    }
    }


    gfile << " imin " << imin << " imax " << imax << " jmin " << jmin << " jmax " << jmax << endl;

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
        morph::Gdisplay mdisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
        mdisp.resetDisplay (fix, eye, rot);
		for (int j=0;j<NUMPOINTS;j++)
		{
          // plot stuff here.
		  int boundaryCount = 0;
		  int internalCount = 0;
		  cout << "size of Hexgrid j "<< j << " is " << S[j].Hgrid->num() << endl;
          for (auto h : S[j].Hgrid->hexen) {
            if (h.boundaryHex())
			{
		      cout << "h.boundaryHex in region " << j <<  " h.vi " << h.vi << endl;
		     // cout << "region " << j <<  " x " << S[j].Hgrid->d_x[h.vi] << " y " << S[j].Hgrid->d_y[h.vi] << endl;
              mdisp.drawHex (h.position(), (h.d/2.0f), cl_a);
			  boundaryCount++;
			}
            else
			{
		      cout << "h.internal hex in region " << j <<  " h.vi " << h.vi << endl;
		      //cout << "region " << j <<  " x " << S[j].Hgrid->d_x[h.vi] << " y " << S[j].Hgrid->d_y[h.vi] << endl;
                mdisp.drawHex (h.position(), offset, (h.d/2.0f), cl_b);
				internalCount++;
            }
        }
		//cout << "boundaryCount 2 "<<boundaryCount<< " internalCount 2 "<< internalCount <<endl;


      }//end of loop on NUMPOINTS to plot boundaries
      // red hex at zero
      //mdisp.drawHex (pos, 0.0005, cl_aa);

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
                    S[j].NN[h.vi] = -(morph::Tools::randDouble()) * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = -(morph::Tools::randDouble()) * aNoiseGain + ccInitialOffset;
				}
		        else
				{
                    S[j].NN[h.vi] = (morph::Tools::randDouble()) * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = (morph::Tools::randDouble()) * aNoiseGain + ccInitialOffset;
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
	  int countHex = 0;
   	  for (int j = 0;j<NUMPOINTS;j++) //loop over regions
	  {
		  //cout << "just before morph 1 step i "<< i <<endl;
     	  S[j].step(dt, Dn, Dchi, Dc);
		  if (i%100 == 0)
		     cout << "just after morph 1 octopus i "<< j << " NN[5] " << S[j].NN[5] << endl;
	  }
	  if (i%numprint == 0)
	  {
          //set up display
          morph::Gdisplay mmdisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
          mmdisp.resetDisplay (fix, eye, rot);
   	      for (int j = 0;j<NUMPOINTS;j++) //loop over regions
	      {
	 	    cout << "in print routine"<<endl;
            unsigned int regsize = S[j].Hgrid->num();
		    cout << "in loop over regions " << " size is " << regsize << endl;
			countHex += regsize;
		    vector<double> regionNN;
			vector<double> tempNN;
		    //for (auto h : S[j].Hgrid->hexen)
			for (int k=0;k<regsize;k++)
			{
			  tempNN.push_back(S[j].NN[k]);
			  cout << "tempNN " << S[j].NN[k] << " k " << k  << endl;
			}
//normalise over the region then write normalised values to normalised NN over hexGrid
            regionNN = L.normalise(tempNN);
		    cout << "total number of hexes counted " << countHex << endl;
			int idx = 0;
		    for (auto &h : S[j].Hgrid->hexen)
		    //for (int k=0;k<regsize;k++)
			{
		      //if (!h.boundaryHex())
			    array<float,3> colour = morph::Tools::getJetColorF(regionNN[idx]);
	            mmdisp.drawHex(h.position(),(h.d/2.0f),colour);
		        cout << "just after drawHex 2 colour "<< " value " << regionNN[idx] << endl;
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
		   cout << "just after close display 1 i " << i << endl;
		 }//end of loop on numprint drawing fields
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


    gfile << " imin " << imin << " imax " << imax << " jmin " << jmin << " jmax " << jmax << endl;

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
        morph::Gdisplay ndisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
        ndisp.resetDisplay (fix, eye, rot);
		for (int j=0;j<NUMPOINTS;j++)
		{
          // plot stuff here.
		  int boundaryCount = 0;
		  int internalCount = 0;
		  cout << "size of Hexgrid j "<< j << " is " << S[j].Hgrid->num() << endl;
          for (auto h : S[j].Hgrid->hexen) {
            if (h.boundaryHex())
			{
		      //cout << "h.boundaryHex in region " << j <<  " h.vi " << h.vi << endl;
		     // cout << "region " << j <<  " x " << S[j].Hgrid->d_x[h.vi] << " y " << S[j].Hgrid->d_y[h.vi] << endl;
              ndisp.drawHex (h.position(), (h.d/2.0f), cl_a);
			  boundaryCount++;
			}
            else
			{
		      //cout << "h.internal hex in region " << j <<  " h.vi " << h.vi << endl;
		      //cout << "region " << j <<  " x " << S[j].Hgrid->d_x[h.vi] << " y " << S[j].Hgrid->d_y[h.vi] << endl;
                ndisp.drawHex (h.position(), offset, (h.d/2.0f), cl_b);
				internalCount++;
            }
        }
		//cout << "boundaryCount 2 "<<boundaryCount<< " internalCount 2 "<< internalCount <<endl;


        // Draw small hex at boundary centroid
        array<float,3> c;
        c[2] = 0;
        c[0] = S[j].Hgrid->boundaryCentroid.first;
        c[1] = S[j].Hgrid->boundaryCentroid.second;
        cout << "d/2: " << S[j].Hgrid->hexen.begin()->d/4.0f << endl;
        //ndisp.drawHex (c, offset, (S[j].Hgrid->hexen.begin()->d/2.0f), cl_a);
        cout << "boundaryCentroid x,y: " << c[0] << "," << c[1] << endl;

      }//end of loop on NUMPOINTS to plot boundaries
      // red hex at zero
      //ndisp.drawHex (pos, 0.005, cl_aa);

      usleep (100000);
	  cout << "before redrawDisplay 2 " << endl;
      ndisp.redrawDisplay();
	  cout << "after redrawDisplay 2" << endl;
      usleep (100000); // one hundred seconds
      ndisp.saveImage(logpath + "/Tesselation3.png");
      ndisp.closeDisplay();
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
			 double choice = morph::Tools::randDouble();
		     if (choice > 0.5)
			 {
                  S[j].NN[h.vi] = -(morph::Tools::randDouble()) * aNoiseGain +nnInitialOffset;
                  S[j].CC[h.vi] = -(morph::Tools::randDouble()) * aNoiseGain + ccInitialOffset;
			 }
		     else
		   	 {
                  S[j].NN[h.vi] = (morph::Tools::randDouble()) * aNoiseGain +nnInitialOffset;
                  S[j].CC[h.vi] = (morph::Tools::randDouble()) * aNoiseGain + ccInitialOffset;
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

     //morph::Gdisplay nndisp(600, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
     //nndisp.resetDisplay (fix, eye, rot);

     //begin second time stepping loop
    for (i=0;i<numsteps;i++)
	{
	  //cout << " head of second time stepping loop i " << i << endl;
   	  for (int j = 0;j<NUMPOINTS;j++) //loop over regions
	  {
     	   S[j].step(dt, Dn, Dchi, Dc);
     	 //  S[j].step(dt, Dn, Dchi, Dc);
	  }

	  if (i%numprint == 0)
	  {
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
		   cout << "just berore redraw display 2 i " << i << endl;
           nndisp.redrawDisplay();
		   cout << "just after redraw display 2 i " <<  i << endl;
           usleep (1000000); // one hundred seconds
		   nndisp.saveImage(logpath + "/nnField3.png");
		   cout << "just after saveImage i " << i  << endl;
           nndisp.closeDisplay();
		   cout << "just after closeDisplay i " << i  << endl;
		 }//end of loop on numprint drawing fields
		 cout << "after numprint loop " << i << endl;
       } // end of second time stepping
	   cout << "after second time stepping loop " << i << endl;
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
