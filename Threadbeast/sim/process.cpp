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
#include "analysis.h"
#include "region.h"
#include <morph/display.h>
using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::runtime_error;
/*
using morph::HexGrid;
using morph::HdfData;
using morph::Tools;
using morph::Display
*/
using namespace morph;
using namespace std;


int main (int argc, char **argv)
{
    if (argc < 9) {
      std::cout << "not enough arguments" << argc << endl;
      return -1;
    }
    //const char* logdir = "cd logs";
   // system (logdir);
    string logpath = argv[1];
    string commandStem;
    bool overwrite_logs = true;
    const char* command;
    double dt = stod(argv[2]); //timesetp passed to M.step
    double Dn = stod(argv[3]); //Dn diffusion passed to M.step
    double Dchi = stod(argv[4]); //Dchi chemotaxis passed to M.step
    double Dc = stod(argv[5]);
    int numsteps = atoi(argv[6]); //length of integration 
    int numprint = atoi(argv[7]); //frequency of printing
    int Lcontinue = atoi(argv[8]); //logical to determine if coldstart
      /*
     * Now create a log directory if necessary, and exit on any
     * failures.
     */
    if (morph::Tools::dirExists (logpath) == false) {
        morph::Tools::createDir (logpath);
        if (morph::Tools::dirExists (logpath) == false) {
            cerr << "Failed to create the logpath directory "
                 << logpath << " which does not exist."<< endl;
            return 1;
        }
    } else {
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


    // DISPLAYS
    //vector<morph::Gdisplay> displays;
    //vector<double>fix(3,0.0);
    //vector<double>eye(3,0.0);
    //vector<double>rot(3,0.0);
    //displays.push_back (morph::Gdisplay (600, "NN field", 0., 0., 0.));
    //displays[0].resetDisplay(fix,eye,rot);
    //displays[0].redrawDisplay();
    //displays.push_back (morph::Gdisplay (600, "CC field", 0., 0., 0.));
    //displays[1].resetDisplay(fix,eye,rot);
    //displays[1].redrawDisplay();
    //bfile << "just after displays" << endl;

// initialise DRegion class setting scale
    DRegion M(9,logpath);
// include the analysis methods
    Analysis L;

    string fname = logpath + "/fileVal.h5";
    cout<< "just before first data read"<< " lcontinue " << Lcontinue <<endl;
// initialise with random field
    if (Lcontinue) {
	    morph::HdfData input (fname,1);
	    cout<< "just after first data read"<< endl;
	    input.read_contained_vals("n",M.NN);
	    input.read_contained_vals("c",M.CC);
	    cout<< "just after input of NN and CC1"<< endl;
//	    input.close();
    }
    else {
	    for (int i=0;i<M.n;i++) {
		    double choice = morph::Tools::randDouble();
		    if (choice > 0.5)
			    M.NN[i]=-(morph::Tools::randDouble())*1.0 + 1.;
		    else
			    M.NN[i]=(morph::Tools::randDouble())*1.0 + 1.;
		    M.CC[i]=(morph::Tools::randDouble())*1.0 + 2.5;
	    } //end of code to set initial random field
    } //end of else on Lcontinue
    cout <<  "just after field creation" << endl;

        vector<double> fix(3, 0.0);
        vector<double> eye(3, 0.0);
        vector<double> rot(3, 0.0);
        double rhoInit = 1.7;
        morph::Gdisplay disp(1200, 1800, 0, 0, "A boundary", rhoInit, 0.0, 0.0);
        disp.resetDisplay (fix, eye, rot);

        // plot stuff here.
        array<float,3> cl_a = morph::Tools::getJetColorF (0.78);
        array<float,3> cl_c = morph::Tools::getJetColorF (0.28);
        array<float,3> cl_b = morph::Tools::getJetColorF (0.58);
        array<float,3> offset = {{0, 0, 0}};
        for (auto h : M.Hgrid->hexen) {
            if (M.Creg[h.vi] ==  1 ) {
                disp.drawHex (h.position(), (h.d/2.0f), cl_a);
			} else if (M.Creg[h.vi] > 1) {
                disp.drawHex (h.position(), (h.d/2.0f), cl_c);
            } else {
                disp.drawHex (h.position(), offset, (h.d/2.0f), cl_b);
            }
        }


        // Draw small hex at boundary centroid
        array<float,3> c;
        c[2] = 0;
        c[0] = M.Hgrid->boundaryCentroid.first;
        c[1] = M.Hgrid->boundaryCentroid.second;
        cout << "d/2: " << M.Hgrid->hexen.begin()->d/4.0f << endl;
        disp.drawHex (c, offset, (M.Hgrid->hexen.begin()->d/2.0f), cl_a);
        cout << "boundaryCentroid x,y: " << c[0] << "," << c[1] << endl;

        // red hex at zero
        array<float,3> cl_aa = morph::Tools::getJetColorF (0.98);
        array<float,3> pos = { { 0, 0, 0} };
        disp.drawHex (pos, 0.05, cl_aa);

      //  usleep (1000000);
        disp.redrawDisplay();
      unsigned int sleep_seconds = 1;
      while (sleep_seconds--) {
         usleep (1000000); // one hundred seconds
      }
      disp.closeDisplay();

      int numcolour = 0;
      for (int i=0;i<numsteps;i++) {
	     //cout << " just before time step " << " i " << i << endl;
         M.step(dt, Dn, Dchi, Dc);
		 /*
		 if (i % numprint == numprint/2) {
          disp.closeDisplay();
	      usleep(1000000); 
		  }
		  */
		 if (i % numprint == 0) {
//          numcolour++;
          morph::Gdisplay disp(400, 800, 0, 0, "A boundary", rhoInit, 0.0, 0.0);
          disp.resetDisplay (fix, eye, rot);
		  cout << "in print routine"<<endl;
		  vector<double> normalNN;
		  normalNN.resize(M.n);
          ofstream afile (logpath + "/NNcompare.txt");
          afile << " NN " << "                normalNN" << endl;
		  int countHex = 0;
		  for (int j=0;j<NUMPOINTS;j++) {
            unsigned int regsize = M.regionIndex[j].size();
		    cout << "in loop over regions " << j << " size is " << regsize << endl;
			countHex += regsize;
		    vector<double> regionNN(regsize);
			vector<double> tempNN(regsize);
			vector<int> regionIdx(regsize);
		    for (unsigned int k=0;k<regsize;k++){
			  int index = M.regionIndex[j][k];
			  tempNN[k] = M.NN[index];
			  regionIdx[k] = index;
			  }
		   //normalise over the region then write normalised values to normalised NN over hexGrid
           regionNN = L.normalise(tempNN);
		   for (unsigned int k=0;k<regsize;k++) {
		     normalNN[regionIdx[k]] = regionNN[k];
             afile << M.NN[regionIdx[k]] << "    " << normalNN[regionIdx[k]] << " " << regionIdx[k] <<endl;
             //afile << regionIdx[k] <<endl;
			 }
		   } //end of loop over regions
		   cout << "total number of hexes counted " << countHex << endl;
		   for (auto h : M.Hgrid->hexen) {
		     cout << "just before drawHex  "<< h.vi << "normalNN " << normalNN[h.vi] << endl;
		     if (M.Cnbr[h.vi] == 6) {
			   array<float,3> colour = morph::Tools::getJetColorF(normalNN[h.vi]);
			  // array<float,3> colour = morph::Tools::getJetColorF(numcolour*0.1);
	           disp.drawHex(h.position(),(h.d/2.0f),colour);
		       cout << "just after drawHex"<<endl;
			  }
			}
		    // else {
			//  disp.drawHex(h.position(),(h.d/2.0f),0.0f);
			// }
        //disp.resetDisplay (fix, eye, rot);
		 cout << "just berore redraw display" << endl;
          disp.redrawDisplay();
          usleep (10000000); // one hundred seconds
		  if (i == (numsteps - numprint)) {
		    disp.saveImage(logpath + "/nnField.png");
			}
          disp.closeDisplay();
         } // end of print on numprint
        } //end of numsteps loop 
         //cout << " just after time step i = " << i << endl;

    //code run at end of timestepping
    //first save the  ofstream outFile;
     morph::HdfData data (fname);
     data.add_contained_vals("c",M.CC);
     data.add_contained_vals("n",M.NN);
     //data.add_contained_vals("X",M.X[0]);
     //data.add_contained_vals("Y",M.X[1]);
     data.add_val ("/Dchi", Dchi);
     data.add_val ("/Dn", Dn);
     data.add_val ("/Dc",Dc);
     
     cout << " just after writing data "  << endl;
     //post run analysis


            vector<std::pair<double,double>> centroids;
            cout << "before dissect_boundary " << endl;
            centroids = M.dissectBoundary(logpath);
            cout << "before correlate_edges" << endl;
	    double avAbsCorrelation;
            avAbsCorrelation = M.correlate_edges(logpath);
            cout << "after correlate_edges" << endl;
          //  cout<<"after correlate_edges" << endl;
            int regionCount = 0;
            int numSectors = 12;
	    vector <int> angleDVector;
            vector <int> radiusDVector;
	    vector <double> angleVector;
            vector <double> radiusVector;
            int degreeRadius;
            int degreeAngle;
	    double tempArea;
	    double tempPerimeter;

for (int j=0;j<NUMPOINTS-1;j++) {
              if (M.regArea(j) != 0){
                  cout<<"in the degree loop" << endl;
                  gfile<<"in the degree loop" << endl;
                  //angle degree
                  tempArea = M.regArea(j)*(5.0/Dn);
                  tempPerimeter = M.regPerimeter(j)*sqrt(5.0/Dn);

                  // digital version
                  int radiusOffset = 0;
                  angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, radiusOffset + 11);
                  degreeAngle = L.find_zeroDAngle(angleDVector);
                  gfile << "region "<< j << " degreeDAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;

                  // analogue version
                  radiusOffset = 0;
                  angleVector = M.sectorize_reg_angle(j,numSectors,radiusOffset, radiusOffset + 11);
                  angleVector = L.meanzero_vector(angleVector);
                  //degreeAngle = M.find_max(angleVector,3);
                  degreeAngle = L.find_zeroAngle(angleVector,3);
                  gfile << "region "<< j << " degreeAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;
                  //radial degree
                  degreeRadius = -100;
                  int newdegreeRadius = 0;
                  for (int angleOffset=0; angleOffset<numSectors - 1; angleOffset +=3){
                  radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + 3);
                  //gfile <<"after sectorize_reg_radius"<<endl;
                  // radiusVector = M.meanzero_vector(radiusVector);

                  newdegreeRadius = L.find_zeroDRadius(radiusDVector);
                  // gfile << " ndR " << newdegreeRadius;
                  if (newdegreeRadius > degreeRadius)
                          degreeRadius = newdegreeRadius;
                  }
                  gfile  << "region "<< j << " degreeDRadius "<< degreeRadius << "  " <<endl ;
		  
                  ///radial degree
                  degreeRadius = -100;
                  newdegreeRadius = 0;
                  for (int angleOffset=0; angleOffset<numSectors -1; angleOffset += 3){
			  radiusVector = M.sectorize_reg_radius(j,numSectors, angleOffset, angleOffset + 3);
			  newdegreeRadius = L.find_zeroRadius(radiusVector,3);
			  if (newdegreeRadius > degreeRadius)
				  degreeRadius = newdegreeRadius;
		  }


                  gfile <<  " region "<< j << " degreeRadius  "<< degreeRadius << "  " <<endl << endl;

//                  W.logfile <<" degreeRadius "<< degreeRadius<<" degreeAngle "<< degreeAngle << " " << tempArea<<"  "<<tempPerimeter<<endl<<flush;

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
	      double totalDegree = 0;
	      int countRegions = 0;
              for (int j=0;j<NUMPOINTS-1;j++) {
	      if (M.regArea(j) != 0){
	          countRegions++;
	          gfile << " fraction of positive NN " << M.regNNfrac(j) << endl;
		  occupancy += M.regNNfrac(j);
                  tempArea = M.regArea(j)*(5.0/Dn);
                  tempPerimeter = M.regPerimeter(j)*sqrt(5.0/Dn);

                  int radiusOffset = 0;
                  angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, radiusOffset + 11);
                  degreeAngle = L.find_zeroDAngle(angleDVector);
		  avDegreeAngle += degreeAngle;
                  //radial degree
                  int holdRadius = 0;
		  degreeRadius = 0;
		  int numiters = 0;
                  for (int angleOffset=0; angleOffset<numSectors -1; angleOffset += 3){
                  radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + 3);
                  holdRadius = L.find_zeroDRadius(radiusDVector);
		  //degreeRadius += holdRadius;
		  //numiters++;
		  if (holdRadius > degreeRadius)
			  degreeRadius = holdRadius;
                  }
                  avDegreeRadius += degreeRadius;

                  gfile << " degreeRadius "<< degreeRadius<<" degreeAngle "<< degreeAngle << " " << tempArea<<"  "<<tempPerimeter<<endl<<flush;
                  regionCount++;
	      } //end of if on non-zero regions
	      } //end of loop on NUMPOINTs
              avDegreeAngle = avDegreeAngle / (1.0 * countRegions);
              avDegreeRadius = avDegreeRadius / (1.0 * countRegions);
	      occupancy = occupancy / (1.0 * countRegions);
	      totalDegree = avDegreeAngle + 4.0*avDegreeRadius;
	     // jfile << avDegreeAngle <<" "<<avDegreeRadius<<endl;
	      jfile <<Dn<<" "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation<<endl;


    return 0;
};
