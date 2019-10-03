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
using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::runtime_error;
using morph::HexGrid;
using morph::HdfData;
using morph::Tools;
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
//    commandStem = "mkdir " + logpath;
//    command = commandStem.c_str();
//    system(command);
     //  cerr << "Error : " << strerror(errno) << endl;
    // else
    //   cout << "Directory created" << endl;
	
    //commandStem = "cd " + logpath;
    //command = commandStem.c_str();
    //system(command);
    // ofstream bfile ( logpath + "/maindebug.out" );
    ofstream gfile ( logpath + "/edges.out");
    ofstream jfile ( logpath + "/results.txt");

    // int mdegree = 0;
    vector <double> rayv;
    vector <double> ringv;

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

// initialise Region class setting scale
    Region M(7,logpath);
// include the analysis methods
    Analysis L;

    string fname = logpath + "/fileVal.h5";
    cout<< "just before first data read"<< " lcontinue " << Lcontinue <<endl;
// initialise with random field
    if (Lcontinue) {
	    morph::HdfData input (fname,1);
	    //cout<< "just after first data read"<< endl;
	    input.read_contained_vals("n",M.NN);
	    input.read_contained_vals("c",M.CC);
	    //cout<< "just after input of NN and CC1"<< endl;
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
      for (int i=0;i<numsteps;i++) {
	     cout << " just before time step " << " i " << i << endl;
         M.step(dt, Dn, Dchi, Dc);
         cout << " just after time step i = " << i << endl;

      }
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
