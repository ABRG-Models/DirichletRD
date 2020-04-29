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
using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::runtime_error;
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
    if (argc < 9) {
      std::cout << "not enough arguments" << argc << endl;
      return -1;
    }
    //const char* logdir = "cd logs";
   // system (logdir);
    string logpath = argv[1];
    string commandStem;
    bool overwrite_logs = true;
    double dt = stod(argv[2]); //timesetp passed to M.step
    double Dn = stod(argv[3]); //Dn diffusion passed to M.step
    double Dchi = stod(argv[4]); //Dchi chemotaxis passed to M.step
    double Dc = stod(argv[5]);
    int numsteps = atoi(argv[6]); //length of integration 
    // int numprint = atoi(argv[7]); //frequency of printing
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
    ofstream gfile ( logpath + "/edges.out",ios::app);
    ofstream jfile ( logpath + "/results.txt",ios::app);
// seed the random number generator
    unsigned int a_seed = time(0);
	srand(a_seed);

// initialise DRegion class setting scale
    DRegion M(8,logpath);
    cout << "before dissect_boundary " << endl;
    vector<std::pair<double,double>> centroids;
    centroids = M.dissectBoundary(logpath); //dissect region boundary
	M.setRadialSegments(); //set the radial segments for regions
// include the analysis methods
    Analysis L;

    string fname = logpath + "/first.h5";
    cout<< "just before first data read"<< " lcontinue " << Lcontinue <<endl;
// load from HDF5 file or initialise with random field
// as a perturbation of n=1.0, c=2.5
    if (Lcontinue) {
	    morph::HdfData input (fname,1);
	    cout<< "just after first data read"<< endl;
	    input.read_contained_vals("n",M.NN);
	    input.read_contained_vals("c",M.CC);
	    cout<< "just after input of NN and CC1"<< endl;
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

      for (int i=0;i<numsteps;i++) 
	  {
	  //cout << " just before time step " << " i " << i << endl;
         M.step(dt, Dn, Dchi, Dc);
      } //end of numsteps loop 
      //cout << " just after time step i = " << i << endl;

    //code run at end of timestepping
    //first save the  ofstream outFile;
     morph::HdfData data (fname);
     data.add_contained_vals("c",M.CC);
     data.add_contained_vals("n",M.NN);
     data.add_val ("/Dchi", Dchi);
     data.add_val ("/Dn", Dn);
     data.add_val ("/Dc",Dc);
     
     cout << " just after writing data "  << endl;
     //post run analysis

     cout << "before correlate_edges" << endl;
	 double avAbsCorrelation = 0;
     avAbsCorrelation = M.correlate_edges(logpath);
     cout << "after correlate_edges" << endl;
     cout<<"after correlate_edges" << endl;
     int regionCount = 0;
     int numSectors = 12;
	 vector <int> angleDVector;
     vector <int> radiusDVector;
	 vector <double> angleVector;
     vector <double> radiusVector;
     int degreeRadius;
     int degreeAngle;
	 double tempArea = 0;
	 double tempPerimeter = 0;

     for (int j=0;j<NUMPOINTS;j++) {
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
            for (int angleOffset=0; angleOffset<numSectors - 1; angleOffset +=3)
			{
               radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + 2);
               //gfile <<"after sectorize_reg_radius"<<endl;
               // radiusVector = M.meanzero_vector(radiusVector);
               newdegreeRadius = L.find_zeroDRadius(radiusDVector);
               // gfile << " ndR " << newdegreeRadius;
               if (newdegreeRadius > degreeRadius)
                  degreeRadius = newdegreeRadius;
            }
            gfile  << "region "<< j << " degreeDRadius "<< degreeRadius << "  " <<endl ;
            //radial degree
            degreeRadius = -100;
            newdegreeRadius = 0;
            for (int angleOffset=0; angleOffset<numSectors -1; angleOffset += 3)
			{
			   radiusVector = M.sectorize_reg_radius(j,numSectors, angleOffset, angleOffset + 2);
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
    for (auto h : M.Hgrid->hexen) 
	{
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
    for (int i=0;i<iextent;i++)
	{
	   for (int j=0;j<jextent;j++)
	   {
		  NNfield[i][j] = 0;
	   }
    }
    //cout << " after filling NNfield " << endl;
    for (auto h : M.Hgrid->hexen) 
	{
       i = h.ri - h.bi - imin;
       j = h.ri + h.bi - jmin;
       //cout << " i " << i << " j " << j << endl;
       NNfield[i][j] = M.NN[h.vi];
    }
    ofstream hfile (logpath + "/NNfield.txt");
    //cout << " after creating NNfield.txt " << endl;

    for (int i=0;i<iextent;i++)
	{ 
	   for (int j=0;j<jextent;j++)
	   {
		  hfile << setprecision(10) << setw(15) << NNfield[i][j];
	   }
    }


//writing out the results file
	double avDegreeAngle = 0;
	double avDegreeRadius = 0;
	double occupancy = 0;
    tempArea = 0;
    tempPerimeter = 0;
	int countRegions = 0;
    for (int j=0;j<NUMPOINTS;j++) 
	{
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
          for (int angleOffset=0; angleOffset<numSectors -1; angleOffset += 3)
		  {
             radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + 3);
             holdRadius = L.find_zeroDRadius(radiusDVector);
		     //degreeRadius += holdRadius;
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
	// jfile << avDegreeAngle <<" "<<avDegreeRadius<<endl;
	jfile <<Dn<<" "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation<<endl; 
    cout << "just before setting curved boundaries" <<endl;
    // section for solving on the curved boundaries
    // first create the vector of ksSolver classes
    vector<ksSolver> S;
    //set the boundaries for the regions, in curvedBoundary a bezier path
    M.populateBoundVector(1);
    cout << "just after setting curved boundaries " << M.curvedBoundary.size()<<endl;
    //use the curved boundaries to create an array of hexgrids for the new regions
    //and also reset the polar coordinates for the local regions.
    for (int j = 0;j<NUMPOINTS;j++)
	{
	   S.push_back(ksSolver(8,logpath,M.curvedBoundary[j],M.centres[j]));
	   cout << "in the loop populating the ksVector"<< j <<endl;
	}	
	cout << "just after populating the ksVector"<<endl;

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
	        for (int i=0;i<S[j].n;i++) {
		        double choice = morph::Tools::randDouble();
		        if (choice > 0.5)
			        S[j].NN[i]=-(morph::Tools::randDouble())*1.0 + 1.;
		        else
			        S[j].NN[i]=(morph::Tools::randDouble())*1.0 + 1.;
		            S[j].CC[i]=(morph::Tools::randDouble())*1.0 + 2.5;
	      } //end of code to set initial random field
		 }
        } //end of else on Lcontinue
     cout <<  "just after field creation first morph" << endl;

      // begin second time stepping loop
	  for (int i = 0; i<numsteps; i++)
	  {
		for (int j = 0;j<NUMPOINTS;j++) //loop over regions
	    {
		  S[j].step(dt, Dn, Dchi, Dc);
		}
	  }	
	    
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
	  // set the segments from the vertices to the seed points
	  M.setRadialSegments();
	  // repopulate the regions
	  for (int j=0;j<NUMPOINTS;j++)
	  {
	    M.renewRegion(j,S[j].Hgrid->hexen);
	  }
	  //determine the new boundary
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
		cout << " correlation region first morph " << j << " = " << correlation;
	  }
	  */

      gfile << endl << "analysis on first morphing iteration " << endl;
      for (int j=0;j<NUMPOINTS;j++) {
              if (M.regArea(j) != 0){
                  cout<<"in the degree loop" << endl;
                  gfile<<"in the degree loop" << endl;
                  //angle degree
                  tempArea = M.regArea(j)*(5.0/Dn);
                  tempPerimeter = M.renewRegPerimeter(j)*sqrt(5.0/Dn);
                  cout << "after area and perimeter first morphing" << endl;
                  // digital version
                  int radiusOffset = 0;
                  angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, radiusOffset + 11);
				  cout << "after sectorization first morphing" << endl;
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
} //end of loop on NUMPOINT

	      avDegreeAngle = 0;
	      avDegreeRadius = 0;
	      occupancy = 0;
	      countRegions = 0;
		  tempArea = 0;
		  tempPerimeter = 0;
		  avAbsCorrelation = 0;
          for (int j=0;j<NUMPOINTS;j++) {
	        if (M.regArea(j) != 0)
			{
	          countRegions++;
	          gfile << " fraction of positive NN " << M.regNNfrac(j) << endl;
		  occupancy += M.regNNfrac(j);
              tempArea = M.regArea(j)*(5.0/Dn);
              tempPerimeter = M.regPerimeter(j)*sqrt(5.0/Dn);

              avAbsCorrelation += M.renewcorrelate_edges(j,logpath);
              int radiusOffset = 0;
              angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, radiusOffset + 11);
              degreeAngle = L.find_zeroDAngle(angleDVector);
		      avDegreeAngle += degreeAngle;
              //radial degree
              int holdRadius = 0;
		      degreeRadius = 0;
              for (int angleOffset=0; angleOffset<numSectors -1; angleOffset += 3)
			  {
                radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + 2);
                holdRadius = L.find_zeroDRadius(radiusDVector);
		        //degreeRadius += holdRadius;
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
        avAbsCorrelation = avAbsCorrelation / (1.0 * countRegions);
	    occupancy = occupancy / (1.0 * countRegions);
	    // jfile << avDegreeAngle <<" "<<avDegreeRadius<<endl;
	    jfile <<Dn<< " "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation<< " after morphing  " << endl;

// now set the boundaries for the regions, stored in curvedBoundary
		M.populateBoundVector(0);
		S.resize(0);
         cout << "just after setting curved boundaries second iteration" << M.curvedBoundary.size()<<endl;
		for (int j = 0;j<NUMPOINTS;j++)
		{
		  S.push_back(ksSolver(8,logpath,M.curvedBoundary[j],M.centres[j]));
		  cout << "in the loop populating the ksVector"<< j <<endl;
		}	
		cout << "just after populating the ksVector"<<endl;
// initialise the fields
// initialise the fields
    string hname = logpath + "/third.h5";
    cout<< "just before second data read"<< " lcontinue " << Lcontinue <<endl;
// initialise with random field
    if (Lcontinue) 
	{
       morph::HdfData hinput(hname,1);
       cout << "just after trying to open ../logs/third.h5" << endl;
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
	        hinput.read_contained_vals(ccst,S[j].CC);
	        hinput.read_contained_vals(nst,S[j].NN); 
		}
	  }
      else 
	  {
        for (unsigned int j=0;j<NUMPOINTS;j++)
		{
	        for (int i=0;i<S[j].n;i++) {
		        double choice = morph::Tools::randDouble();
		        if (choice > 0.5)
			        S[j].NN[i]=-(morph::Tools::randDouble())*1.0 + 1.;
		        else
			        S[j].NN[i]=(morph::Tools::randDouble())*1.0 + 1.;
		            S[j].CC[i]=(morph::Tools::randDouble())*1.0 + 2.5;
	      } //end of code to set initial random field
		 }
        } //end of else on Lcontinue
     cout <<  "just after field creation first morph" << endl;
       
      // begin third time stepping loop
	  for (int i = 0; i<numsteps; i++)
	  {
		for (int j = 0;j<NUMPOINTS;j++) //loop over regions
	    {
		  S[j].step(dt, Dn, Dchi, Dc);
		}
	  }	

    //code run at end of third timestepping
    //first save the  ofstream outFile;
    morph::HdfData hdata(hname);
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
        hdata.add_contained_vals(ccst,S[j].CC);
        hdata.add_contained_vals(nst,S[j].NN);
        //data.add_contained_vals("X",M.X[0]);
        //data.add_contained_vals("Y",M.X[1]);
     }	
     hdata.add_val ("/Dchi", Dchi);
     hdata.add_val ("/Dn", Dn);
     hdata.add_val ("/Dc",Dc);
     
     cout << " just after writing data "  << endl;
	    
	  // M.setRadialSegments();
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
            cout<<"in the degree loop" << endl;
            gfile<<"in the degree loop" << endl;
            //angle degree
            tempArea = M.regArea(j)*(5.0/Dn);
            tempPerimeter = M.renewRegPerimeter(j)*sqrt(5.0/Dn);
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
            for (int angleOffset=0; angleOffset<numSectors - 1; angleOffset +=3)
			{
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
            for (int angleOffset=0; angleOffset<numSectors -1; angleOffset += 3)
			{
			   radiusVector = M.sectorize_reg_radius(j,numSectors, angleOffset, angleOffset + 3);
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
    avAbsCorrelation = 0;
    tempPerimeter = 0;
    tempArea = 0;
	countRegions = 0;
    for (int j=0;j<NUMPOINTS;j++) 
	{
	   if (M.regArea(j) != 0)
	   {
	      countRegions++;
          avAbsCorrelation += M.renewcorrelate_edges(j,logpath);
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
          for (int angleOffset=0; angleOffset<numSectors -1; angleOffset += 3)
		  {
             radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + 2);
             holdRadius = L.find_zeroDRadius(radiusDVector);
		     //degreeRadius += holdRadius;
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
		avAbsCorrelation = avAbsCorrelation/(1.0 * countRegions);
	    // jfile << avDegreeAngle <<" "<<avDegreeRadius<<endl;
	    jfile <<Dn<<"  "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation<< " after morphing  2" <<endl;
    return 0;
};
