/* 
 * ksSolver class
 * Author: John Brooke
 *
 * Date 2019/10
 *
 * creates a hexGrid given a boundary curve and solves 
 * the KS equations
 *
 */
#include <morph/tools.h>
#include <morph/HexGrid.h>
#include <morph/ReadCurves.h>
#include <morph/HdfData.h>
#include <morph/BezCurve.h>
#include <morph/BezCurvePath.h>
#include <morph/BezCoord.h>
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
//#include "hexGeometry.h"
// #include <boost/math/special_functions/bessel.hpp>
#define PI 3.1415926535897932
using std::vector;
using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::runtime_error;
using morph::HexGrid;
using morph::HdfData;
using morph::ReadCurves;
using morph::Tools;
using morph::BezCurve;
using morph::BezCurvePath;
using morph::BezCoord;
using namespace std;

class ksSolver
{
public:
    int scale;
    int n;
    double ds;
	//double overds;

  // list of objects visible to member functions
  vector<double> rad;
  vector<vector<int> > N; // hex neighbourhood 
  vector<double> psi;
  vector<double> NN, CC; //hold the field values for each he
  morph::HexGrid* Hgrid;
  // Class constructor
    ksSolver (int scale, string logpath, BezCurvePath bound) {
    ofstream afile (logpath + "/ksdebug.out" );
    this->scale = scale;
    double s = pow(2.0, scale-1);
	ds = 1.0/s;
    n = 0;
    Hgrid = new HexGrid(this->ds, 2.0, 0.0, morph::HexDomainShape::Boundary);
    n = Hgrid->num();
    afile << " max x " << Hgrid->getXmax(0.0) << " min x " << Hgrid->getXmin(0.0) << endl; 
    afile << "before filling H " << Hgrid->num() << endl;
    afile << "after creating HexGrid"<<endl;
    Hgrid->setBoundaryDRegion (bound);
    afile << "after setting boundary on  H " << Hgrid->num() << endl;
    n = Hgrid->num();
    afile << " max x " << Hgrid->getXmax(0.0) << " min x " << Hgrid->getXmin(0.0) << endl; 
    afile << "after  filling H " << " n = " << n <<endl;
    N.resize(n);
    this->psi.resize(n,0.0); // vector of azimuths
    this->rad.resize(n,0.0); // vector of radii
// making a neighbour array for convenience
   for (int idx = 0; idx < n; idx++) {
     N[idx].resize(6);
     N[idx][0] = Hgrid->d_ne[idx];
     N[idx][1] = Hgrid->d_nne[idx];
     N[idx][2] = Hgrid->d_nnw[idx];
     N[idx][3] = Hgrid->d_nw[idx];
     N[idx][4] = Hgrid->d_nsw[idx];
     N[idx][5] = Hgrid->d_nse[idx];
  }
    NN.resize(n);
    CC.resize(n);
	afile << "after alloc NN and CC" <<endl;
 } // end of ksSolver constructor 


    vector<double> getLaplacian(vector<double> Q, double dx) {
        double overds = 1./(1.5*29.0*29.0*dx*dx);
        vector<double> L(n,0.);
        for(auto h : this->Hgrid->hexen){
         int i = int(h.vi);
            L[i]=(Q[N[i][0]]+Q[N[i][1]]+Q[N[i][2]]+Q[N[i][3]]+Q[N[i][4]]+Q[N[i][5]]-6.*Q[i])*overds;
        }
        return L;
    }
        
		vector<double> chemoTaxis(vector<double> Q, vector<double> P, double dx) {
		vector<double> cT(n,0.);
          double overds = 1./(1.5*29.0*29.0*dx*dx);

        for (auto h : Hgrid->hexen) {
		  unsigned int i = h.vi;
        // finite volume method Lee et al. https://doi.org/10.1080/00207160.2013.864392
	      double dr0Q = (Q[N[i][0]]+Q[i])/2.;
	      double dg0Q = (Q[N[i][1]]+Q[i])/2.;
	      double db0Q = (Q[N[i][2]]+Q[i])/2.;
	      double dr1Q = (Q[N[i][3]]+Q[i])/2.;
	      double dg1Q = (Q[N[i][4]]+Q[i])/2.;
	      double db1Q = (Q[N[i][5]]+Q[i])/2.;
	  //double ncentre = Q[i];

          double dr0P = P[N[i][0]]-P[i];
          double dg0P = P[N[i][1]]-P[i];
          double db0P = P[N[i][2]]-P[i];
	      double dr1P = P[N[i][3]]-P[i];
          double dg1P = P[N[i][4]]-P[i];
          double db1P = P[N[i][5]]-P[i];


	  //finite volume for NdelC, h = s/2
	  cT[i] = (dr0Q*dr0P+dg0Q*dg0P+db0Q*db0P+dr1Q*dr1P+dg1Q*dg1P+db1Q*db1P)*overds;

	  //G[i] = (drN*drC+dgN*dgC+dbN*dbC)/(6.*ds*ds) + NN[i]*lapC[i];

        } //matches for on i
		return cT;
	} //end of function chemoTaxis

  //function to timestep coupled equations
    void step(double dt, double Dn, double Dchi, double Dc) {
      dt = dt * 2.5 / Dn;


       // Set up boundary conditions with ghost points
      //cout << " in time step before ghost points" << endl;
      for(auto h : Hgrid->hexen){
	   // cout << "top of ghost loop hex " << h.vi << " x " << h.x << " y " << h.y << endl;
        if(!h.onBoundary()){
          for(int j=0;j<6;j++){
		     int i = int(h.vi);
		 // cout << "in ghost loop j = " << j << " i= " << i << " Nbr " << N[h.vi][j] << endl;
             if(N[h.vi][j] == i) {
       	       NN[N[h.vi][j]] = NN[h.vi];
			//   cout << " NN " << NN[N[h.vi][j]] << " NN central " << NN[h.vi] << endl;
	           CC[N[h.vi][j]] = CC[h.vi];
	      }
	     }
	   }
      }



        double beta = 5.;
        double a = 1., b = 1., mu = 1;
		//cout << " before calls to Laplacian " << endl;
        vector<double> lapN = getLaplacian(NN,ds);
        vector<double> lapC = getLaplacian(CC,ds);
        vector<double> cTaxis = chemoTaxis(NN,CC,ds);

        // step N
        for (auto h : Hgrid->hexen) {
          NN[h.vi]+=dt*( a-b*NN[h.vi] + Dn*lapN[h.vi] - Dchi*cTaxis[h.vi]);
        }

        // step C
        double N2;
        for(auto h : Hgrid->hexen){
		    unsigned int i = h.vi;
            N2 = NN[i]*NN[i];
            CC[i]+=dt*( beta*N2/(1.+N2) - mu*CC[i] + Dc*lapC[i] );
        }
    }//end step

}; //end of class KSsolver
