/*
  ksSolver class
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
// list of objects visible to member functions
    int scale;
	double xspan;
    int n;
    double ds;
	double nnInitialOffset = 1.0;
	double ccInitialOffset = 2.5;
	double boundaryFalloffDist = 0.024;
    pair<double,double> seedPoint;
	BezCurvePath<float> bound;
	string logpath;
    vector<vector<int> > N; // hex neighbourhood
    vector<double> NN, CC; //hold the field values for each he
    morph::HexGrid* Hgrid;
// empty constructor
    ksSolver(){}
// Constructor with boundary passed in
    ksSolver (int scale, double xspan, string logpath, BezCurvePath<float> bound, pair<double,double> seedPoint) {
        this->scale = scale;
        this->xspan = xspan;
        this->logpath = logpath;
        this->bound = bound;
        this->seedPoint = seedPoint;
        ofstream afile (this->logpath + "/ksdebug.out",ios::app );
        double s = pow(2.0, this->scale-1);
        this->ds = 1.0/s;
        n = 0;
        Hgrid = new HexGrid(this->ds, this->xspan, 0.0, morph::HexDomainShape::Boundary);
        n = Hgrid->num();
        afile << " max x " << Hgrid->getXmax(0.0) << " min x " << Hgrid->getXmin(0.0) << endl;
        afile << "before filling H " << Hgrid->num() << endl;
        afile << "after creating HexGrid"<<endl;
        Hgrid->setBoundary(bound,false);
        afile << "after setting boundary on  H " << Hgrid->num() << endl;
        n = Hgrid->num();
        afile << " max x " << Hgrid->getXmax(0.0) << " min x " << Hgrid->getXmin(0.0) << endl;
        afile << "after  filling H in boundary constructor" << " n = " << n <<endl;
        N.resize(n);
        this->setHexType();
        Hgrid->computeDistanceToBoundary();
        for (auto &h : Hgrid->hexen) {
            afile  << "dist to bdry " << h.distToBoundary << " for " << h.vi << endl;
        }
        NN.resize(n);
        CC.resize(n);
        afile << "after alloc NN and CC" <<endl;
        pair<double, double> centroid = set_kS_polars(this->seedPoint);
        cout << " end of ksSolver constructor " << " x seedPoint " << seedPoint.first << " y seedPoint " << seedPoint.second << endl;
        cout << " end of ksSolver constructor " << " centroid x " << centroid.first << " centroid y " << centroid.second << endl;
    } // end of ksSolver constructor

// constructor with radius passed in for solving on radial boundaries
    ksSolver (int scale, double xspan, string logpath, float radius, pair<float, float> seedPoint) {
        this->scale = scale;
        this->xspan = xspan;
        this->logpath = logpath;
        this->bound = bound;
        this->seedPoint = seedPoint;
        ofstream afile (this->logpath + "/ksdebug.out",ios::app );
        double s = pow(2.0, this->scale-1);
        this->ds = 1.0/s;
        n = 0;
        Hgrid = new HexGrid(this->ds, this->xspan, 0.0, morph::HexDomainShape::Boundary);
        n = Hgrid->num();
        afile << " max x " << Hgrid->getXmax(0.0) << " min x " << Hgrid->getXmin(0.0) << endl;
        afile << "before filling H " << Hgrid->num() << endl;
        afile << "after creating HexGrid with radius "<< radius << endl;
        Hgrid->setCircularBoundary(radius, seedPoint, false);
        afile << "after setting boundary on  H " << Hgrid->num() << endl;
        n = Hgrid->num();
        afile << " max x " << Hgrid->getXmax(0.0) << " min x " << Hgrid->getXmin(0.0) << endl;
        afile << "after  filling H in circular constructor" << " n = " << n <<endl;
        N.resize(n);
         for (auto &h : this->Hgrid->hexen){
            this->N[h.vi].resize(6);
            if (!HAS_NE(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][0] = h.vi;
            }
            else {
                this->N[h.vi][0] = Hgrid->d_ne[h.vi];
            }

            if (!HAS_NNE(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][1] = h.vi;
            }
            else {
                this->N[h.vi][1] = Hgrid->d_nne[h.vi];
            }

            if (!HAS_NNW(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][2] = h.vi;
            }
            else {
                this->N[h.vi][2] = Hgrid->d_nnw[h.vi];
            }

            if (!HAS_NW(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][3] = h.vi;
            }
            else {
                this->N[h.vi][3] = Hgrid->d_nw[h.vi];
            }

            if (!HAS_NSW(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][4] = h.vi;
            }
            else {
                this->N[h.vi][4] = Hgrid->d_nsw[h.vi];
            }

            if (!HAS_NSE(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][5] = h.vi;
            }
            else {
                this->N[h.vi][5] = Hgrid->d_nse[h.vi];
            }
        } //end of loop over HexGri this->setHexType();
        /*
        Hgrid->computeDistanceToBoundary();
        for (auto &h : Hgrid->hexen) {
            afile  << "dist to bdry " << h.distToBoundary << " for " << h.vi << endl;
        }
        */
        NN.resize(n);
        CC.resize(n);
        afile << "after alloc NN and CC" <<endl;
        pair<double, double> centroid = set_kS_polars(this->seedPoint);
        cout << " end of ksSolver constructor " << " x seedPoint " << seedPoint.first << " y seedPoint " << seedPoint.second << endl;
        cout << " end of ksSolver constructor " << " centroid x " << centroid.first << " centroid y " << centroid.second << endl;
    } // end of ksSolver constructor


// this determines the type of hex
    void setHexType() {
        for (auto &h : this->Hgrid->hexen){
            this->N[h.vi].resize(6);
            if (!HAS_NE(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][0] = h.vi;
            }
            else {
                this->N[h.vi][0] = Hgrid->d_ne[h.vi];
            }

            if (!HAS_NNE(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][1] = h.vi;
            }
            else {
                this->N[h.vi][1] = Hgrid->d_nne[h.vi];
            }

            if (!HAS_NNW(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][2] = h.vi;
            }
            else {
                this->N[h.vi][2] = Hgrid->d_nnw[h.vi];
            }

            if (!HAS_NW(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][3] = h.vi;
            }
            else {
                this->N[h.vi][3] = Hgrid->d_nw[h.vi];
            }

            if (!HAS_NSW(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][4] = h.vi;
            }
            else {
                this->N[h.vi][4] = Hgrid->d_nsw[h.vi];
            }

            if (!HAS_NSE(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][5] = h.vi;
            }
            else {
                this->N[h.vi][5] = Hgrid->d_nse[h.vi];
            }
        } //end of loop over HexGrid
    } //end of method setHexType

// method to calculate the Laplacian
    vector<double> getLaplacian(vector<double> Q, double dx) {
        double overds = 1./(1.5*29.0*29.0*dx*dx);
        vector<double> L(n,0.);
        for(auto &h : this->Hgrid->hexen){
         int i = int(h.vi);
            L[i]=(Q[N[i][0]]+Q[N[i][1]]+Q[N[i][2]]+Q[N[i][3]]+Q[N[i][4]]+Q[N[i][5]]-6.*Q[i])*overds;
        }
        return L;
    }

		vector<double> chemoTaxis(vector<double> Q, vector<double> P, double dx) {
		vector<double> cT(n,0.);
          double overds = 1./(1.5*29.0*29.0*dx*dx);

        for (auto &h : Hgrid->hexen) {
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

  // function to compute the derivative
     void compute_dNNdt(vector<double>& inN, vector<double>& dNdt, double Dn, double Dchi) {
        vector<double> lapN(this->n,0);
		vector<double> cTaxis(this->n,0);
		//cout << "in compute_dNN just before laplacian" << endl;
		lapN = getLaplacian(inN,this->ds);
        cTaxis = chemoTaxis(inN,CC,this->ds);
        double a = 1., b = 1.;
		//cout << "in compute_dNN just before the loop" << endl;
        for (int h=0; h < this->n; h++) {
          dNdt[h] = a-b*inN[h] + Dn*lapN[h] - Dchi*cTaxis[h];
        }
	}

	void compute_dCCdt(vector<double>& inC, vector<double>&  dCdt, double Dc) {
        double beta = 5.;
        double mu = 1;
		//cout << " before calls to Laplacian " << endl;
        vector<double> lapC(this->n,0);
		lapC = getLaplacian(inC,this->ds);
        double N2;
        for(int h=0; h < this->n;h++){
            N2 = this->NN[h]*this->NN[h];
            dCdt[h] =  beta*N2/(1.+N2) - mu*inC[h] + Dc*lapC[h];
        }
    }

  //function to timestep coupled equations solely b.c. on the flux
    void step(double dt, double Dn, double Dchi, double Dc)
	{
        dt = dt * 2.5 / Dn;

        // cout  << "value of NN[5] start Runge " << this->NN[5] << endl;


        // Set up boundary conditions with ghost points
        //cout << " in time step before ghost points" << endl;
        for(auto &h : Hgrid->hexen)
		{
	   // cout << "top of ghost loop hex " << h.vi << " x " << h.x << " y " << h.y << endl;
            if(h.boundaryHex())
			{
                for(int j=0;j<6;j++)
				{
		            int i = int(h.vi);
		 // cout << "in ghost loop j = " << j << " i= " << i << " Nbr " << N[h.vi][j] << endl;
                    if(N[h.vi][j] == i)
					{
       	                NN[N[h.vi][j]] = NN[h.vi];
			//   cout << " NN " << NN[N[h.vi][j]] << " NN central " << NN[h.vi] << endl;
	                    CC[N[h.vi][j]] = CC[h.vi];
	                }
	            }
	         }
          }

        // 2. Do integration of NN
        {
            // Runge-Kutta integration for A. This time, I'm taking
            // ownership of this code and properly understanding it.

            // Ntst: "A at a test point". Ntst is a temporary estimate for A.
            vector<double> Ntst(this->n, 0.0);
            vector<double> dNdt(this->n, 0.0);
            vector<double> K1(this->n, 0.0);
            vector<double> K2(this->n, 0.0);
            vector<double> K3(this->n, 0.0);
            vector<double> K4(this->n, 0.0);

            /*
             * Stage 1
             */
			 //cout << "in ksSolver before compute_dNNdt" << endl;
            this->compute_dNNdt (this->NN, dNdt, Dn, Dchi);
			 //cout << "in ksSolver after compute_dNNdt" << endl;
            for (int h=0; h< this->n; ++h) {
                K1[h] = dNdt[h] * dt;
                Ntst[h] = this->NN[h] + K1[h] * 0.5 ;
            }

            /*
             * Stage 2
             */
            this->compute_dNNdt (Ntst, dNdt, Dn, Dchi);
            for (int h=0; h< this->n; ++h) {
                K2[h] = dNdt[h] * dt;
                Ntst[h] = this->NN[h] + K2[h] * 0.5;
            }

            /*
             * Stage 3
             */
            this->compute_dNNdt (Ntst, dNdt, Dn, Dchi);
            for (int h=0; h < this->n; ++h) {
                K3[h] = dNdt[h] * dt;
                Ntst[h] = this->NN[h] + K3[h];
            }

            /*
             * Stage 4
             */
            this->compute_dNNdt (Ntst, dNdt, Dn, Dchi);
            for (int h=0; h < this->n; ++h) {
                K4[h] = dNdt[h] * dt;
            }

            /*
             * Final sum together. This could be incorporated in the
             * for loop for Stage 4, but I've separated it out for
             * pedagogy.
             */
            for (int h=0;h<this->n;h++) {
                this->NN[h] += ((K1[h] + 2.0 * (K2[h] + K3[h]) + K4[h])/(double)6.0);
				//this->NN[i] = i * 1.0;
            }
        }

        // 3. Do integration of B
        {
            // Ctst: "B at a test point". Ctst is a temporary estimate for B.
            vector<double> Ctst(this->n, 0.0);
            vector<double> dCdt(this->n, 0.0);
            vector<double> K1(this->n, 0.0);
            vector<double> K2(this->n, 0.0);
            vector<double> K3(this->n, 0.0);
            vector<double> K4(this->n, 0.0);

            /*
             * Stage 1
             */
            this->compute_dCCdt (this->CC, dCdt, Dc);
            for (int h=0; h < this->n; ++h) {
                K1[h] = dCdt[h] * dt;
                Ctst[h] = this->CC[h] + K1[h] * 0.5 ;
            }

            /*
             * Stage 2
             */
		    this->compute_dCCdt (Ctst, dCdt, Dc);
            for (int h=0; h < this->n; ++h) {
                K2[h] = dCdt[h] * dt;
                Ctst[h] = this->CC[h] + K2[h] * 0.5;
            }

            /*
             * Stage 3
             */
            this->compute_dCCdt (Ctst, dCdt, Dc);
            for (int h=0; h < this->n; ++h) {
                K3[h] = dCdt[h] * dt;
                Ctst[h] = this->CC[h] + K3[h];
            }

            /*
             * Stage 4
             */
            this->compute_dCCdt (Ctst, dCdt, Dc);
            for (int h=0; h < this->n; ++h) {
                K4[h] = dCdt[h] * dt;
            }

            /*
             * Final sum together. This could be incorporated in the
             * for loop for Stage 4, but I've separated it out for
             * pedagogy.
             */
            for (int h=0; h < this->n; ++h) {
                this->CC[h] += ((K1[h] + 2.0 * (K2[h] + K3[h]) + K4[h])/(double)6.0);
            }
        }
        //cout  << "value of NN[5] end Runge " << this->NN[5] <<  " number of hexes " << this->n << endl;
    }//end step

  //function to timestep coupled equations option to set boundary to constant value
    void step(double dt, double Dn, double Dchi, double Dc, int steps, int numAdjust)
	{
        dt = dt * 2.5 / Dn;


      if ((steps%numAdjust == 0) && (steps/numAdjust != 0))
	  {
	     cout << "in numAdjust if step " << steps << endl;
	     for (auto &h : this->Hgrid->hexen)
	     {
		     //cout << "dist to bdry" << h.distToBoundary << endl;
	         if (h.distToBoundary > -0.5)
		     { // It's possible that distToBoundary is set to -1.0
                double bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary- this->boundaryFalloffDist)) );
				//cout << "bSig " << bSig << " for hex " << h.vi << endl;
                this->NN[h.vi] = (this->NN[h.vi] - this->nnInitialOffset) * bSig + this->nnInitialOffset;
                this->CC[h.vi] = (this->CC[h.vi] - this->ccInitialOffset) * bSig + this->ccInitialOffset;
		     } //end of if on boundary distance
	     }//end of loop over hexGrid
	  } //end of code applied to keep boundary conditions static

        // Set up boundary conditions with ghost points
        //cout << " in time step before ghost points" << endl;
        for(auto &h : Hgrid->hexen)
		{
	   // cout << "top of ghost loop hex " << h.vi << " x " << h.x << " y " << h.y << endl;
            if(h.boundaryHex())
			{
                for(int j=0;j<6;j++)
				{
		            int i = int(h.vi);
		 // cout << "in ghost loop j = " << j << " i= " << i << " Nbr " << N[h.vi][j] << endl;
                    if(N[h.vi][j] == i)
					{
       	                NN[N[h.vi][j]] = NN[h.vi];
			//   cout << " NN " << NN[N[h.vi][j]] << " NN central " << NN[h.vi] << endl;
	                    CC[N[h.vi][j]] = CC[h.vi];
	                }
	            }
	         }
          }
          void step(double dt, double Dn, double Dchi, double Dc);
    }//end step

    void reverse_y ()
	{
	  for (auto &h : this->Hgrid->hexen)
	    {
	      cout << " in reverse_y " << h.vi << endl;
	      //int index = h.vi;
	      cout << " in y reversing loop " << h.vi << endl;
	      double temp = double(this->Hgrid->d_y[h.vi]);
	      cout << " after getting y " << temp << endl;
		  if (temp != 0.0)
		  {
	          cout << " in y reversing loop " << endl;
	          this->Hgrid->d_y[h.vi] = -temp;
		  }
		}
	}

// function to give r and theta relative to region centre
    pair <double,double> set_kS_polars(pair<double,double> centre){
        pair <double, double> result;
		result.first = 0.0;
		result.second = 0.0;
        double xav=0;
        double yav = 0;
        int hexcount = 0;
	    this->reverse_y();
	    cout <<"in set polars ksSolver xcentre" << centre.first << " y_centre  " << centre.second <<endl;
        for (auto &h : this->Hgrid->hexen) {
            hexcount++;
            xav += this->Hgrid->d_x[h.vi];
            yav += this->Hgrid->d_y[h.vi];
        }
		//cout << "after set centres " << endl;
		if (hexcount != 0) {
        xav = xav / (hexcount*1.0);
        yav = yav / (hexcount*1.0);
		}
		else {
		  //cout << " in set_polars no hexes in region "<<endl;
		  }
//go over the region and put the hexes into bins then average
        for (auto&  h : this->Hgrid->hexen) {
            int index = h.vi;
			double angle = 0;
			double dx = this->Hgrid->d_x[index];
			double dy = this->Hgrid->d_y[index];
			//cout <<"in set polars ksSolver index " << index << " i " << h.vi <<endl;
			cout << "d_x " << dx << " dy " << dy <<endl;
            h.r = sqrt((dx - centre.first)*(dx - centre.first)
			+ (dy - centre.second)*(dy - centre.second));
            if (dy >= centre.second) {
              angle =   atan2((dy - centre.second), (dx - centre.first));
			  h.phi = angle;
			  cout<< " setPhi if 1 " << h.phi<<  " index " << h.vi << endl;
			  }
            else {
              angle =  2*PI + atan2((dy - centre.second), (dx - centre.first));
			  h.phi = angle;
			  cout<< " setPhi if 2 " << h.phi<<  " index " << h.vi << endl;
			  }
        }
/*
        result.first = centre.first - xcentre; //diff between seed point and barycentre
        result.second = centre.second - ycentre;
*/
        result.first = xav - centre.first ; // barycentre
        result.second = yav - centre.second ;
		cout << "centre x "<< centre.first << " centre y " << centre.second << " result.first " << result.first << " result.second " << result.second <<endl;
        return result;
    } //end of function set_polars

}; //end of class KSsolver
