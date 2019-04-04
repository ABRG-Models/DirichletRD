#include "morph/world.h"
#include "morph/sockserve.h"
#include "morph/display.h"
#include "morph/tools.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <math.h>
#include <boost/math/special_functions/bessel.hpp>
#define PI 3.14159265

using namespace std;

class Erm2009
{
public:
    int scale;
    int offset;
    int n;
    double ds;
    struct extremum {
    int radialIndex;
    double radialValue;
     };

    vector<vector<double> > X; //
    vector<vector<double> > H; // hex-grid info
    vector<vector<double> > N; // hex neighbourhood
    vector<int> C;
    vector<extremum> turnVal;

    vector<double> NN, CC;
       double maxRadius = sqrt(0.75) * 29.0;
  //class constructor
  Erm2009 (int scale, int offset, double z) {



        this->scale = scale;
        this->offset = offset;

        H.resize(6);
        double s = pow(2.0, scale-1); //make S bigger to incorporate the circle.
	ds = 29.0/s;
        n = 0;
	double x = 0;
	double y = 0;
	double radius = 0;
	double theta = 0;
	for (int r=-s; r<0; r++) {
	  for (int b=-s; b<=s; b++) {
	 	      x = (r-b)*sqrt(3.)/(2.*s);
	 	      y = (1.*(r+b))/(2.*s);
		      radius = 29.0*sqrt(x*x+y*y);
		       if (radius == 0)
			theta = 0;
		      else
		        if (y >= 0)
			  theta = atan2(y,x);
		        else
			  theta = 2*PI + atan2(y,x);


		             if (x*x + y*y <= 0.75) { // circular
                                 //if (x*x/0.1 + y*y/0.8 <= 1) { // ellipse 1
                                 //if (x*x/0.8 + y*y/0.1 <= 1) { // ellipse 2
                                 //if (x*x/0.8 + y*y/0.01 <= 0.9) { // ellipse 3
                                 //if (x*x/0.01 + y*y/0.8 <= 0.9) { // ellipse 4
                                 //if (x*x < 0.2 && y*y < 0.45) { // square 1
                                 //if (x*x < 0.2 && y*y < 0.2) { // square 2
                                 //if (x*x/0.05 - y*y/0.05 <= 0.7) { // hyperbolic 1
                              H[0].push_back(x);                      //X
                              H[1].push_back(y);                      //Y
                              H[2].push_back(r);                      //r
                              H[3].push_back(b);                      //b
                              H[4].push_back(radius);                 //r
                              H[5].push_back(theta);                  //theta
                              n++;
                          } //if on shape determination
	  } //if on b
	} //if on r


	for (int r=0; r<=s; r++) {
            for (int b=-s; b<=s; b++) {
		      x = (r-b)*sqrt(3.)/(2.*s);
	 	      y = (1.*(r+b))/(2.*s);
		      radius = 29.0*sqrt(x*x+y*y);
		      if (radius == 0)
			theta = 0;
		      else
		        if (y >= 0)
			  theta = atan2(y,x);
		        else
			  theta = 2*PI + atan2(y,x);
                      if (x*x + y*y <= 0.75) { // circular
                                  //if (x*x/0.1 + y*y/0.8 <= 1) { // ellipse 1
                                  //if (x*x/0.8 + y*y/0.1 <= 1) { // ellipse 2
                                  //if (x*x/0.8 + y*y/0.01 <= 0.9) { // ellipse 3
                                  //if (x*x/0.01 + y*y/0.8 <= 0.9) { // ellipse 4
                                  //if (x*x < 0.4 && y*y < 0.4) { // square 1
                          // if (x*x < 0.2 && y*y < 0.45) { // square 2
                                  //if (x*x/0.05 - y*y/0.05 <= 0.7) { // hyperbolic 1
                               H[0].push_back(x);                      //X
                               H[1].push_back(y);                      //Y
                               H[2].push_back(r);                      //r
                               H[3].push_back(b);                      //b
                               H[4].push_back(radius);                      //r
                               H[5].push_back(theta);                  //theta
                               n++;
                           } //if on shape determination
	   	} //if on b
	    } //if on r





        // get neighbours
        N.resize(n);
        C.resize(n,0);
        for(int i=0;i<n;i++){
            N[i].resize(6,i); // CONNECT ALL TO 'boundary' UNIT AT N+1
            for(int j=0;j<n;j++){
                int dr = H[2][j]-H[2][i];
		// int dg = H[3][j]-H[3][i];
                int db = H[3][j]-H[3][i];
		// if(max(max(abs(dr),abs(dg)),abs(db))==1){
		   if(max(abs(dr),abs(db))==1){
                    //anticlockwise from east
                      if(db==0&&dr==+1){N[i][0]=j;C[i]++;}
                      if(db==+1&&dr== +1){N[i][1]=j;C[i]++;}
                      if(db== +1&&dr==0){N[i][2]=j;C[i]++;}
                      if(db==0&&dr==-1){N[i][3]=j;C[i]++;}
                      if(db==-1&&dr==-1){N[i][4]=j;C[i]++;}
                      if(db== -1&&dr==0){N[i][5]=j;C[i]++;}

                     //previous definition
		     // if(dg==0&&dr==+1){N[i][0]=j;C[i]++;}
                     // if(dg==1&&dr== 0){N[i][1]=j;C[i]++;}
                     // if(dg==0&&dr==-1){N[i][2]=j;C[i]++;}
                     // if(dg==-1&&dr==-1){N[i][3]=j;C[i]++;}
                     // if(dg==-1&&dr== 0){N[i][4]=j;C[i]++;}
                     // if(dg==-1&&dr==+1){N[i][5]=j;C[i]++;}
		    }
            }
        }
        X.resize(n);
        for(int i=0;i<n;i++){
            X[i].resize(3,0.);
            X[i][0] = H[0][i];
            X[i][1] = H[1][i];
        }
        NN.resize(n);
        CC.resize(n);
  } //end of constructor

    vector<double> getLaplacian(vector<double> Q, double dx) {
      double overdxSquare = 1./(1.5*dx*dx);
        vector<double> L(n,0.);
        for(int i=0;i<n;i++){
            L[i]=(Q[N[i][0]]+Q[N[i][1]]+Q[N[i][2]]+Q[N[i][3]]+Q[N[i][4]]+Q[N[i][5]]-6.*Q[i])*overdxSquare;
        }
        return L;
    }

    void step(double dt, double Dn, double Dc) {


       // Set up boundary conditions with ghost points

        for(int i=0;i<n;i++){
            if(C[i]<6){
                for(int j=0;j<6;j++){
		  if(N[i][j] == i) {
		    // int temp = (j+3) % 6;
		    //   NN[N[i][j]] = NN[N[i][temp]];
		    //   CC[N[i][j]] = CC[N[i][temp]];
		    // I think this is the correct ghost point to enforce no-flux
		     NN[N[i][j]] = NN[i];
		     CC[N[i][j]] = CC[i];
		  }
                }
	    }
	}



        double beta = 5.;
        double a = 1., b = 1., mu = 1., chi = Dn;
        vector<double> lapN = getLaplacian(NN,ds);
        vector<double> lapC = getLaplacian(CC,ds);


        vector<double> G(n,0.);

        for (int i=0;i<n;i++) {
        // finite volume method Lee et al. https://doi.org/10.1080/00207160.2013.864392
	  double dr0N = (NN[N[i][0]]+NN[i])/2.;
	  double dg0N = (NN[N[i][1]]+NN[i])/2.;
	  double db0N = (NN[N[i][2]]+NN[i])/2.;
	  double dr1N = (NN[N[i][3]]+NN[i])/2.;
	  double dg1N = (NN[N[i][4]]+NN[i])/2.;
	  double db1N = (NN[N[i][5]]+NN[i])/2.;
	  //double ncentre = NN[i];

          double dr0C = CC[N[i][0]]-CC[i];
          double dg0C = CC[N[i][1]]-CC[i];
          double db0C = CC[N[i][2]]-CC[i];
	  double dr1C = CC[N[i][3]]-CC[i];
          double dg1C = CC[N[i][4]]-CC[i];
          double db1C = CC[N[i][5]]-CC[i];


	  //finite volume for NdelC, h = s/2
	  G[i] = (dr0N*dr0C+dg0N*dg0C+db0N*db0C+dr1N*dr1C+dg1N*dg1C+db1N*db1C)/(1.5*ds*ds);

	  //G[i] = (drN*drC+dgN*dgC+dbN*dbC)/(6.*ds*ds) + NN[i]*lapC[i];

        } //matches for on i

        // step N
        for(int i=0;i<n;i++){
            NN[i]+=dt*( a-b*NN[i] + Dn*lapN[i] - chi*G[i]);
        }

        // step C
        double N2;
        for(int i=0;i<n;i++){
            N2 = NN[i]*NN[i];
            CC[i]+=dt*( beta*N2/(1.+N2) - mu*CC[i] + Dc*lapC[i] );
        }
    }//end step

//function to smooth a vector by moving average
    vector <double> smooth_vector(vector<double> invector, int window) {
        vector<double> outvector;
        int size = invector.size();
        outvector.resize(size);
        for (int i=1; i<size+1; i++) {
            outvector[i%size] = (invector[(i -1)%size] + invector[i%size] + invector[(i + 1)%size])/3.0;
        }
        return outvector;
    }

// to find the maximum value of a vector
    double maxVal( vector<double> invector) {
            double result = -1.0e7;
            for (unsigned int i = 0; i < invector.size(); i++) {
                    if (invector[i] > result)
                            result = invector[i];
            }
            return result;
    }


    // to find the minimum value of a vector
    double minVal( vector<double> invector) {
            double result = 1.0e7;
            for (unsigned int i = 0; i < invector.size(); i++) {
                    if (invector[i] < result)
                            result = invector[i];
            }
            return result;
    }


    // transform vector so its mean is zero
    vector<double> meanzero_vector(vector<double> invector) {
        ofstream meanzero ("meanzero.txt",ios::app);
        vector <double> result;
        int size = invector.size();
        meanzero << "size " << size << endl;
        double sum = 0;
        double absSum = 0;
        for (int i=0; i <size; i++) {
            sum += invector[i];
            absSum += fabs(invector[i]);
	    meanzero << " i " << invector[i] << " sum " << sum << endl;
        }
        sum = sum / (1.0 * size );
        absSum = absSum / (1.0 * size);
        meanzero << " mean  " << sum << endl;
        meanzero << " absolute mean  " << absSum << endl;
        for (int i=0; i < size; i++) {
            result.push_back(invector[i] - sum);
            //meanzero << " i " << result[i] << endl;
        }
        return result;
    }

  //function find_max to find turning points both values and indices.

  int find_maxR(vector<double> ray) {

    int iend = ray.size();
    //int iend = 0;
    turnVal.resize(1000);
    cout <<"ray size ="<<iend<<iend<<flush;
    double old_slope = 0;
    double new_slope = 0;
    int count = 0;
    old_slope = ray[1] - ray[0];
    for (int i =2; i<iend;i++){
      new_slope = ray[i]-ray[i-1];
      if (new_slope*old_slope < 0.) {
	turnVal[count].radialIndex = i;
	turnVal[count].radialValue = ray[i]; //should really interpolate
        count++;

      old_slope = new_slope;
      }
    }
    return count;
  }

  int find_maxTh(vector<double> ray) {
    ofstream cfile ( "logsWhole/angleTurn.txt" );
    int iend = ray.size();
    turnVal.resize(1000);
    cout <<"ray size ="<<iend<<iend<<flush;
    double old_slope = 0;
    double new_slope = 0;
    int count = 0;
    old_slope = (ray[1] - ray[0]) / fabs(ray[1] -ray[0]);
    for (int i =2; i<=iend+1;i++){

      new_slope = (ray[i%iend]-ray[(i-1)%iend]) / fabs((ray[i%iend]-ray[(i-1)%iend]));
      if (new_slope*old_slope < 0.000001) {
	turnVal[count].radialIndex = i;
	turnVal[count].radialValue = ray[i%iend]; //should really interpolate
 cfile << " index "<< turnVal[count].radialIndex << " value " << turnVal[count].radialValue <<endl;
       cfile << "new slope " << new_slope << " old slope" << old_slope<<endl;
       count++;
      }
	old_slope = new_slope;
    }

    return count;
  }
 
  // functions to find zeros
// find the zeros in a ray angular analogue version
    int find_zeroAngle(vector<double> ray, int window) {
    int size = ray.size();
    ofstream zerofile ("zero.txt",ios::app);
    vector<double> smoothRay;
    //smoothRay = this->smooth_vector(ray, window);
    smoothRay = ray;
    zerofile <<"size = " << size <<endl;
    turnVal.resize(1000);
    double old_val = 0;
    double new_val = 0;
    int count = 0;
    old_val = smoothRay[0];
    for (int i =1; i<size+1;i++){
      new_val = smoothRay[i%size];
      zerofile << " " << i%size << " " << old_val << " "<<new_val <<endl;
      if (new_val*old_val < 0.) {
        turnVal[count].radialIndex = i;
        turnVal[count].radialValue = smoothRay[i]; //should really interpolate
	zerofile << "in zeroAngle" << endl;
        zerofile << "turn index " << turnVal[count].radialIndex << " turn value " << turnVal[count].radialValue << endl;
        count++;
      }

    old_val = new_val;
    }
    return count;
  }


      // find the zeros in a ray angular digital version
    int find_zeroDAngle(vector<int> ray) {
    int size = ray.size();
    ofstream zerofile ("zero.txt",ios::app);
    int old_val = ray[0];
    int new_val;
    int count = 0;
    for (int i = 1 ; i<size+1;i++){
        if (ray[i%size] != 0) {
           new_val = ray[i%size];
           if (old_val*new_val == -1) {
                   count++;
           }
        old_val = new_val;
        }
	zerofile << "in zeroDangle" << endl;
      zerofile << " radius " << i%size << " " << old_val << " "<<new_val <<endl;
      }
    return count;
  }

     // find the zeros in a ray angular
    int find_zeroDRadius(vector<int> ray) {
    int size = ray.size();
    ofstream zerofile ("zero.txt",ios::app);
    int old_val = ray[0];
    int new_val;
    int count = 0;
    for (int i = 1 ; i<size;i++){
        if (ray[i%size] != 0) {
           new_val = ray[i%size];
           if (old_val*new_val == -1) {
                   count++;
           }
        old_val = new_val;
        }
	zerofile<<"in zeroDRadius "<<endl;
      zerofile << " angle " << i%size << " " << old_val << " "<<new_val <<endl;
      }
    return count;
  }

   // find the zeros in a ray angular
    int find_zeroRadius(vector<double> ray, int window) {
    int size = ray.size();
    ofstream zerofile ("zero.txt",ios::app);
    vector<double> smoothRay;
    smoothRay = ray;
    zerofile <<"size = " << size <<endl;
    turnVal.resize(1000);
    double old_val = 0;
    double new_val = 0;
    int count = 0;
    old_val = smoothRay[0];
    for (int i =1; i<size;i++){
      new_val = smoothRay[i%size];
      zerofile << " " << i%size << " " << old_val << " "<<new_val <<endl;
      if (new_val*old_val < 0.) {
        turnVal[count].radialIndex = i;
        turnVal[count].radialValue = smoothRay[i]; //should really interpolate
	zerofile<<"in zeroRadius"<<endl;
        zerofile << "turn index " << turnVal[count].radialIndex << " turn value " << turnVal[count].radialValue << endl;
        count++;
      }
      old_val = new_val;
    }
    return count;
    }


  //sectorize over radius
    //adapted for digital output
    vector <int> sectorize_reg_Dradius (int numSectors, int beginAngle, int endAngle) {
        ofstream dfile ( "logs/sectorRadius.txt",ios::app);
    vector <int>  radiusNN;
    radiusNN.resize(numSectors,0);
    vector <int> radiusCount;
    vector <double> radiusHold;
    vector <double> normalNN;
    radiusCount.resize(numSectors,0);
    radiusHold.resize(numSectors,0.0);
    double startRadius, finishRadius, radiusInc; //sector radii
    dfile << " maxRadius used " << this->maxRadius << endl;
    radiusInc = this->maxRadius /(1.0*numSectors);
    double startAngle, finishAngle, angleInc; //sector angles
    angleInc = 2*PI/(1.*numSectors);
    startAngle = beginAngle*angleInc;
    finishAngle = endAngle*angleInc;
    int size = this->n;
    dfile << "size = " << size <<endl;
    double sum = 0.0;
    for (int i=0;i<size;i++){
          normalNN.push_back(this->NN[i]);
	  sum += this->NN[i];
      }
    // to normalise the NN field
    normalNN = meanzero_vector(normalNN);
    dfile << " sum " << sum << endl;
    double epsilon = 0.0001*(maxVal(normalNN) - minVal(normalNN));
    dfile << "after the call to maxVal and minVal" <<endl;
      //for (int i=0;i<size;i++) {
        //   dfile << "in the i loop" << endl;
          // dfile << " i " << i << " normalNN[i] " << normalNN[i] << endl;
    //  }
    for (int k=0;k<numSectors;k++) {
      startRadius = (k*radiusInc);
      finishRadius = (k+1)*radiusInc;
      for (int i=0; i< size;i++) {
        if (this->H[5][i]>=startAngle && this->H[5][i]<finishAngle) {
          if (this->H[4][i]>=startRadius && this->H[4][i]<finishRadius) {
              radiusCount[k]++;
              radiusHold[k] += normalNN[i];
	  } //end of if on radius
	} //end of if on angleSector
      } //end of loop over i
     } //end of loop over k

     // radiusHold = meanzero_vector(radiusHold);

     for (int k=0;k<numSectors;k++) {
        startRadius = (k*radiusInc);
        finishRadius = (k+1)*radiusInc;
	     if (radiusCount[k] == 0) {
		     radiusNN[k] = 2;
		     continue;
        }
    //    radiusHold[k]  = radiusHold[k]  / (1.*radiusCount[k]);
        if (radiusHold[k] > epsilon)
                radiusNN[k] = 1;
        else if (radiusHold[k] < -epsilon)
                radiusNN[k] = -1;
        else
                radiusNN[k] = 0;

      dfile << "startRadius "<<startRadius << " radiusNN " << radiusNN[k] << " radiusHold  " << radiusHold[k] << endl;
    }//end loop on k
     return radiusNN;

  } //end of function sectorize_radius



     vector <double> sectorize_reg_radius ( int numSectors, int beginAngle, int endAngle) {
        ofstream dfile ( "logsWhole/sectorRadius.txt",ios::app);
    vector <double>  radiusNN;
    radiusNN.resize(numSectors,0);
    vector <int> radiusCount;
    vector <double> normalNN;
    radiusCount.resize(numSectors,0);
    double startRadius, finishRadius, radiusInc; //sector radii
    radiusInc = this->maxRadius /(1.0*numSectors);
    double startAngle, finishAngle, angleInc; //sector angles
    angleInc = 2*PI/(1.*numSectors);
    startAngle = beginAngle*angleInc;
    finishAngle = endAngle*angleInc;
    int size = this->n;
    // to normalise the NN field
     for (int i=0;i<size;i++){
          normalNN.push_back(this->NN[i]);
      }
      normalNN = meanzero_vector(normalNN);
      //for (int i=0;i<size;i++)
          // dfile << " i " << i << " normalNN[i] " << normalNN[i] << endl;

    for (int k=0;k<numSectors;k++) {
      startRadius = (k*radiusInc);
      finishRadius = (k+1)*radiusInc;

      for (int i=0; i< size;i++) {
        if (this->H[5][i]>=startAngle && this->H[5][i]<finishAngle) {
          if (this->H[4][i]>=startRadius && this->H[4][i]<finishRadius) {
              radiusCount[k]++;
              radiusNN[k] += normalNN[i];
          } //end of if on radius
        } //end of if on angleSector
      } //end of loop over i
if (radiusCount[k] != 0)
        //radiusNN[k] = radiusNN[k] / (1.*radiusCount[k]);
        radiusNN[k] = radiusNN[k]; 
      else
        radiusNN[k] = 0.0;
      dfile << "startRadius "<<startRadius<<"  finishRadius "<<finishRadius<< " radiusNN " << radiusNN[k] << endl;
    }//end loop on k
 return radiusNN;

  } //end of function sectorize_radius
 //function to count the hexes in sectors of a region via angular sectors
 //analogue version
    vector <double> sectorize_reg_angle (int numSectors, int beginRadius, int endRadius) {
     ofstream cfile ("logs/sectorAngle.txt",ios::app);
    //std::pair<double,double> diff; //difference between seed point and CoG of region
    // diff = this->set_polars(regNum);
    vector <double> angleNN; //average value of CC in each sector
    vector <double> normalNN;
    vector <double>  angleVal;
    vector <int> count; //number of hexes in each sector
    angleNN.resize(numSectors,0);
    count.resize(numSectors,0);
    double startAngle, endAngle, angleInc; //sector angles
    double startRadius, finishRadius,radiusInc;
    // double minRadius = min_radius(regNum);
    radiusInc = this->maxRadius/ (1.0*numSectors);
    startRadius = beginRadius*radiusInc;
    finishRadius = endRadius*radiusInc;
    angleInc = 2*PI/(1.*numSectors);
    // to normalise the NN field
    int size = this->n;
     for (int i=0;i<size;i++){
          normalNN.push_back(this->NN[i]);
      }
      normalNN = meanzero_vector(normalNN);
      //for (int i=0;i<size;i++)
        //   cfile << " i " << i << " normalNN[i] " << normalNN[i] << endl;

    for (int k=0;k<numSectors;k++) {
      //double angle;
      startAngle = k*angleInc;
      endAngle = (k+1)*angleInc;
                     if ((k+1) == numSectors)
         endAngle = 2*PI;


      for (int i=0; i< size;i++) {
          if (this->H[4][i]>=startRadius && this->H[4][i]<finishRadius) {
              if (this->H[5][i]>=startAngle && this->H[5][i]<endAngle) {
                  angleVal.push_back(H[5][i]);
                  count[k]++;
          //angle[k] += this->[regionIndex[regNum][i]];
                  angleNN[k] += normalNN[i];
              }
        //cfile << setw(5) << angleVal[i]  <<"  ";
          }
      }
    }

    angleNN = meanzero_vector(angleNN);
    for (int k=0;k<numSectors;k++) {
      startAngle = k*angleInc;
      endAngle = (k+1)*angleInc;
      cfile << endl;
      // cfile << "diff x " << diff.first << "diff y " << diff.second << endl;
      if (count[k] != 0)
        angleNN[k] = angleNN[k] / (1.*count[k]);
      else
        angleNN[k] = -999.999;
      cfile << "  start angle "  << startAngle << "  end angle "<< endAngle << " NN field " << angleNN[k] << endl;
    }  //end loop on k
     return angleNN;

  } //end of function sectorize_region


 //function to count the hexes in sectors of a region via angular sectors
 //Produces the angular degree, digitized version
    vector <int> sectorize_reg_Dangle (int numSectors, int beginRadius, int endRadius) {
        ofstream cfile ("logsWhole/sectorAngle.txt",ios::app);
    //std::pair<double,double> diff; //difference between seed point and CoG of region
    // diff = this->set_polars(regNum);
    vector <int> angleNN; //digitized value of NN in each sector
    vector <double> angleHold;
    vector <double> normalNN;
    vector <double>  angleVal;
    vector <int> angleCount; //number of hexes in each sector
    angleNN.resize(numSectors,0);
    angleHold.resize(numSectors,0.0);
    angleCount.resize(numSectors);
    double startAngle, endAngle, angleInc; //sector angles
    double startRadius, finishRadius,radiusInc;
    // double minRadius = min_radius(regNum);
    radiusInc = this->maxRadius/ (1.0*numSectors);
    startRadius = beginRadius*radiusInc;
    finishRadius = endRadius*radiusInc;
    angleInc = 2*PI/(1.*numSectors);
    // to normalise the NN field
    int size = this->n;
     for (int i=0;i<size;i++){
          normalNN.push_back(this->NN[i]);
      }
      normalNN = meanzero_vector(normalNN);
      double epsilon = 0.0001*(maxVal(normalNN) - minVal(normalNN));
     // for (int i=0;i<size;i++)
       //    cfile << " i " << i << " normalNN[i] " << normalNN[i] << endl;

    for (int k=0;k<numSectors;k++) {
      startAngle = k*angleInc;
      endAngle = (k+1)*angleInc;
      if ((k+1) == numSectors)
         endAngle = 2*PI;

      for (int i=0; i< size;i++) {
          if (this->H[4][i]>=startRadius && this->H[4][i]<finishRadius) {
              if (this->H[5][i]>=startAngle && this->H[5][i]<endAngle) {
                  angleVal.push_back(H[5][i]);
                  angleCount[k]++;
          //angle[k] += this->[regionIndex[regNum][i]];
                  angleHold[k] += normalNN[i];
              }
        //cfile << setw(5) << angleVal[i]  <<"  ";
          }
      }
    }
    angleHold = meanzero_vector(angleHold);
      cfile << endl;
      // cfile << "diff x " << diff.first << "diff y " << diff.second << endl;

      for (int k=0;k<numSectors;k++){
      startAngle = k*angleInc;
      endAngle = (k+1)*angleInc;
        if (angleCount[k] == 0) {
           angleNN[k] = 2;
           continue;
        }
        //angleHold[k] = angleHold[k] / (1.*angleCount[k]);
        if (angleHold[k] > epsilon)
                angleNN[k] = 1;
        else if (angleHold[k] < epsilon)
                angleNN[k] = -1;
        else
                angleNN[k] = 0;

       cfile << "  start angle "  << startAngle << "  end angle "<< endAngle << " NN field " << angleNN[k] << endl;
    } //end loop on k
     return angleNN;

  } //end of function sectorize_region


  //function bessel_ray for computing bessel functions along a ray
  // takes a vector of doubles representing the radius
  // returns the Bessel function values
  // vector<double> bessel_ray (int v, vector<double> ray) {
  //   vector <double> result(n,0.);
  //   result = cyl_bessel_j(v, ray);
  //   return result;
  // }


}; // Erm2009


int main (int argc, char **argv)
{
    srand(atoi(argv[3]));

    double dt = .0001;            // integration timestep // SHOULD BE 0.001

    // INITIALIZATION
    morph::World W(argv[1],argv[2],atoi(argv[3]),atoi(argv[4]),dt);
    //W.waitForConnection();

    // DISPLAYS
    vector<morph::Gdisplay> displays;
    vector<double>fix(3,0.);
    vector<double>eye(3,0.);
    eye[2] = -0.4;
    vector<double>rot(3,0.);
    displays.push_back (morph::Gdisplay (600, "NN field", 0., 0., 0.));
    displays[0].resetDisplay(fix,eye,rot);
    displays[0].redrawDisplay();
    displays.push_back (morph::Gdisplay (600, "CC field", 0., 0., 0.));
    displays[1].resetDisplay(fix,eye,rot);
    displays[1].redrawDisplay();


    Erm2009 M(7,0,0.);
   double  value;
   double factor = 7.0 / 29.0;
  vector <double> ray;
   for (int i=0;i<M.n;i++) {
     value  = boost::math::cyl_bessel_j(0,M.H[4][i]*factor);
       ray.push_back(value);
  }




   ofstream bfile ("logsWhole/sectorRadius.out");
   ofstream afile ("logsWhole/sectorAngle.out");


   for (int i=0;i<M.n;i++) {
  //    M.NN[i]=(morph::Tools::randDouble())*0.01 +1.;
    //  M.CC[i]=(morph::Tools::randDouble())*0.01+2.5;
     // M.CC[i]=sin(i*3.142)*0.1+2.5;
     //	 M.CC[i]=sin(i*3.142)*0.1+2.5;
     }

      for (int i=0;i<M.n;i++) {
     //   // M.NN[i]= 1. + sin(i*3.142)*sin(i*3.142)*0.001;
        M.NN[i]= 1. + (morph::Tools::randDouble())*0.1;
        M.CC[i]= 1. + (morph::Tools::randDouble())*0.1;
//        M.CC[i]=ray[i]*0.05+2.5;
       }

    unsigned int frameN = 0;
    unsigned int frameM = 0;
    //initial values for Dn Dc
    double Dn = 100.;
    double Dc = Dn*0.3;
    double TIME=0.;
    vector<double*> f;
    f.push_back(&TIME);







    vector<vector<double*> > S;
    {
        vector<double*> s;
        for(unsigned int i=0;i<M.CC.size();i++){
            s.push_back(&M.CC[i]);
        } S.push_back(s);
    }
    {
        vector<double*> s;
        for(unsigned int i=0;i<M.NN.size();i++){
            s.push_back(&M.NN[i]);
        } S.push_back(s);
    }

    bool doing = true;
    while (doing) {

        std::stringstream TIMEss;
        TIMEss<<setw(10)<<setfill('0')<<TIME;
        const char* TIMEcs = TIMEss.str().c_str();

        std::stringstream out;
        out.clear();
        out.setf(ios::fixed,ios::floatfield);

        // DEFINE OUTPUT MESSAGE
        out<<TIME<<",";

        vector <string> command;
        string messageI=W.master.exchange(out.str().c_str());
        stringstream ss(messageI);
        while (ss.good()){
            string substr;
            getline(ss,substr,',');
            command.push_back(substr);
        } ss.clear();

        // Interpret commands:
        switch (stoi(command[0])) {

        case 0: // *** QUIT ***
        {
            W.logfile<<W.processName<<"@"<<TIMEcs<<": 0=QUIT"<<endl<<flush;

            W.logfile.close();
            for(unsigned int i=0;i<displays.size();i++){
                displays[i].closeDisplay();
            }
            W.master.closeSocket();
            for(unsigned int i=0;i<W.ports.size();i++){
                W.ports[i].closeSocket();
            }
            doing=false;
            break;
        }

        case 1: // *** STEP ***
        {
	  //W.logfile<<W.processName<<"@"<<TIMEcs<<": 1=STEP"<<endl<<flush;

            // Exchange comms
            vector< vector <string> > commands(W.ports.size());
            for(unsigned int i=0;i<W.ports.size();i++){
                string messageI=W.ports[i].exchange(out.str().c_str());
                stringstream ss(messageI);
                while (ss.good()){
                    string substr;
                    getline(ss,substr,',');
                    commands[i].push_back(substr);
                } ss.clear();
            }

            // DO STUFF HERE

            M.step(dt, Dn, Dc);


            // END STEP
            TIME++;
            break;
        }

        case 2: // *** PLOT ***
        {
	  //W.logfile<<W.processName<<"@"<<TIMEcs<<": 2=PLOT"<<endl<<flush;
            displays[0].resetDisplay(fix,eye,rot);
	    displays[1].resetDisplay(fix,eye,rot);

            vector<double> plt = M.NN;
            double maxV = -1e7;
            double minV = +1e7;
            for (int i=0;i<M.n;i++) {
                if (M.C[i]==6) {
                    if (plt[i]>maxV) { maxV = plt[i]; }
                    if (plt[i]<minV) { minV = plt[i]; }
                }
            }
            double scaleV = 1./(maxV-minV);
            vector<double> P(M.n,0.);
            for (int i=0;i<M.n;i++) {
                P[i] = fmin (fmax (((plt[i])-minV)*scaleV,0.),1.);
                // M.X[i][2] = P[i];
            }

            for (int i=0;i<M.n;i++) {
                vector <double> cl = morph::Tools::getJetColor(P[i]);
                displays[0].drawTriFill(M.X[i],M.X[M.N[i][0]],M.X[M.N[i][1]],cl);
                displays[0].drawTriFill(M.X[i],M.X[M.N[i][3]],M.X[M.N[i][4]],cl);
            }
            displays[0].redrawDisplay();

	    vector<double> plu = M.CC;
            double maxV1 = -1e7;
            double minV1 = +1e7;
            for (int i=0;i<M.n;i++) {
                if (M.C[i]==6) {
                    if (plu[i]>maxV1) { maxV1 = plu[i]; }
                    if (plu[i]<minV1) { minV1 = plu[i]; }
                }
            }
            double scaleV1= 1./(maxV1-minV1);
            vector<double> P1(M.n,0.);
            for (int i=0;i<M.n;i++) {
                P1[i] = fmin (fmax (((plu[i])-minV1)*scaleV1,0.),1.);
                // M.X[i][2] = P[i];
            }

            for (int i=0;i<M.n;i++) {
                vector <double> cl1 = morph::Tools::getJetColor(P1[i]);
                displays[1].drawTriFill(M.X[i],M.X[M.N[i][0]],M.X[M.N[i][1]],cl1);
                displays[1].drawTriFill(M.X[i],M.X[M.N[i][3]],M.X[M.N[i][4]],cl1);
            }
            displays[1].redrawDisplay();
            break;
        }

        case 3:
        {

            std::stringstream frameFile1;
            frameFile1<<"logsWhole/"<<W.processName<<"frameNN";
            frameFile1<<setw(5)<<setfill('0')<<frameN;
            frameFile1<<".png";
            displays[0].saveImage(frameFile1.str());
            frameN++;
	    std::stringstream frameFile2;
            frameFile2<<"logsWhole/"<<W.processName<<"frameCC";
            frameFile2<<setw(5)<<setfill('0')<<frameM;
            frameFile2<<".png";
            displays[1].saveImage(frameFile2.str());
            frameM++;





            break;
        }

        case 4:
        {
            W.logfile<<W.processName<<"@"<<TIMEcs<<": 4=RECD"<<endl<<flush;
            switch(stoi(command[1]))
            {
            case 0: Dn = stod(command[2]); break;
            case 1: Dc = stod(command[2]); break;
            }
            break;
        }

        case 5: // *** SAVE ***
        {
            W.logfile<<W.processName<<"@"<<TIMEcs<<": 5=SAVE"<<endl<<flush;


            if (command.size()==2) {
                ofstream outFile;
                std::stringstream oFile; oFile<<command[1];
                outFile.open(oFile.str().c_str(),ios::out|ios::binary);
                double n=(double)S.size();double* N=&n;
                outFile.write((char*)N,sizeof(double));
                for (unsigned int i=0;i<S.size();i++) {
                    double m=(double)S[i].size();double* M=&m;
                    outFile.write((char*)M,sizeof(double));
                    for (unsigned int j=0;j<S[i].size();j++) {
                        outFile.write((char*)S[i][j],sizeof(double));
                    }
                }
                outFile.close();
            } else {
                W.logfile<<"No output filename."<<endl<<flush;
            }
                  
	         vector <int> angleDVector;
                 vector <int> radiusDVector;
                 vector <double> angleVector;
                 vector <double> radiusVector;
                 int degreeRadius;
                 int degreeAngle;
                 int numSectors = 12;
	          int radiusOffset = 11;
		  ofstream gfile ("edges.out");
                  angleDVector = M.sectorize_reg_Dangle(numSectors,radiusOffset, radiusOffset + 1);
                  //degreeAngle = M.find_max(angleVector,3);
                  degreeAngle = M.find_zeroDAngle(angleDVector);
                  gfile <<  " degreeDAngle "<< degreeAngle <<endl<<flush;

                  // analogue version
                  radiusOffset = 11;
                  angleVector = M.sectorize_reg_angle(numSectors,radiusOffset, radiusOffset + 1);
                  degreeAngle = M.find_zeroAngle(angleVector,3);
                  gfile <<  " degreeAngle "<< degreeAngle <<endl<<flush;
		  degreeAngle = M.find_maxTh(angleVector);
		  gfile <<  " degreeAngle via turnval "<< degreeAngle <<endl<<flush;
                  //radial degree
                  degreeRadius = -100;
                  int newdegreeRadius = 0;
                  for (int angleOffset=0; angleOffset<numSectors; angleOffset += 1){
                  radiusDVector = M.sectorize_reg_Dradius(numSectors, angleOffset, angleOffset + 3);
                  //gfile <<"after sectorize_reg_radius"<<endl;
                  // radiusVector = M.meanzero_vector(radiusVector);

                  newdegreeRadius = M.find_zeroDRadius(radiusDVector);
                  if (newdegreeRadius > degreeRadius)
                          degreeRadius = newdegreeRadius;
                  gfile << " ndR " << newdegreeRadius <<endl;
                  }
                  gfile  << " degreeRadius analogue  "<< degreeRadius << "  " <<endl << endl;
                  degreeRadius = -100;
                  newdegreeRadius = 0;
                  for (int angleOffset=0; angleOffset<numSectors; angleOffset += 1){
                  radiusVector = M.sectorize_reg_radius(numSectors, angleOffset, angleOffset + 3);
                  //gfile <<"after sectorize_reg_radius"<<endl;
                  // radiusVector = M.meanzero_vector(radiusVector);

                  newdegreeRadius = M.find_maxR(radiusVector);
                  // gfile << " ndR " << newdegreeRadius;
                  if (newdegreeRadius > degreeRadius)
                          degreeRadius = newdegreeRadius;
                  gfile << " ndR " << newdegreeRadius <<endl;
                  //degreeRadius = sumRadius / numSectors;
		  }
		  gfile  << " degreeRadius via max  "<< degreeRadius << "  " <<endl << endl;
                  ///radial degree
                  degreeRadius = -100;
                  newdegreeRadius = 0;
                  for (int angleOffset=0; angleOffset<numSectors; angleOffset += 1){
                  radiusVector = M.sectorize_reg_radius(numSectors, angleOffset, angleOffset + 3);
                  //gfile <<"after sectorize_reg_radius"<<endl;
                  // radiusVector = M.meanzero_vector(radiusVector);

                  newdegreeRadius = M.find_zeroRadius(radiusVector,3);
                  gfile << " ndR " << newdegreeRadius << endl;
                  if (newdegreeRadius > degreeRadius)
                                  degreeRadius = newdegreeRadius;
                  }


                  gfile <<  " degreeRadius "<< degreeRadius << "  " <<endl;


                  for (int k=0; k < (int) radiusVector.size();k++){
                      bfile<< setw(5) << " radius "<< radiusVector[k];
                      bfile<< setw(5) << " angle "<< angleVector[k];
                      bfile <<endl;
                  } //end of for on printing sector values

            break;
        }

        case 6: // *** LOAD ***
        {
            W.logfile<<W.processName<<"@"<<TIMEcs<<": 6=LOAD"<<endl<<flush;
	    W.logfile<<W.processName<<"value of n "<<M.n<<endl<<flush;

            if (command.size()==2) {
                double dummy;
                ifstream inFile;
                std::stringstream iFile;
                iFile<<command[1];
                inFile.open(iFile.str().c_str(),ios::in|ios::binary);
                inFile.read((char*)&dummy,sizeof(double));
                int I=(int)dummy;
                if (I==static_cast<int>(S.size())) {
                    for (int i=0;i<I;i++) {
                        inFile.read((char*)&dummy,sizeof(double));
                        int J=(int)dummy;
                        if (J==static_cast<int>(S[i].size())) {
                            for (int j=0;j<J;j++) {
                                inFile.read((char*)S[i][j],sizeof(double));
                            }
                        } else {
                            W.logfile<<"Wrong dims I."<<endl<<flush;
                        }
                    }
                } else {
                    W.logfile<<"Wrong dims J."<<endl<<flush;
                }
                inFile.close();
            } else {
                W.logfile<<"No input filename."<<endl<<flush;
            }



            break;
        }
	}
    }

    return 0;
};
