#include "morph/world.h"
#include "morph/sockserve.h"
#include "morph/display.h"
#include "morph/tools.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <math.h>
#include <random>
#include <algorithm>
#include <iomanip>
#include <boost/math/special_functions/bessel.hpp> 
#define PI 3.14159265
#define NUMPOINTS 82

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
  struct point {
    double xval;
    double yval;
  };

  // list of objects visible to member functions
  vector<vector<double> > X; //cartesian coordinates of each hex
  vector<vector<double> > H; // hex-grid info
  vector<vector<double> > N; // hex neighbourhood
  vector<vector<int> > region; //indices of all the regions
  vector<vector<double> >::iterator mydoubleiterator; //for iterating vectors of vectors
  vector<vector<int> >::iterator myintiterator;
  
  point  centres[NUMPOINTS]; //seed points
  vector<vector<double> > regionDist; //distances to each seed point
  vector<vector<int> > regionIndex; //indexes of each region
  vector<int> C; //count of neighbours
  vector<extremum> turnVal; //radial turning points

  vector<double> NN, CC; //hold the field values for each hex
  vector <bool> boundary;
  vector <bool> outerBoundary;
  //class constructor
  Erm2009 (int scale, int offset, double z) {

        ofstream afile ( "debug.out" );

	
	srand(time(NULL));
        this->scale = scale;
        this->offset = offset;


        H.resize(6);
	//centres.resize(NUMPOINTS);

        double s = pow(2.0, scale-1);
	ds = 29.0/s;
	//02-08-18 keep ds fixed now as region changes size.
	//ds = 0.453125;
	int inring = s;

       	// std::default_random_engine generator;
        // std::uniform_int_distribution<int> distribution(0,inring);
	// generator.min(0);
	// generator.max(inring);
	int red = 0; int blue = 0;
	double x = 0;
	double y = 0;
	double radius = 0;
	double theta = 0;
     for (int i=0; i<NUMPOINTS; i++) {
	    red = rand()%inring;
	    blue = rand()%inring;
	    if (blue%2 == 0){
	    x = (red-blue)*sqrt(3.)/(2.*s);
	    y = (1.*(red+blue))/(2.*s);
	    }
	    else{
	     x = -(red-blue)*sqrt(3.)/(2.*s);

	     y =  -(1.*(red+blue))/(2.*s);
	    }
	    centres[i].xval = x;
          centres[i].yval = y;
	 }


      vector<double> xy (2,0.);

     ///////////////

//     double sc= 1.6;
//     centres[0].xval=-0.06748*sc; centres[0].yval=0.3829*sc;
//     centres[1].xval=0.0589*sc; centres[1].yval=0.4221*sc;
//     centres[2].xval=0.1718*sc; centres[2].yval=0.4248*sc;
//     centres[3].xval=0.2699*sc; centres[3].yval=0.4171*sc;
//     centres[4].xval=0.2158*sc; centres[4].yval=0.3605*sc;
//     centres[5].xval=0.1081*sc; centres[5].yval=0.3532*sc;
//     centres[6].xval=0.0027*sc; centres[6].yval=0.327*sc;
//     centres[7].xval=-0.1679*sc; centres[7].yval=0.2928*sc;
//     centres[8].xval=-0.1041*sc; centres[8].yval=0.2645*sc;
//     centres[9].xval=0.0073*sc; centres[9].yval=0.2438*sc;
//     centres[10].xval=0.115*sc; centres[10].yval=0.2747*sc;
//     centres[11].xval=0.217*sc; centres[11].yval=0.2803*sc;
//     centres[12].xval=0.2239*sc; centres[12].yval=0.2094*sc;
//     centres[13].xval=0.12*sc; centres[13].yval=0.2024*sc;
//     centres[14].xval=-0.0*sc; centres[14].yval=0.1718*sc;
//     centres[15].xval=-0.1184*sc; centres[15].yval=0.1793*sc;
//     centres[16].xval=-0.2066*sc; centres[16].yval=0.1974*sc;
//     centres[17].xval=-0.2406*sc; centres[17].yval=0.0977*sc;
//     centres[18].xval=-0.1329*sc; centres[18].yval=0.1067*sc;
//     centres[19].xval=-0.0007*sc; centres[19].yval=0.113*sc;
//     centres[20].xval=0.1206*sc; centres[20].yval=0.1322*sc;
//     centres[21].xval=0.2155*sc; centres[21].yval=0.1462*sc;
//     centres[22].xval=0.2178*sc; centres[22].yval=0.0675*sc;
//     centres[23].xval=0.1678*sc; centres[23].yval=0.0911*sc;
//     centres[24].xval=0.155*sc; centres[24].yval=0.0307*sc;
//     centres[25].xval=0.207*sc; centres[25].yval=0.0081*sc;
//     centres[26].xval=0.2006*sc; centres[26].yval=-0.0459*sc;
//     centres[27].xval=0.1323*sc; centres[27].yval=-0.0934*sc;
//     centres[28].xval=0.0941*sc; centres[28].yval=-0.0589*sc;
//     centres[29].xval=0.0704*sc; centres[29].yval=-0.0162*sc;
//     centres[30].xval=0.0417*sc; centres[30].yval=0.0184*sc;
//     centres[31].xval=0.0208*sc; centres[31].yval=0.062*sc;
//     centres[32].xval=0.1067*sc; centres[32].yval=0.0615*sc;
//     centres[33].xval=-0.1276*sc; centres[33].yval=0.0408*sc;
//     centres[34].xval=-0.0974*sc; centres[34].yval=-0.0033*sc;
//     centres[35].xval=-0.0534*sc; centres[36].xval=-0.0286*sc;
//     centres[36].xval=-0.0107*sc; centres[36].yval=-0.0531*sc;
//     centres[37].xval=-0.0107*sc; centres[38].yval=-0.0531*sc;
//     centres[38].xval=0.0266*sc; centres[38].yval=-0.081*sc;
//     centres[39].xval=0.0619*sc; centres[39].yval=-0.1174*sc;
//     centres[40].xval=0.1274*sc; centres[40].yval=-0.1564*sc;
//     centres[41].xval=0.1692*sc; centres[41].yval=-0.167*sc;
//     centres[42].xval=0.1091*sc; centres[42].yval=-0.2025*sc;
//     centres[43].xval=0.0707*sc; centres[43].yval=-0.166*sc;
//     centres[44].xval=-0.0119*sc; centres[45].yval=-0.1357*sc;
//     centres[45].xval=-0.0663*sc; centres[45].yval=-0.12*sc;
//     centres[46].xval=-0.1154*sc; centres[46].yval=-0.1104*sc;
//     centres[47].xval=-0.169*sc; centres[48].yval=-0.1001*sc;
//     centres[48].xval=-0.1252*sc; centres[48].yval=-0.1853*sc;
//     centres[50].xval=-0.066*sc; centres[50].yval=-0.1891*sc;
//     centres[51].xval=0.0035*sc; centres[51].yval=-0.1907*sc;
//     centres[52].xval=0.0596*sc; centres[52].yval=-0.2268*sc;
//     centres[53].xval=0.1044*sc; centres[53].yval=-0.2541*sc;
//     centres[54].xval=0.0685*sc; centres[54].yval=-0.2867*sc;
//     centres[55].xval=0.1052*sc; centres[55].yval=-0.3137*sc;
//     centres[56].xval=0.0965*sc; centres[56].yval=-0.3544*sc;
//     centres[57].xval=0.0173*sc; centres[57].yval=-0.317*sc;
//     centres[58].xval=0.0337*sc; centres[58].yval=-0.3674*sc;
//     centres[60].xval=-0.0242*sc; centres[60].yval=-0.3959*sc;
//     centres[61].xval=-0.0388*sc; centres[61].yval=-0.3435*sc;
//     centres[62].xval=-0.0593*sc; centres[62].xval=-0.2815*sc;
//     centres[63].xval=0.0177*sc; centres[63].yval=-0.2543*sc;
//     centres[64].xval=-0.0351*sc; centres[64].yval=-0.2405*sc;
//     centres[65].xval=-0.1236*sc; centres[65].yval=-0.2457*sc;
//     centres[66].xval=-0.2093*sc; centres[66].yval=-0.2598*sc;
//     centres[67].xval=-0.1615*sc; centres[67].yval=-0.2907*sc;
//     centres[68].xval=-0.1096*sc; centres[68].yval=-0.3126*sc;
//     centres[69].xval=-0.1013*sc; centres[69].yval=-0.3663*sc;
//     centres[70].xval=-0.1664*sc; centres[70].yval=-0.3576*sc;
//     centres[71].xval=-0.2292*sc; centres[71].yval=-0.3302*sc;
//     centres[72].xval=-0.2699*sc; centres[72].yval=-0.2974*sc;
//     centres[73].xval=-0.2614*sc; centres[73].yval=-0.3801*sc;
//     centres[74].xval=-0.2185*sc; centres[74].yval=-0.4057*sc;
//     centres[75].xval=-0.1629*sc; centres[75].yval=-0.4205*sc;
//     centres[76].xval=-0.0968*sc; centres[76].yval=-0.4248*sc;
//     centres[77].xval=-0.2009*sc; centres[77].yval=-0.0728*sc;
//     centres[78].xval=-0.2245*sc; centres[78].yval=-0.0243*sc;
//     centres[79].xval=-0.2499*sc; centres[79].yval=0.0177*sc;
//     centres[80].xval=-0.239*sc; centres[80].yval=-0.1865*sc;
//     centres[81].xval=-0.1812*sc; centres[81].yval=-0.183*sc;
//	
//	 centres[4].xval = 0.0; centres[0].yval = 0.0;
//	 centres[3].xval = -0.1; centres[3].yval = 0.5;
	// centres[2].xval = -0.5; centres[2].yval = 0.5;
	// centres[1].xval = -0.5; centres[1].yval = -0.5;
	// centres[0].xval = -0.1; centres[0].yval = -0.5;
	//for (int j=0;j<NUMPOINTS;j++)
        //  afile << "x = " << centres[j] .xval << "  y = " << centres[j].yval <<endl;

        n = 0;
	x = 0;
        y = 0;
        radius = 0;
        theta = 0;
	for (int r=-3*s/2; r<0; r++) {
	  for (int b=-3*s/2; b<=3*s/2; b++) {
	 	      x = (r-b)*sqrt(3.)/(2.*s);
	 	      y = (1.*(r+b))/(2.*s);
		      radius = 29.0*sqrt(x*x+y*y);
		      if (y >= 0) 
		        theta = atan2(y,x);	
		      else
			theta = 2*PI + atan2(y,x);
		      
	  
		      if (x*x + y*y <= 0.75) { // circular
                      //if (x*x/0.1 + y*y/0.8 <= 1) { // ellipse 1
                                 //if (x*x/0.8 + y*y/0.1 <= 1) { // ellipse 2
                                 //if (x*x/0.8 + y*y/0.01 <= 0.9) { // ellipse 3
                                 //if (x*x/0.01 + y*y/0.8 <= 0.9) { // ellipse 4
                    //  if (x*x < 0.25 && y*y < 0.4) { // square 1
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
	 

	for (int r=0; r<=3*s/2; r++) {
	  for (int b=-3*s/2; b<=3*s/2;b++){
	 	      y = (1.*(r+b))/(2.*s);
		      x = (r-b)*sqrt(3.)/(2.*s);
		      
		      radius = 29.0*sqrt(x*x+y*y);
		      if (radius == 0)
			theta = 0;
		      else
		        if (y >= 0) 
			  theta = atan2(y,x);  
		        else
			  theta = 2*PI + atan2(y,x);
		      if (x*x + y*y <= 0.75) { // circular
                      //    if (x*x/0.1 + y*y/0.8 <= 1) { // ellipse 1
                                  //if (x*x/0.8 + y*y/0.1 <= 1) { // ellipse 2
                                  //if (x*x/0.8 + y*y/0.01 <= 0.9) { // ellipse 3
                                  //if (x*x/0.01 + y*y/0.8 <= 0.9) { // ellipse 4
                     // if (x*x < 0.25 && y*y < 0.4) { // square 1
                                  //if (x*x < 0.2 && y*y < 0.2) { // square 2
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
	 
	

//these are the vectors of vectors for the regions
	regionDist.resize(n);
	region.resize(n);

        

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
                    //anticlockwise froafile<<endl;m east
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
	int boundaryCount = 0;
        for (int i=0;i<n;i++){
	 if (C[i] == 6) {
	 outerBoundary.push_back(0);
	 }
	 else {
	 outerBoundary.push_back(1);
	 boundaryCount++;
	 }
       }

       cout <<"number of boundary elements  " <<boundaryCount<<endl;
     

	//X is used in plotting, potentially also for stacked grids hence 3-vector
         X.resize(n);
        for(int i=0;i<n;i++){
	  X[i].resize(3,0.);
            X[i][0] = H[0][i];
            X[i][1] = H[1][i];
	}
	
	//arrays to hold the field values
        NN.resize(n);
        CC.resize(n);
	

    for (int i=0;i<n;i++){
      point hexcen;
      hexcen.xval = H[0][i];
      hexcen.yval = H[1][i];
    for (int j=0;j<NUMPOINTS;j++){
      double temp = this->getdist(hexcen,centres[j]);
      // afile << "   " <<hexcen.xval <<"  "<<hexcen.yval;
      regionDist[i].push_back(temp);
      // afile << temp <<"  ";
     }
    //afile << endl;
   }
    
    // to produce a list of the sorted distances of each hex from the seed points
    vector <vector <double> > sortedDist;
    sortedDist.resize(n);

    for (int i=0;i<n;i++) {
     	vector <double> tempvector1;
        tempvector1 = regionDist[i];
     	std::sort(tempvector1.begin(),tempvector1.end());
     	sortedDist[i] = tempvector1;
        vector <int> tempint = sort_indexes(regionDist[i]);
        region[i] = tempint;
       
    }
     afile <<"list of all the regions"<<endl;
     for(int i=0;i<n;i++){
      
       if (sortedDist[i][0] == sortedDist[i][1]){
	    afile<<"i ="<<i;
     	   afile << "  "<<region[i][0]<<"  "<< region[i][1];
       }
         afile<<endl;
      }

     // this determines the internal boundaries
     for(int i=0;i<n;i++){
        for(int j=0;j<6;j++){
	   if ((region[i][0] !=region[N[i][j]][0])) {
	      N[i][j] = i;
	      C[i]--;
	   }
	}
     }

    //information about the indexes of hexes in regions  
     regionIndex.resize(NUMPOINTS);
      for (int j=0;j<NUMPOINTS;j++) {
        for (int i=0;i<this->n;i++){
	  if (region[i][0] == j)
      	   regionIndex[j].push_back(i);
        }
      }

      //this marks the internal boundaries
     for (int i=0;i<n;i++){
       if (C[i] == 6)
	 boundary.push_back(0);
       else
	 boundary.push_back(1);
     }
    
    
     
  } //end of constructor

  //function, gets distance between points
  double getdist(point a, point b) {
    double result;
    result = sqrt((a.xval-b.xval)*(a.xval-b.xval) + (a.yval-b.yval)*(a.yval-b.yval));
    return result;
  }

  //function, finds the indexes from original list in sorted order
template <typename T>
vector<int> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<int> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](int i1, int i2) {return v[i1] < v[i2];});

  return idx;
}

  //function derives Laplacian for a given field 
    vector<double> getLaplacian(vector<double> Q, double dx) {
      double overdxSquare = 1./(1.5*dx*dx);
        vector<double> L(n,0.);
        for(int i=0;i<n;i++){
            L[i]=(Q[N[i][0]]+Q[N[i][1]]+Q[N[i][2]]+Q[N[i][3]]+Q[N[i][4]]+Q[N[i][5]]-6.*Q[i])*overdxSquare;
        }
        return L;
    }

  //function to timestep coupled equations
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


  //function find_max to find turning points both values and indices.
  int find_max(vector<double> ray) {
    int iend = ray.size();
    ofstream dfile ("turn.txt");
    dfile <<"iend = " << iend <<endl;
    
    turnVal.resize(1000);
    cout <<" "<<iend<<iend<<flush;
    double old_slope = 0;
    double new_slope = 0;
    int count = 0;
    old_slope = ray[1] - ray[0];
    for (int i =2; i<=iend+1;i++){
      
      new_slope = ray[i%iend]-ray[(i-1)%iend];
      dfile << " " << i%iend << " " << old_slope << " "<<new_slope <<endl;
      if (new_slope*old_slope < 0.) {
	turnVal[count].radialIndex = i;
	turnVal[count].radialValue = ray[i]; //should really interpolate
        count++;
      }
      old_slope = new_slope;
    }
    return count;
  }				  
    

  //function bessel_ray for computing bessel functions along a ray
  // takes a vector of doubles representing the radius
  // returns the Bessel function values
  // vector<double> bessel_ray (int v, vector<double> ray) {
  //   vector <double> result(n,0.);
  //   result = boost::cyl_bessel_j(v, ray);
  //   return result;
  // }

 
  //function to return area of a region
  double regArea (int regNum) {
    double area = 0;
    for (int i=0;i < (int) regionIndex[regNum].size();i++){
      if (outerBoundary[regionIndex[regNum][i]] == 1){
        area = 0;
        break;
      	}
        else {
	  area += 1.;
	}
    }
    return area;
  } //end of funtion regArea

  //function to return perimeter of a region
  double regPerimeter (int regNum) {
    double perimeter = 0;
    for (int i=0;i < (int) regionIndex[regNum].size();i++)
      if (this->boundary[regionIndex[regNum][i]] == 1)
	perimeter += 1.0;
    return perimeter;
  } //end of function regPerimeter


// function to give r and theta relative to region centre
  vector <double> set_polars(int regNum){
     vector<double> result;
     result.resize(2);
     double xcentre = this->centres[regNum].xval;
     double ycentre = this->centres[regNum].yval;
     double xav=0;
     double yav = 0;
     int hexcount = 0;

     for (int i=0;i< (int) this->regionIndex[regNum].size();i++) {
       hexcount++;
       xav += this->H[0][regionIndex[regNum][i]];
       yav += this->H[1][regionIndex[regNum][i]];
     }
     xav = xav / hexcount;
     yav = yav / hexcount;
       
     //go over the region and put the hexes into bins then average 
     for (int i=0;i< (int) this->regionIndex[regNum].size();i++) {
       int index = regionIndex[regNum][i];
       this->H[4][index] = sqrt((this->H[0][index]-xav)*(this->H[0][index]-xav) + \
				(this->H[1][index]-yav)*(this->H[1][index]-yav));
       if (this->H[1][i] >= 0)
	  this->H[5][index] = atan2((this->H[1][index]-yav), (this->H[0][index]-xav));
       else
	 this->H[5][index] = PI + atan2((this->H[1][index]-yav), (this->H[0][index]-xav));
     }

     result[0] = xav - xcentre; //diff between seed point and barycentre
     result[1] = yav - ycentre;
     return result;
  } //end of function set_polars


  //function to count the hexes in sectors of a region
  vector <double> sectorize_region (int regNum) {
    ofstream cfile ("logs/sector.txt");
    vector <double> diff; //difference between seed point and CoG of region
    diff.resize(2);
    diff = set_polars(regNum);
    vector <double> totalCC; //average value of CC in each sector
    vector <double>  angleVal;
    vector <int> count; //number of hexes in each sector
    int numSectors = 12; // number of sectors
    totalCC.resize(numSectors);
    count.resize(numSectors);
    double startAngle, endAngle, angleInc; //sector angles
    angleInc = 2*PI/(1.*numSectors);
    for (int k=0;k<numSectors;k++) {
      //double angle;
      startAngle = k*angleInc;
      endAngle = (k+1)*angleInc;
      if ((k+1) == numSectors)
	 endAngle = 2*PI;
      int size = (int) regionIndex[regNum].size();
      
      for (int i=0; i< size;i++) {
	if (this->H[5][regionIndex[regNum][i]]>=startAngle && this->H[5][regionIndex[regNum][i]]<endAngle) {
	  angleVal.push_back(H[5][regionIndex[regNum][i]]);
	  count[k]++;
	  totalCC[k] += this->CC[regionIndex[regNum][i]];
	}
	//cfile << setw(5) << angleVal[i]  <<"  ";
      }
      // if (k == 0) {
      //   std::sort(angleVal.begin(),angleVal.end());
      //   //for (int i=0; i< (int) angleVal.size();i++) {
      //  	for (int i=0; i<2;i++) {
      //   cfile << setw(5) << angleVal[i]  <<endl;
      // 	//cfile << setw(5) << this->H[5][regionIndex[regNum][5]]  <<endl;
      // 	}
      // }
      cfile << endl;
      cfile << "diff x " << diff[0] << "diff y " << diff[1] << endl;
      cfile << " angleinc " <<angleInc << "  start angle "  << startAngle << "  end angle "<< endAngle <<endl;
      if (count[k] != 0)
	totalCC[k] = totalCC[k] / (1.*count[k]) - 2.5;
      else
	totalCC[k] = -999.999;
    }//end loop on k
     return totalCC;
    
  } //end of function sectorize_region
    
       
}; // Erm2009


int main (int argc, char **argv)
{
    srand(atoi(argv[3]));
    ofstream bfile ( "logs/sector.out" );
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


    Erm2009 M(8,0,0.);
   double  value;
  vector <double> ray;
   for (int i=0;i<M.n;i++) {
     value  = boost::math::cyl_bessel_j(0,M.H[4][i]);
       ray.push_back(value);
  }

 

   

   
 
   for (int i=0;i<M.n;i++) {
     M.NN[i]=(morph::Tools::randDouble())*0.1 +1.;
     M.CC[i]=(morph::Tools::randDouble())*0.1+2.5;
     //        M.CC[i]=sin(i*3.142)*0.1+2.5;
     //	 M.CC[i]=sin(i*3.142)*0.1+2.5;
     }

     // for (int i=0;i<M.n;i++) {
     //   // M.NN[i]= 1. + sin(i*3.142)*sin(i*3.142)*0.001;
     //    M.NN[i]= 1.;
     //   M.CC[i]=ray[i]*0.1+2.5;
     //  }
   
    unsigned int frameN = 0;
    unsigned int frameM = 0;
    //initial values for Dn Dc
    double Dn = 100.;
    double Dc = Dn*0.3;
    double TIME=0.;
    vector<double*> f;
    f.push_back(&TIME);
    
    int vdegree = 0;
    // int mdegree = 0;
     vector <double> rayv;
     vector <double> ringv;
     
   
 
  
  

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
	  //W.logfile<<W.processName<<"@"<<TIMEcs<<": 0=QUIT"<<endl<<flush;
	  // W.logfile.close();
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
//start of code for determining colour scale
	     for (int j=0;j<NUMPOINTS;j++){
	       vector<double> plt;
	       int index;
	       double tempNN;
	       double maxV = -1e7;
	       double minV = +1e7;
	       int size= (int) M.regionIndex[j].size();
	       for (int i=0;i<size ;i++){
	         index = M.regionIndex[j][i];
		 tempNN = M.NN[index];
		 plt.push_back(tempNN);
           
		 //for (int i=0;i<M.n;i++) {
                if (M.C[index]==6) {
		  if (tempNN>maxV) { maxV = tempNN; }
                    if (tempNN<minV) { minV = tempNN; }  
		}//if on M.C
	       }//if on i, going through the region
	       
            double scaleV = 1./(maxV-minV);
            vector<double> P(size,0.);
            for (int i=0;i<size;i++) {
                P[i] = fmin (fmax (((plt[i])-minV)*scaleV,0.),1.);
                // M.X[i][2] = P[i];
            }

            for (int i=0;i<size;i++) {
	      index = M.regionIndex[j][i];
                vector <double> cl = morph::Tools::getJetColor(P[i]);
                displays[0].drawTriFill(M.X[index],M.X[M.N[index][0]],M.X[M.N[index][1]],cl);
		displays[0].drawTriFill(M.X[index],M.X[M.N[index][3]],M.X[M.N[index][4]],cl);
            } //end of loop on i
            
	    }//end of loop on j over regions
	    displays[0].redrawDisplay();//this MUST be outside the j loop so the whole image is
	                                 //build BEFORE it is displayed

	    //start of code for determining colour scale CC field
	    for (int j=0;j<NUMPOINTS;j++){
	       vector<double> plu;
	       int index1;
	       double tempCC;
	       double maxV1 = -1e7;
	       double minV1 = +1e7;
	       int size1= (int) M.regionIndex[j].size();
	       for (int i=0;i<size1 ;i++){
	         index1 = M.regionIndex[j][i];
		 tempCC = M.CC[index1];
		 plu.push_back(tempCC);
           
		 //for (int i=0;i<M.n;i++) {
                if (M.C[index1]==6) {
		  if (tempCC>maxV1) { maxV1 = tempCC; }
                    if (tempCC<minV1) { minV1 = tempCC; }    
                }//if on M.C
	       }//if on i, going through the region
	       
            double scaleV1 = 1./(maxV1-minV1);
            vector<double> P1(size1,0.);
            for (int i=0;i<size1;i++) {
                P1[i] = fmin (fmax (((plu[i])-minV1)*scaleV1,0.),1.);
                // M.X[i][2] = P[i];
            }

            for (int i=0;i<size1;i++) {
	      index1 = M.regionIndex[j][i];
                vector <double> cl1 = morph::Tools::getJetColor(P1[i]);
                displays[1].drawTriFill(M.X[index1],M.X[M.N[index1][0]],M.X[M.N[index1][1]],cl1);
		displays[1].drawTriFill(M.X[index1],M.X[M.N[index1][3]],M.X[M.N[index1][4]],cl1);
            } //end of loop on 
	    }//end of loop on j over regions
            displays[1].redrawDisplay(); //this MUST be outside the j loop so the whole image is
	                                 //build BEFORE it is displayed
	    
	    // vector<double> plu = M.CC;
            // double maxV1 = -1e7;
            // double minV1 = +1e7;
            // for (int i=0;i<M.n;i++) {
            //     if (M.C[i]==6) {
            //         if (plu[i]>maxV1) { maxV1 = plu[i]; }
            //         if (plu[i]<minV1) { minV1 = plu[i]; }
            //     }
            // }
            // double scaleV1= 1./(maxV1-minV1);
            // vector<double> P1(M.n,0.);
            // for (int i=0;i<M.n;i++) {
            //     P1[i] = fmin (fmax (((plu[i])-minV1)*scaleV1,0.),1.);
            //     // M.X[i][2] = P[i];
            // }
	    // //end of code for determining colour scale

	    // //code for plotting, note this is on X not on H
            // for (int i=0;i<M.n;i++) {
            //     vector <double> cl1 = morph::Tools::getJetColor(P1[i]);
            //     displays[1].drawTriFill(M.X[i],M.X[M.N[i][0]],M.X[M.N[i][1]],cl1);
            //     displays[1].drawTriFill(M.X[i],M.X[M.N[i][3]],M.X[M.N[i][4]],cl1);
            // }
            // displays[1].redrawDisplay();
	    
            break;
        }

        case 3:
        {
	  
            std::stringstream frameFile1;
            frameFile1<<"logs/"<<W.processName<<"frameNN";
            frameFile1<<setw(5)<<setfill('0')<<frameN;
            frameFile1<<".png";
            displays[0].saveImage(frameFile1.str());
            frameN++;
	    std::stringstream frameFile2;
            frameFile2<<"logs/"<<W.processName<<"frameCC";
            frameFile2<<setw(5)<<setfill('0')<<frameM;
            frameFile2<<".png";
            displays[1].saveImage(frameFile2.str());
            frameM++;


            break;
        }

        case 4:
        {
	  //W.logfile<<W.processName<<"@"<<TIMEcs<<": 4=RECD"<<endl<<flush;
	   switch(stoi(command[1]))
            {
            case 0: Dn = stod(command[2]); break;
            case 1: Dc = stod(command[2]); break;
            }
            break;
        }

        case 5: // *** SAVE ***
        {
	  //W.logfile<<W.processName<<"@"<<TIMEcs<<": 5=SAVE"<<endl<<flush;
	  //W.logfile<<W.processName<<"@"<<TIMEcs<<": 6=LOAD"<<endl<<flush;
	  //W.logfile<<W.processName<<"value of n "<<M.n<<endl<<flush;
	    // for (int i=0; i<M.n;i++)
	    // W.logfile<<"i =  "<<i<<"region = "<<M.region[0][i]<<endl<<flush;

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
	      //    W.logfile<<"No output filename."<<endl<<flush;
            }

	    
            int regionCount = 0;
	    vector <double> tempvector;
              for (int j=0;j<NUMPOINTS;j++) {
	      if (M.regArea(j) != 0){
	      tempvector = M.sectorize_region(j);
	      vdegree = M.find_max(tempvector);
	      W.logfile<<W.processName<<" "<<vdegree<<"  "<<M.regArea(j)<<"  "<<M.regPerimeter(j)<<endl<<flush;
	      regionCount++;
		
	      //bfile <<"region "<<j<<" ";
		for (int k=0; k <= (int) tempvector.size();k++){
                bfile<< setw(5) << " "<< tempvector[k%tempvector.size()];
		bfile <<endl;
		} //end of for on printing sector values
	      } //end of if on non-zero regions
		bfile <<regionCount<<endl;
	      } //end of loop on NUMPOINTs
	    
	    
            // mdegree = M.find_max(ringv);
	    
            // W.logfile<<W.processName<<"  mdegree= "<<mdegree<<endl<<flush;
	    // W.logfile<<W.processName<<"  rad  "<<rad<<endl<<flush;
            // W.logfile<<W.processName<<"  ringvsize= "<<temp<<endl<<flush;
	    // W.logfile<<W.processName<<"  ds "<<dstemp<<M.ds<<endl<<flush;

            break;
        }

        case 6: // *** LOAD ***
        {
	  //W.logfile<<W.processName<<"@"<<TIMEcs<<": 6=LOAD"<<endl<<flush;
	  // W.logfile<<W.processName<<"value of n "<<M.n<<endl<<flush;
	    // for (int i=0; i<M.n;i++)
	    //   W.logfile<<"i =  "<<i<<"region = "<<M.region[0][i]<<endl<<flush;
	    
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
	
	
	    
	    //  rayv.resize(0);
	    //  ringv.resize(0);
	    //  double dstemp = M.ds;

            //    for (int i=0;i<M.n;i++)
            // // //  //if ((M.H[3][i] == 0) && (M.H[2][i] >= 0)) //positive r axis
            //      if (M.H[3][i] == 0) //positive r axis
            //         rayv.push_back(abs(M.NN[i]));
   

            //    vdegree = M.find_max(rayv);
  
           
            //  for (int i=0;i<M.n;i++) 
            //    if ((M.H[4][i] <= rad + dstemp/2) && (M.H[4][i] >= rad - dstemp/2))
            //       ringv.push_back(abs(M.NN[i]));
	    //  int temp = ringv.size();
	    

            // mdegree = M.find_max(ringv)
	      ;
	    //W.logfile<<W.processName<<"  vdegree= "<<vdegree<<endl<<flush;
            // W.logfile<<W.processName<<"  mdegree= "<<mdegree<<endl<<flush;
	    // W.logfile<<W.processName<<"  rad  "<<rad<<endl<<flush;
            // W.logfile<<W.processName<<"  ringvsize= "<<temp<<endl<<flush;
	    // W.logfile<<W.processName<<"  dstemp "<<dstemp<<dstemp<<endl<<flush;


	    
            break;
        }
	}
    }

    return 0;
};
