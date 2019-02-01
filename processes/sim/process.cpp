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
#include <map>
#include <boost/math/special_functions/bessel.hpp>
#define PI 3.14159265
#define NUMPOINTS 78  //this is after deleting point 73

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
    //   vector<vector<double> > Nold; // holds values before resetting b.c.s
  vector<vector<int> > region; //indices of all the regions
  vector<vector<int>> hexRegionList; //list of the different regions around a hex
  vector<vector<int>> regionList; // list of the regions around a region
  std::map<int,vector<int>> edges; //map of (i,j) edges, uses pair<->int converters
    // std::map<int,double> > corrAdjacent;

  point  centres[NUMPOINTS]; //seed points
  vector<vector<double> > regionDist; //distances to each seed point
  vector<vector<int> > regionIndex; //indexes of hexes in each region
  vector<int> C; //count of neighbours
  vector<int> Creg; //count of regions
  vector<extremum> turnVal; //radial turning points
  vector<std::pair<double,double>> diff; //difference between seed point and CoG of region

  vector<double> NN, CC; //hold the field values for each hex
  vector <bool> boundary;
  vector <bool> outerBoundary;
  int base = 1000;
  //class constructor
  Erm2009 (int scale, int offset, double z) {
        ofstream afile ( "debug.out" );
        ofstream jfile("correlateVector.out");

	srand(time(NULL));
        this->scale = scale;
        this->offset = offset;
// #include "barrel.h"

        H.resize(6);
        double s = pow(2.0, scale-1);
	ds = 29.0/s;
	//02-08-18 keep ds fixed now as region changes size.
	//ds = 0.453125;
	int inring = s;

       	// std::default_random_engine generator;
        // std::uniform_int_distribution<int> distribution(0,inring);
        // this->generator.min(0);
        // this->generator.max(inring);
	// int red = 0; int blue = 0;
	 double x = 0;
	 double y = 0;
	 double radius = 0;
	 double theta = 0;
        // for (int i=0; i<NUMPOINTS; i++) {
     	//     red = rand()%inring;
     	//     blue = rand()%inring;
     	//     if (blue%2 == 0){
        //         blue = -blue;
     	//     }
     	//     if (red%2 == 0) {
        //         red = -red;
     	//     }
        //     x = (red-blue)*sqrt(3.)/(2.*s);
     	//     y = (1.*(red+blue))/(2.*s);
     	//     centres[i].xval = x;
        //     centres[i].yval = y;
     	//  }

        for (int i=0; i<NUMPOINTS;i++) {
            centres[i].xval = 0.88*(rand()/(RAND_MAX +1.0) - 0.5);
            centres[i].yval = 1.34*(rand()/(RAND_MAX + 1.0) - 0.5);
        }



      vector<double> xy (2,0.);

     ///////////////


//	 centres[4].xval = 0.0; centres[0].yval = 0.0;
//	 centres[3].xval = -0.1; centres[3].yval = 0.5;
	// centres[2].xval = -0.5; centres[2].yval = 0.5;
	// centres[1].xval = -0.5; centres[1].yval = -0.5;
	// centres[0].xval = -0.1; centres[0].yval = -0.5;
	for (int j=0;j<NUMPOINTS;j++)
          afile << "j = " << j <<" x = " << centres[j] .xval  << "  y = " << centres[j].yval <<endl;

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


		      //if (x*x + y*y <= 0.75) { // circular
                      //if (x*x/0.1 + y*y/0.8 <= 1) { // ellipse 1
                                 //if (x*x/0.8 + y*y/0.1 <= 1) { // ellipse 2
                                 //if (x*x/0.8 + y*y/0.01 <= 0.9) { // ellipse 3
                                 //if (x*x/0.01 + y*y/0.8 <= 0.9) { // ellipse 4
                      if (x*x < 0.2 && y*y < 0.45) { // square 1
                                 //if (x*x < 0.2 && y*y < 0.2) { // square 2
		      //if (x*x/0.05 - y*y/0.05 <= 0.7) { // hyperbolic 1
                              H[0].push_back(x);                      //X
                              H[1].push_back(y);                      //Y
                              H[2].push_back(r);                      //r
                              H[3].push_back(b);                      //b
                              H[4].push_back(0);                 //r
                              H[5].push_back(0);                  //theta

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
		      // if (x*x + y*y <= 0.75) { // circular
                      //    if (x*x/0.1 + y*y/0.8 <= 1) { // ellipse 1
                                  //if (x*x/0.8 + y*y/0.1 <= 1) { // ellipse 2
                                  //if (x*x/0.8 + y*y/0.01 <= 0.9) { // ellipse 3
                                  //if (x*x/0.01 + y*y/0.8 <= 0.9) { // ellipse 4
                       if (x*x < 0.2 && y*y < 0.45) { // square 1
                                  //if (x*x < 0.2 && y*y < 0.2) { // square 2
		      //if (x*x/0.05 - y*y/0.05 <= 0.7) { // hyperbolic 1
                               H[0].push_back(x);                      //X
                               H[1].push_back(y);                      //Y
                               H[2].push_back(r);                      //r
                               H[3].push_back(b);                      //b
                               H[4].push_back(0);                      //r
                               H[5].push_back(0);                  //theta

                               n++;
		      } //if on shape determination
	  } //if on b

	} //if on r



//these are the vectors of vectors for the regions
	regionDist.resize(n);
	region.resize(n);



        // get neighbours
        N.resize(n); //neighbouring hexes
        //Nold.resize(n); //holds values before b.c.s
        C.resize(n,0); //count of neighbouring hexes
        Creg.resize(n,0);
        hexRegionList.resize(n); //neighbouring regions for a hex
        regionList.resize(n); //neighbouring regions for a region

        for(int i=0;i<n;i++){
            N[i].resize(6,i); // CONNECT ALL TO 'boundary' UNIT AT N+1
            for(int j=0;j<n;j++){
                int dr = H[2][j]-H[2][i];
		// int dg = H[3][j]-H[3][i];
                int db = H[3][j]-H[3][i];
		// if(max(max(abs(dr),abs(dg)),abs(db))==1){
		   if(max(abs(dr),abs(db))==1){
                    //anticlockwise froafile<<endl;m east
                       if(db==0&&dr==+1){N[i][0]=j;C[i]++;} //Nold[i][0]=j;C[i]++;}
                       if(db==+1&&dr== +1){N[i][1]=j;C[i]++;}//Nold[i][1];C[i]++;}
                       if(db== +1&&dr==0){N[i][2]=j;C[i]++;}//Nold[i][2];C[i]++;}
                       if(db==0&&dr==-1){N[i][3]=j;C[i]++;}//Nold[i][3];C[i]++;}
                       if(db==-1&&dr==-1){N[i][4]=j;C[i]++;}//Nold[i][4];C[i]++;}
                       if(db== -1&&dr==0){N[i][5]=j;C[i]++;}//Nold[i][5];C[i]++;}
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

       afile <<"number of boundary elements  " <<boundaryCount<<endl;


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

	// get a list of distances from each Dirichlet point for each hex.
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
    // afile << endl;
    }

    // to produce a list of the sorted distances of each hex from the seed points
    // Note the process of sorting the indexes in distance order is completely independent
    // from the production of Sorted list but we assume that the similar sorting used for
    // each means that they are compatible.
    vector <vector <double> > sortedDist;
    sortedDist.resize(n);

    for (int i=0;i<n;i++) {
     	vector <double> tempvector1;
        tempvector1 = regionDist[i];
     	std::stable_sort(tempvector1.begin(),tempvector1.end());
     	sortedDist[i] = tempvector1;
        vector <int> tempint = sort_indexes(regionDist[i]);
        region[i] = tempint;

    }

     afile <<"list of all the regions"<<endl;
     for(int i=0;i<n;i++){

       if (sortedDist[i][0] == sortedDist[i][1]){
           //afile<<"i ="<<i;
     	   //afile << "  "<<region[i][0]<<"  "<< region[i][1];
       }
     }


     // this determines the internal boundaries
     // for(int i=0;i<n;i++){
     //    for(int j=0;j<6;j++){
     //       if ((region[i][0] !=region[N[i][j]][0])) {
     //          N[i][j] = i;
     //          C[i]--;
     //       }
     //    }
     // }

     // this determines the type of hex
      for (int i=0;i < n ;i++){
          if(outerBoundary[i] == true)
              continue;
          int centralRegion = region[i][0];
          int oldRegion = centralRegion;
          int newRegion = 0;

          for(int j=0;j<6;j++) {
              hexRegionList[i].push_back(region[N[i][j]][0]); //push back the region of each neighbour
              newRegion =  region[N[i][j]][0];
              if (centralRegion != newRegion){ //its a boundary hex
                  C[i]--;
                  N[i][j] = i;
                  if (oldRegion != newRegion){ // its a vertex
                      Creg[i]++;
                      oldRegion = newRegion;
                  }
              }
          }
      } //end of logic determining hex types






    //information about the indexes of hexes in regions
     regionIndex.resize(NUMPOINTS);
      for (int j=0;j<NUMPOINTS;j++) {
        for (int i=0;i<this->n;i++){
            if (region[i][0] == j)
      	   regionIndex[j].push_back(i);
        }
      }



      //set the difference between seed point and CoG of region
      diff.resize(NUMPOINTS);
      for (int j=0;j<NUMPOINTS;j++){
          diff[j] = this->set_polars(j);
          afile << "diff seed-centre" << diff[j].first << " " << diff[j].second<<endl;
      }
  } //end of constructor

  //function, gets distance between points
  double getdist(point a, point b) {
    double result;
    result = sqrt((a.xval-b.xval)*(a.xval-b.xval) + (a.yval-b.yval)*(a.yval-b.yval));
    return result;
  }

  //function, finds the indexes from original list in sorted order
  // now using stable_sort to ensure consistency in case of hexes equidistant between
  // centres
template <typename T>
vector<int> sort_indexes(const vector<T> &v) {

  // initialize original index locations, iota assigns integer values from 0 onwards
  // to the vector of ints idx
  vector<int> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  stable_sort(idx.begin(), idx.end(),
       [&v](int i1, int i2) {return v[i1] < v[i2];});

  return idx;
}

     int pair2int (pair<int,int> d, int ibase) {
        int result;
        result = d.first*ibase + d.second;
        return result;
    }

    std::pair<int,int> int2pair(int c, int ibase) {
        std::pair<int, int> result;
        result.first = c / ibase;
        result.second = c%ibase;
        return result;
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
      dt = dt * 2.5 / Dn;


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

    // find the zeros in a ray
    int find_zero(vector<double> ray) {
    int iend = ray.size();
    ofstream zerofile ("zero.txt");
    zerofile <<"iend = " << iend <<endl;
    double sum = 0;
    // to normalise the ray round 0
    for (int i=0;i<iend;i++)
        sum += ray[i];
    sum = sum /(iend*1.0);
    for (int i=0;i<iend;i++)
        ray[i] -= sum;

    turnVal.resize(1000);
    cout <<" "<<iend<<iend<<flush;
    double old_val = 0;
    double new_val = 0;
    int count = 0;
    old_val = ray[0];
    for (int i =1; i<iend+1;i++){
      new_val = ray[i%iend];
      zerofile << " " << i%iend << " " << old_val << " "<<new_val <<endl;
      if (new_val*old_val < 0.) {
	turnVal[count].radialIndex = i;
	turnVal[count].radialValue = ray[i]; //should really interpolate
        count++;
      }
      old_val = new_val;
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
      if (outerBoundary[regionIndex[regNum][i]] == true){
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
       if (Creg[regionIndex[regNum][i]] > 0)
	perimeter += 1.0;
    return perimeter;
  } //end of function regPerimeter

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
        }
        sum = sum/(1.0*size);
        absSum = absSum / (1.0 * size);
        meanzero << " mean  " << sum << endl;
        meanzero << " absolute mean  " << absSum << endl;
        for (int i=0; i < size; i++) {
            result.push_back(invector[i] - sum);
            meanzero << " i " << result[i] << endl;
        }
        return result;
    }


// function to give r and theta relative to region centre
    pair <double,double> set_polars(int regNum){
        pair <double, double> result;
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

        result.first = xav - xcentre; //diff between seed point and barycentre
        result.second = yav - ycentre;
        return result;
    } //end of function set_polars

    void shift_polars (vector<double> hexVector, double angle) {
        for (auto i = hexVector.begin(); i != hexVector.end();++i)
            hexVector[*i] += angle;
    }

    //function to find all the edges and vertices of the internal boundary
    vector <std::pair<double,double>> dissectBoundary(void) {
        ofstream hfile ( "dissectdebug.out" );
        ofstream ifile ( "regionList.out" );
        ofstream kfile ( "edgesList.out" );
        hfile<<"just in dissectBoundary"<<endl;
        vector<std::pair<double,double>> result;
        std::pair<double,double> holdpair;
        vector<int> regionBoundary;
        int kcount = 0;
        int lcount = 0;
        bool exit_loops;
        for (int i=0; i < NUMPOINTS; i++) {
            hfile <<" kcount "<<kcount<<endl;
            exit_loops = false;
            regionBoundary.resize(0);

            for (unsigned int j = 0; j < this->regionIndex[i].size(); ++j) {
                // hfile<< "index j "<< j <<" hex angle = " <<H[5][regionIndex[i][j]]<<endl;

                if (this->outerBoundary[regionIndex[i][j]] == true){
                     hfile<<"boundary region = "<<i <<endl;
                     exit_loops = true;
                     hfile <<" lcount "<<lcount<<endl;
                     lcount++;
                     break;
                }

                if (this->Creg[this->regionIndex[i][j]] >0){
                    // hfile << "boundary i = "<< regionIndex[i][j] << " Creg = "<< this->Creg[this->regionIndex[i][j]]<<endl;
                    regionBoundary.push_back(this->regionIndex[i][j]);
                }
            } //end of loop on a single region

            if (exit_loops == true){ //any hexes outer boundary? cycle and zerosize regionBoundary
                 kcount++;
                 continue;
             }

            hfile<<"after the filling of regionEdge and regionVertex" <<endl;
            hfile<<"regionBoundary.size "<<regionBoundary.size() << endl;


            vector<double> rB; //contains the boundary thetas
            vector<int> irB; //holds the sorted boundary indicies
            for (unsigned int j = 0; j < regionBoundary.size();j++){
                //  hfile<< "index "<< j <<" vertex angle = " <<H[5][regionIndex[i][j]]<<endl;
                rB.push_back(H[5][regionBoundary[j]]);
            }
            irB = sort_indexes(rB); //indices after sort on theta

             // double angleOffset = -rB[irB[0]];
             // this->shift_polars(rB, angleOffset);
             // hfile<<"after shift_polars call"<<endl;

             unsigned int idissect = 0;
             int count = 0;
             int Vcount = 0;
             int Ecount = 0;
             unsigned int offset = 0;
             bool lvertex;
             vector<int> ihE; //contains the sorted indicies of each edge
             while (offset < irB.size()) {
                 if  (Creg[regionBoundary[irB[offset]]] > 1) {
                     Vcount++; //its a vertex
                     idissect++;
                      break;
                 }
                 offset++; //holds the number of edge indices before first vertex
              }
             hfile<<"after offset loop" << " offset " << offset << " idissect " << idissect << endl;
              unsigned int Size = irB.size();
              // if (Size == 0) {
              //     hfile << "region i " << region[i][0] << " irB size " << Size << endl;
              //     continue;
              // }
              while ((idissect < Size) && (count < 1000)) {
                 count++;
                 lvertex = false;
                 Ecount = 0;
                 ihE.resize(0);
                 while (Creg[regionBoundary[irB[(idissect+offset) % Size]]] > 1) {
                     lvertex = true;
                     idissect++;
                     continue;
                 }
                 hfile << "after vertex loop" <<endl;
                 if (lvertex)
                     Vcount++;

                 while (Creg[regionBoundary[irB[(idissect + offset)%Size]]] == 1) {
                     ihE.push_back(regionBoundary[irB[(idissect+offset)%Size]]);
                     Ecount++;
                     idissect++;
                 }
                 hfile << "after edge loop" <<endl;
                 // if (Ecount == 0){
                 //   hfile<<"***************************************************************"<<endl;
                 //continue;
                 //}
                  hfile<<"after filling tempvector"<<endl;
                  hfile<<"ihE Size "<<ihE.size()<<endl;
                  // hfile<<"ihE[1] "<<ihE[1]<<endl;
                  hfile<<"Vcount "<< Vcount << " Ecount "<< Ecount << " idissect = " << idissect << endl;
                  int regMiddle = region[ihE[0]][0];
                  int edgeOuter = -1;
                  for (int ihex = 0; ihex<6; ihex++){ //find the first region not the same as the central
                      if (regMiddle != hexRegionList[ihE[0]][ihex]){
                              edgeOuter = hexRegionList[ihE[0]][ihex];
                              break;
                      }
                  }
                   hfile<<"after edgeOuter assignment"<<endl;
                   hfile<<"edgeOuter "<<edgeOuter<< " edgeInner " << i << endl;
                   std::pair <int,int> keypair(i,edgeOuter);
                   int keyint = this->pair2int(keypair,this->base);
                   std::pair <int, vector<int>> p1(keyint,ihE);
                   hfile <<"after pair set"<<endl;
                   this->edges.insert(p1);
                   hfile << "after edges insert" << endl;
                   hfile <<"=================================================="<<endl;
                   this->regionList[i].push_back(edgeOuter);
             }

              for (unsigned int iregion = 0; iregion < regionList[i].size(); iregion++)
                  ifile << " r " << i << " rNbr " << regionList[i][iregion];
              ifile << endl;
              ifile << "---------------------------------------"<< endl;



        } //end of loop on regions
        int difference = 0;
        for (int irlook = 0; irlook < NUMPOINTS; irlook++)
            for (int jrlook = 0; jrlook < irlook; jrlook++) {
                std::pair <int,int> klook(irlook,jrlook);
                int k = pair2int(klook,this->base);
                if (edges.count(k) == 0)
                    continue;
                else {
                    int sizeij = edges[k].size();
                    std::pair <int, int> koolk(jrlook,irlook);
                    k = pair2int(koolk,this->base);
                    int sizeji = edges[k].size();
                    difference = (sizeij - sizeji)*(sizeij - sizeji);
                    if ((sizeij*sizeji != 0) && (difference < 1000))
                        kfile << irlook << " " << jrlook << " sizeij " << sizeij << " "<< jrlook << " " << irlook << " sizeji " << sizeji << endl;
                }
            }

        return result;
    }//end of function dissectBoundary

    // function to correlate matching edges
     void correlate_edges()  {
        ofstream edgefile("logs/edgeCorrelations.txt");
        vector<double> tempvect1;
        vector<double> tempvect2;
        vector<double> tempvect3;
        //debug edgefile << " In correlate_edges " << endl;
        // iterate over regions
        //for each region iterate over region edges
        // for each edge pair (i,j), (j,i) call correlate_vectors
        // write the i,j and correlation to a file
        // what happens if there are multiple entries in regionList at the start and the end?
        // there are some regions that have this
        for (int i = 0; i <NUMPOINTS; i++) {
            //debug edgefile << " i iteration " << i << endl;
            for (auto j = this->regionList[i].begin(); j != this->regionList[i].end();j++) {
                //debug edgefile << " j iteration " << *j << endl;
                tempvect1.resize(0);
                tempvect2.resize(0);
                tempvect3.resize(0);
                std::pair<int,int> edgePair1(i,*j);
                int edgeIndex1 = this->pair2int(edgePair1,this->base);
                std::pair<int,int>  edgePair2(*j,i);
                int edgeIndex2 = this->pair2int(edgePair2,this->base);
                int count1 = 0;
                int count2 = 0;
                for (auto itr = this->edges[edgeIndex1].begin(); itr != this->edges[edgeIndex1].end();itr++){
                    tempvect1.push_back(this->NN[*itr] - 1.0);
                    count1++;
                }
                // edgefile << " after filling tempvector1 " << endl;
                for (auto itr = this->edges[edgeIndex2].begin(); itr != this->edges[edgeIndex2].end();itr++) {
                    tempvect2.push_back(this->NN[*itr] - 1.0);
                    count2++;
                }
                   std::reverse(tempvect2.begin(),tempvect2.end());
                // edgefile << " after filling tempvector2 " << endl;
                //int correlationValue = this->correlate_vector(tempvect1,tempvect2);
                //debug edgefile << i << " tv1 " << tempvect1.size() << " j " << *j << " tv22  " <<  tempvect2.size() << endl;
                if (tempvect1.size() == tempvect2.size() && tempvect1.size()*tempvect2.size() != 0){
                    double correlationValue = this->correlate_Eqvector(tempvect1,tempvect2);
                    edgefile << " i " <<  i << " c1 " << count1 << " j " << *j << " c2  " <<  count2 << " cV " << correlationValue << endl;
                    //edgefile << i << " Size1 " << edges[edgeIndex1].size() << " j " << *j << " Size2  " << edges[edgeIndex2].size() << endl;
                } //end of code if both edges are equal
                else if (tempvect1.size() > tempvect2.size() && tempvect1.size()*tempvect2.size() != 0) {
                    tempvect3 = this->equalize_vector(tempvect2,tempvect1);
                    //debug if (tempvect3.size() == 0)
                    //debug edgefile << " i " << i << " count1 " << tempvect1.size() << " j " << *j << " tempvect3 is zero  "   << endl;
                    double correlationValue = this->correlate_Eqvector(tempvect1,tempvect3);
                    edgefile << " i " << i << " c1 " << count1 << " j " << *j << " c2  " <<  count2 << " cV " << correlationValue << endl;
                }
                else if (tempvect1.size() < tempvect2.size() && tempvect1.size()*tempvect2.size() != 0) {
                    tempvect3 = this->equalize_vector(tempvect1,tempvect2);
                    //debug if (tempvect3.size() == 0)
                    //debug edgefile << " i " << i << " count2 " << tempvect2.size() << " j " << *j << " tempvect3 is zero  "   << endl;
                    double correlationValue = this->correlate_Eqvector(tempvect3,tempvect2);
                    edgefile << " i = " << i << " c1 " << count1 << " j " << *j << " c2  " <<  count2 << " cV " << correlationValue << endl;
                }
            } //end of single edge comparison
        } // end of loop on regions
     } //end of function correlate_edges

    //function to return the correlation of two vectors
    double correlate_Eqvector(vector<double> vector1, vector<double> vector2) {
        ofstream jfile;
        jfile.open("correlateEqvector.out",ios::app);
        double result;
        jfile << " In correlateEqvector " << endl;
        if (vector1.size() != vector2.size()){
            jfile << "error: vectors must be same length" << endl;
            return -2;
        }
        double vector1Norm = 0;
        double vector2Norm = 0;
        double vector12Product = 0;
        unsigned int vectSize = vector1.size();
        // for (auto i = vector1.begin(); i != vector1.end(); i++) {
        for (unsigned int  i = 0; i != vectSize; i++) {
            vector1Norm +=  vector1[i]*vector1[i];
            vector2Norm += vector2[i]*vector2[i];
            vector12Product += vector1[i]*vector2[i];
        }
            vector1Norm = sqrt(vector1Norm);
            vector2Norm = sqrt(vector2Norm);
            result = vector12Product / (vector1Norm*vector2Norm);
            jfile << "result = " << result << endl;
            return result;

    } //end of function correlateEqvector

    vector <double> equalize_vector (vector<double> svector, vector<double> lvector) {
        ofstream kfile;
        kfile.open("eqVector.out",ios::app);
        vector <double> result;
        int sSize, lSize = 0;
        double sStep, lStep = 0;
        // double delta = 0.000001;
        result.resize(0);
        sSize = svector.size();
        lSize = lvector.size();
        sStep = 1.0 / (1.0 * (sSize-1.0));
        lStep = 1.0 / (1.0 * (lSize-1.0));
        double start = 0;
        double finish = 0;
        double value = 0;
        int count = 1;
        result.push_back(svector[0]);
        for (int i = 0; i < sSize-1; i++) { // walk along the short vector
            start = i*sStep;
            finish = (i+1)*sStep;
            for (int j = 0; j < lSize; j++) { //walk along the long
                if ((j*lStep > start) && (j*lStep) <= finish) {
                    value = svector[i]*(j*lStep - start) + svector[i+1]*((j+1)*(finish - j*lStep));
                    result.push_back(value);
                    count++;
                }
            }
        }
        if (count != lSize){
            kfile << " count " << count << " lSize " << lSize << " sSize " << sSize <<  " not filled" << endl;
            result.resize(0);
            return result;
        }
        else {
            kfile << " count " << count << " l size " << lSize << " sSize " << sSize << " filled" << endl;
            // result.push_back(svector[sSize - 1]);
            return result;
        }
    } //end of function equalize_vector

    double min_radius(int regNum) {
        point  barycentre;
        point boundHex;
        barycentre.xval = diff[regNum].first + centres[regNum].xval;
        barycentre.yval = diff[regNum].second + centres[regNum].yval;
        double minRadius = 100000.0;
        for (unsigned int i = 0; i < regionIndex[regNum].size();i++) {
            if (Creg[i] > 0) {
                boundHex.xval = H[0][regionIndex[regNum][i]],
                boundHex.yval = H[1][regionIndex[regNum][i]];
                double boundDist = getdist(boundHex,barycentre);

                if (boundDist < minRadius)
                    minRadius = boundDist;
            }
        }
        return minRadius;
    }



    //sectorize over radius
    vector <double> sectorize_reg_radius (int regNum, int numSectors, int beginAngle, int endAngle) {
    ofstream dfile ( "logs/sectorRadius.txt" );
    vector <double>  radiusCC;
    radiusCC.resize(numSectors);
    vector <int> radiusCount;
    radiusCount.resize(numSectors);
    double startRadius, endRadius, radiusInc; //sector radii
    double maxRadius = min_radius(regNum);
    radiusInc = maxRadius /(1.*numSectors);
    double startAngle, finishAngle, angleInc; //sector angles
    angleInc = 2*PI/(1.*numSectors);
    startAngle = beginAngle*angleInc;
    finishAngle = endAngle*angleInc;
    int size = (int) regionIndex[regNum].size();

    for (int k=0;k<numSectors;k++) {
      startRadius = (k*radiusInc);
      endRadius = (k+1)*radiusInc;

      for (int i=0; i< size;i++) {
	if (this->H[5][regionIndex[regNum][i]]>=startAngle && this->H[5][regionIndex[regNum][i]]<finishAngle) {
	  if (this->H[4][regionIndex[regNum][i]]>=startRadius && this->H[4][regionIndex[regNum][i]]<endRadius) {
	      radiusCount[k]++;
	      radiusCC[k] += this->CC[regionIndex[regNum][i]];
	  } //end of if on radius
	} //end of if on angleSector
      } //end of loop over i

      dfile <<" radiusInc "<<radiusInc <<"  startRadius "<<startRadius<<"  endRadius "<<endRadius<<endl;
      if (radiusCount[k] != 0)
	radiusCC[k] = radiusCC[k] / (1.*radiusCount[k]) - 2.5;
      else
	radiusCC[k] += -999.999;
    }//end loop on k
     return radiusCC;

  } //end of function sectorize_radius

 //function to count the hexes in sectors of a region via angular sectors
    vector <double> sectorize_reg_angle (int regNum, int numSectors, int beginRadius, int endRadius) {
    ofstream cfile ("logs/sectorAngle.txt");
    //std::pair<double,double> diff; //difference between seed point and CoG of region
    // diff = this->set_polars(regNum);
    vector <double> angleCC; //average value of CC in each sector
    vector <double>  angleVal;
    vector <int> count; //number of hexes in each sector
    angleCC.resize(numSectors);
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
	  angleCC[k] += this->CC[regionIndex[regNum][i]];
	}
	//cfile << setw(5) << angleVal[i]  <<"  ";
      }

      cfile << endl;
      // cfile << "diff x " << diff.first << "diff y " << diff.second << endl;
      cfile << " angleinc " <<angleInc << "  start angle "  << startAngle << "  end angle "<< endAngle <<endl;
      if (count[k] != 0)
	angleCC[k] = angleCC[k] / (1.*count[k]) - 2.5;
      else
	angleCC[k] = -999.999;
    }//end loop on k
     return angleCC;

  } //end of function sectorize_region


}; // Erm2009


int main (int argc, char **argv)
{
    srand(atoi(argv[3]));
    ofstream bfile ( "logs/sector.out" );
    ofstream gfile ( "logs/edges.out");
    double dt = 0.00005;            // integration timestep // SHOULD BE 0.001

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
       double choice = morph::Tools::randDouble();
       if (choice > 0.5)
           M.NN[i]=-(morph::Tools::randDouble())*1.0 + 1.;
       else
           M.NN[i]=(morph::Tools::randDouble())*1.0 + 1.;
       M.CC[i]=(morph::Tools::randDouble())*1.0 + 2.5;
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
        // const char* TIMEcs = TIMEss.str().c_str();

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

            vector<std::pair<double,double>> centroids;
            centroids = M.dissectBoundary();
            gfile<<"before correlate_edges" << endl;
            M.correlate_edges();
            gfile<<"after correlate_edges" << endl;
            cout<<"after correlate_edges" << endl;
            int regionCount = 0;
            int numSectors = 12;
	    vector <double> angleVector;
            vector <double> radiusVector;
            int degreeRadius;
            int degreeAngle;
	    double tempArea;
	    double tempPerimeter;
            int angleOffset = 0;
            int radiusOffset = 7;
              for (int j=0;j<NUMPOINTS-1;j++) {
	      if (M.regArea(j) != 0){
                  //gfile<<"in the degree loop" << endl;
                  angleVector = M.sectorize_reg_angle(j,numSectors,radiusOffset, radiusOffset + 4);
                  angleVector = M.meanzero_vector(angleVector);
                  radiusVector = M.sectorize_reg_radius(j,numSectors, angleOffset, angleOffset + 5);
                  // radiusVector = M.meanzero_vector(radiusVector);
                  // gfile << "tempvector size " << tempvector.size() << endl;
                  //for (unsigned int i = 0; i<tempvector.size();i++)
                  //   gfile << "i = 0 "<< i << " tempvector = " << tempvector[i] << endl;
                  tempArea = M.regArea(j)*(5.0/Dn);
                  tempPerimeter = M.regPerimeter(j)*sqrt(5.0/Dn);
                  //gfile << "just before degree" <<endl;
                  //  degreeTurn = M.find_max(tempvector);
                  degreeRadius = M.find_max(radiusVector);
                  degreeAngle = M.find_zero(angleVector);
                  W.logfile <<" degreeRadius "<< degreeRadius<<" degreeAngle "<< degreeAngle << " " << tempArea<<"  "<<tempPerimeter<<endl<<flush;
                  gfile << "region "<< j << " degreeRadius "<< degreeRadius<<" degreeAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;
                  regionCount++;
                  bfile <<"region "<<j<<" ";
                  for (int k=0; k <= (int) radiusVector.size();k++){
                      bfile<< setw(5) << " "<< radiusVector[k%radiusVector.size()];
                      bfile <<endl;
                  } //end of for on printing sector values
	      } //end of if on non-zero regions
              bfile <<regionCount<<endl;
	      } //end of loop on NUMPOINTs



            break;
        }

        case 6: // *** LOAD ***
        {
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
