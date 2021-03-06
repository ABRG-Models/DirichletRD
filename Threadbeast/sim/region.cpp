#include <morph/tools.h>
#include <morph/HexGrid.h>
#include <morph/ReadCurves.h>
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
// #include <boost/math/special_functions/bessel.hpp>
#define PI 3.1415926535897932
#define NUMPOINTS 25 //just the A-E rows.

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
using namespace std;

#define NE(hi) (this->Hgrid->d_ne[hi])
#define HAS_NE(hi) (this->Hgrid->d_ne[hi] == -1 ? false : true)

#define NW(hi) (this->Hgrid->d_nw[hi])
#define HAS_NW(hi) (this->Hgrid->d_nw[hi] == -1 ? false : true)

#define NNE(hi) (this->Hgrid->d_nne[hi])
#define HAS_NNE(hi) (this->Hgrid->d_nne[hi] == -1 ? false : true)

#define NNW(hi) (this->Hgrid->d_nnw[hi])
#define HAS_NNW(hi) (this->Hgrid->d_nnw[hi] == -1 ? false : true)

#define NSE(hi) (this->Hgrid->d_nse[hi])
#define HAS_NSE(hi) (this->Hgrid->d_nse[hi] == -1 ? false : true)

#define NSW(hi) (this->Hgrid->d_nsw[hi])
#define HAS_NSW(hi) (this->Hgrid->d_nsw[hi] == -1 ? false : true)
  
class Region
{
public:
    int scale;
    int n;
    double ds;
  struct point {
    double xval;
    double yval;
  };

  // list of objects visible to member functions
  vector<vector<int> > N; // hex neighbourhood 
  vector<vector<int> > region; //indices of all the regions
  vector<vector<int>> hexRegionList; //vector of the different regions around a hex
  vector<vector<int>> regionList; // vector of the regions around a region
  vector<vector<int>> regionVertex; //vectoir of vertices around a region
  std::map<int,vector<int>> edges; //map of (i,j) edges, uses pair<->int converters
    // std::map<int,double> > corrAdjacent;

  pair <float, float> centres[NUMPOINTS]; //seed points
  vector<vector<double> > regionDist; //distances to each seed point
  vector<vector<int> > regionIndex; //indexes of hexes in each region
  vector<int> C; //count of neighbours
  vector<int> Creg; //count of regions

  vector<double> NN, CC; //hold the field values for each hex
  vector<std::pair<double,double>> diff; //difference between seed point and CoG of region
  morph::HexGrid* Hgrid;
  vector<morph::Hex> *H;
  int base = 1000;
  //class constructor
  Region (int scale, string logpath) {
//        ofstream //afile (logpath + "/debug.out" );
//        ofstream jfile (logpath + "/correlateVector.out");

	srand(time(NULL));
    this->scale = scale;

	//centres.resize(NUMPOINTS);

    double s = pow(2.0, scale-1);
	ds = 1.0/s;
	//02-08-18 keep ds fixed now as region changes size.
	//ds = 0.453125;
	//int inring = s;

       	// std::default_random_engine generator;
        // std::uniform_int_distribution<int> distribution(0,inring);
	// generator.min(0);
	// generator.max(inring);
	//int red = 0; int blue = 0;
     // for (int i=0; i<NUMPOINTS; i++) {
     //     centres[i].xval = 0.88*(rand()/(RAND_MAX +1.0) - 0.5);
     //     centres[i].yval = 1.34*(rand()/(RAND_MAX + 1.0) - 0.5);
     // }


     ///////////////
#include "centresAE.h"

	
//	for (int j=0;j<NUMPOINTS;j++)
//          afile << "j = " << j <<" x = " << centres[j] .xval  << "  y = " << centres[j].yval <<endl;
//	afile << " s = " << s << endl;

  n = 0;
  Hgrid = new HexGrid(this->ds, 0.45, 0.0, morph::HexDomainShape::Boundary);
  morph::ReadCurves r("./barrelAE.svg");
  Hgrid->setBoundary (r.getCorticalPath());
  cout << "before filling H " << Hgrid->num() << endl;
  H->resize(Hgrid->num());
  for (auto h : Hgrid->hexen) {
    cout << " in H fill loop" << endl;
    H->push_back(h);
  }	
  cout << "after  filling H " << endl;
  n = H->size();
//these are the vectors of vectors for the regions
  regionDist.resize(n);
  region.resize(n);  
  C.resize(n,0); //count of neighbouring hexes
  N.resize(n);
  Creg.resize(n,0);
  hexRegionList.resize(n); //neighbouring regions for a hex
  regionList.resize(NUMPOINTS); //neighbouring regions for a region
  regionVertex.resize(NUMPOINTS); //vertices for a region
    
// making a neighbour array for convenience
  for (unsigned int idx = 0; idx < H->size(); idx++) {
    N[idx].resize(6);
    N[idx][0] = Hgrid->d_ne[idx];
    N[idx][1] = Hgrid->d_nne[idx];
    N[idx][2] = Hgrid->d_nnw[idx];
    N[idx][3] = Hgrid->d_nw[idx];
    N[idx][4] = Hgrid->d_nsw[idx];
    N[idx][5] = Hgrid->d_nse[idx];
  }

	//arrays to hold the field values
    NN.resize(n);
    CC.resize(n);

	// get a list of distances from each Dirichlet point for each hex.
    for (auto h : *H) {
      for (int j=0;j<NUMPOINTS;j++) {
        double temp = h.distanceFrom(centres[j]);
        // afile << "   " <<hexcen.xval <<"  "<<hexcen.yval;
        regionList[h.vi].push_back(temp);
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

     //afile <<"list of all the regions"<<endl;
     for(int i=0;i<n;i++){

       if (sortedDist[i][0] == sortedDist[i][1]){
           //afile<<"i ="<<i;
     	//   afile << "  "<<region[i][0]<<" is equidistant from "<< region[i][1];
       }
     }


     // this determines the type of hex
      for (auto h : *H){

          if(h.boundaryHex == true) {
            if (!HAS_NE(h.vi)) {
				hexRegionList[h.vi][0] = -1;
            } 
            if (!HAS_NNE(h.vi)) {
				hexRegionList[h.vi][1] = -1;
            } 
            if (!HAS_NNW(h.vi)) {
				hexRegionList[h.vi][2] = -1;
            } 
            if (!HAS_NW(h.vi)) {
				hexRegionList[h.vi][3] = -1;
            } 
            if (!HAS_NSW(h.vi)) {
				hexRegionList[h.vi][4] = -1;
            } 
            if (!HAS_NSE(h.vi)) {
				hexRegionList[h.vi][5] = -1;
            } 
		  } //end of if dealing with outer boundary

		  // processing of internal boundaries
          int centralRegion = region[h.vi][0];
          int oldRegion = centralRegion;
          int newRegion = 0;

          for(int j=0;j<6;j++) {
              hexRegionList[h.vi].push_back(region[N[h.vi][j]][0]); //push back the region of each neighbour
              newRegion =  region[N[h.vi][j]][0];
              if (centralRegion != newRegion){ //its a boundary hex
                  C[h.vi]--;
              //    N[i][j] = h.vi; pre-morphologica 
                  if (oldRegion != newRegion){ //logic to test if vertex
                    oldRegion = newRegion;
                    Creg[h.vi]++;
                    }
                  }
              }
      } //end of logic determining hex types



      //afile << "after creating internal boundaries" << endl;




    //information about the indexes of hexes in regions
     regionIndex.resize(NUMPOINTS);
      for (int j=0;j<NUMPOINTS;j++) {
        for (auto h : *H){
	      if (region[h.vi][0] == j)
      	   regionIndex[j].push_back(h.vi);
        }
      }


      diff.resize(NUMPOINTS);
      for (int j=0;j<NUMPOINTS;j++){
          diff[j] = this->set_polars(j);
          //afile << "diff seed-centre" << diff[j].first << " " << diff[j].second<<endl;
      }
  } //end of constructor

  //function, gets distance between points now supseded by morphologica function
  //float distanceFrom (const pair<float, float> cartesianPoint) const {

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

    vector<double> getLaplacian(vector<double> Q, double dx) {
      double overdxSquare = 1./(1.5*dx*dx);
        vector<double> L(n,0.);
        for(auto h : *this->H){
         int i = int(h.vi);
            L[i]=(Q[N[i][0]]+Q[N[i][1]]+Q[N[i][2]]+Q[N[i][3]]+Q[N[i][4]]+Q[N[i][5]]-6.*Q[i])*overdxSquare;
        }
        return L;
    }

  //function to timestep coupled equations
    void step(double dt, double Dn, double Dchi, double Dc) {
      dt = dt * 2.5 / Dn;


       // Set up boundary conditions with ghost points

      for(auto h : *this->H){
        if(C[h.vi]<6){
          for(int j=0;j<6;j++){
		     int i = int(h.vi);
             if(N[h.vi][j] == i) {
       	       NN[N[h.vi][j]] = NN[h.vi];
	           CC[N[h.vi][j]] = CC[h.vi];
	      }
	     }
	   }
      }


        double beta = 5.;
        vector<double> lapN = getLaplacian(NN,ds);
        vector<double> lapC = getLaplacian(CC,ds);
        vector<double> G(n,0.);
        double a = 1., b = 1., mu = 1;

        for (auto h : *H) {
		  unsigned int i = h.vi;
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
        for (auto h : *H) {
          NN[h.vi]+=dt*( a-b*NN[h.vi] + Dn*lapN[h.vi] - Dchi*G[h.vi]);
        }

        // step C
        double N2;
        for(auto h : *H){
		    unsigned int i = h.vi;
            N2 = NN[i]*NN[i];
            CC[i]+=dt*( beta*N2/(1.+N2) - mu*CC[i] + Dc*lapC[i] );
        }
    }//end step

  //function to return area of a region
  double regArea (int regNum) {
    double area = 0;
    for (int i=0;i < (int) regionIndex[regNum].size();i++){
	  area += 1.;
	}
    return area;
  } //end of funtion regArea

  //function to return perimeter of a region
  double regPerimeter (int regNum) {
    // cout << "in regPerimenter " << endl;
    double perimeter = 0;
    for (int i=0;i < (int) regionIndex[regNum].size();i++)
      if (Creg[regionIndex[regNum][i]] > 0)
	perimeter += 1.0;
    return perimeter;
  } //end of function regPerimeter

    // transform vector so its mean is zero
    vector<double> meanzero_vector(vector<double> invector) {
        //ofstream meanzero ("meanzero.txt",ios::app);
        vector <double> result;
        int size = invector.size();
        //meanzero << "size " << size << endl;
        double sum = 0;
        double absSum = 0;
        for (int i=0; i <size; i++) {
            sum += invector[i];
            absSum += fabs(invector[i]);
        }
        sum = sum/(1.0*size);
        absSum = absSum / (1.0 * size);
        //meanzero << " mean  " << sum << endl;
        //meanzero << " absolute mean  " << absSum << endl;
        for (int i=0; i < size; i++) {
            result.push_back(invector[i] - sum);
            //meanzero << " i " << result[i] << endl;
        }
        return result;
    }

    // return the mean of the absolute values of a  vector
        double absmean_vector(vector<double> invector) {
        //ofstream meanzero ("meanzero.txt",ios::app);
          double result = 0;
	      double sum = 0;
          int size = invector.size();
        //meanzero << "size " << size << endl;
          for (int i=0; i <size; i++) {
            sum += fabs(invector[i]);
        }
        sum = sum/(1.0*size);
        return result;
    } 



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
  //function to return fraction of area with NN positive
  double regNNfrac (int regNum) {
    double area = 0;
    double positive_area = 0;
    int size = regionIndex[regNum].size();
    vector<double> normalNN;
    // to normalise the NN field
     for (int i=0;i<size;i++){
          normalNN.push_back(this->NN[regionIndex[regNum][i]]);
      }
      normalNN = this->meanzero_vector(normalNN);
      for (int i=0;i < (int) regionIndex[regNum].size();i++){
	  area += 1.;
	  if ((normalNN[i]) > 0) {
		  positive_area += 1.0;
	  }
    }
    return positive_area / area;
  } //end of function regNNfrac

// function to give r and theta relative to region centre
    pair <double,double> set_polars(int regNum){
        pair <double, double> result;
        double xcentre = this->centres[regNum].first;
        double ycentre = this->centres[regNum].second;
        double xav=0;
        double yav = 0;
        int hexcount = 0;
        for (int i=0;i< (int) this->regionIndex[regNum].size();i++) {
            hexcount++;
            xav += this->H->at(regionIndex[regNum][i]).x;
            yav += this->H->at(regionIndex[regNum][i]).y;
        }
        xav = xav / hexcount;
        yav = yav / hexcount;

//go over the region and put the hexes into bins then average
        for (int i=0;i< (int) this->regionIndex[regNum].size();i++) {
            int index = regionIndex[regNum][i];
            this->H->at(index).r = sqrt((this->H->at(index).x-xav)*(this->H->at(index).x-xav) + \
                                     (this->H->at(index).y-yav)*(this->H->at(index).y-yav));
            if (this->H->at(i).x >= 0)
                this->H->at(index).phi = atan2((this->H->at(index).y-yav), (this->H->at(index).x-xav));
            else
                this->H->at(index).phi = PI + atan2((this->H->at(index).y-yav), (this->H->at(index).x-xav));
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
    vector <std::pair<double,double>> dissectBoundary(string logpath) {
        ofstream hfile ( logpath + "/dissectdebug.out" );
        ofstream ifile ( logpath + "/regionList.out" );
        ofstream kfile ( logpath + "/edgesList.out" );
        ofstream lfile ( logpath + "/vertices.List.out" );
        hfile<<"just in dissectBoundary"<<endl;
        vector<std::pair<double,double>> result;
        // std::pair<double,double> holdpair;
        vector<int> regionBoundary;
        for (int i=0; i < NUMPOINTS; i++) {
          regionBoundary.resize(0);
          // fill the regionBoundary
          for (unsigned int j = 0; j < this->regionIndex[i].size(); ++j) {

            if (this->Creg[this->regionIndex[i][j]] >0){
                    // hfile << "boundary i = "<< regionIndex[i][j] << " Creg = "<< this->Creg[this->regionIndex[i][j]]<<endl;
              regionBoundary.push_back(this->regionIndex[i][j]);
            }
          } //end of loop on a single region


         hfile<<"after the filling of regionEdge and regionVertex" <<endl;
         hfile<<"regionBoundary.size "<<regionBoundary.size() << endl;


         vector<double> rB; //contains the boundary thetas
         vector<int> irB; //holds the sorted boundary indicies
         for (unsigned int j = 0; j < regionBoundary.size();j++){
                //  hfile<< "index "<< j <<" vertex angle = " <<H[5][regionIndex[i][j]]<<endl;
           rB.push_back(H->at(regionBoundary[j]).phi);
         }
         irB = sort_indexes(rB); //indices after sort on theta

         // double angleOffset = -rB[irB[0]];
         // this->shift_polars(rB, angleOffset);
         // hfile<<"after shift_polars call"<<endl;
         // we now walk round the boundary in theta order
         unsigned int count = 0;
         int Vcount = 0; //count of the vertices
         int Ecount = 0; //count of the edges
         unsigned int offset = 0; //number of boundary hexes before the first vertex
         bool lvertex;
         vector<int> ihE; //contains the sorted indicies of each edge
         while (offset < irB.size()) {
               if  (Creg[regionBoundary[irB[offset]]] > 1) {
                 Vcount++; //its a vertex
				 regionVertex[i].push_back(irB[offset]);
                 break;
                }
                offset++; //holds the number of edge indices before first vertex
         }
         hfile<<"after offset loop" << " offset " << offset  << endl;
         unsigned int Size = irB.size();
              // if (Size == 0) {
              //     hfile << "region i " << region[i][0] << " irB size " << Size << endl;
              //     continue;
              // }
         while (count < Size+1) {
              count++;
              lvertex = false;
              Ecount = 0;
              ihE.resize(0);
              while (Creg[regionBoundary[irB[(count + offset) % Size]]] > 1) {
                   lvertex = true;
                   continue;
              }
              hfile << "after vertex loop" <<endl;
              if (lvertex) {
                Vcount++;
				regionVertex[i].push_back(irB[(count+offset)%Size]);
				}
              //walk along the edge until the next vertex
              while (Creg[regionBoundary[irB[(count + offset)%Size]]] == 1) {
                   ihE.push_back(regionBoundary[irB[(count+offset)%Size]]);
                   Ecount++;
                 }
                 hfile << "after edge loop" <<endl;
                 // if (Ecount == 0){
                 //   hfile<<"***************************************************************"<<endl;
                 //continue;
                 //}
               hfile<<"after filling tempvector"<<endl;
               hfile<<"ihE Size "<<ihE.size()<<endl;
			   for (unsigned int ii = 0; ii<ihE.size();ii++){
               hfile<<ihE[ii]<<endl;
			   }
               hfile<<"Vcount "<< Vcount << " Ecount "<< Ecount << " vertices counted = " << count << endl;
               int regMiddle = region[ihE[0]][0];
               int edgeOuter = -2;
               for (int ihex = 0; ihex<6; ihex++){ //find the first region not the same as the central, since Creg = 1 there can only be one such region.
                  if (regMiddle != hexRegionList[ihE[0]][ihex]){
                     edgeOuter = hexRegionList[ihE[0]][ihex];
                     break;
                  }
                }
                hfile<<"after edgeOuter assignment"<<endl;
                hfile<<"edgeOuter "<<edgeOuter<< " edgeInner " << i << endl;
				if (edgeOuter > -1) {
                  std::pair <int,int> keypair(i,edgeOuter);
                  int keyint = this->pair2int(keypair,this->base);
                  std::pair <int, vector<int>> p1(keyint,ihE);
                  hfile <<"after pair set"<<endl;
                  this->edges.insert(p1);
                  hfile << "after edges insert" << endl;
                  hfile <<"=================================================="<<endl;
                  this->regionList[i].push_back(edgeOuter);
				}
             }

             for (unsigned int iregion = 0; iregion < regionList[i].size(); iregion++){
                ifile << " r " << i << " rNbr " << regionList[i][iregion];
                ifile << endl;
                ifile << "---------------------------------------"<< endl;
             }

             for (unsigned int iregion = 0; iregion < regionVertex[i].size(); iregion++){
                lfile << " r " << i << " vertex " << regionList[i][iregion];
                lfile << endl;
                lfile << "---------------------------------------"<< endl;
             }

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
     double correlate_edges(string logpath)  {
	double result = 0;     
        ofstream edgefile(logpath + "/edgeCorrelations.txt");
        vector<double> tempvect1;
        vector<double> tempvect2;
        vector<double> tempvect3;
        edgefile << " In correlate_edges " << endl;
        // iterate over regions
        //for each region iterate over region edges
        // for each edge pair (i,j), (j,i) call correlate_vectors
        // write the i,j and correlation to a file
        // what happens if there are multiple entries in regionList at the start and the end?
        // there are some regions that have this
        int countResult = 0;
        for (int i = 0; i <NUMPOINTS; i++) {
		vector<double> normalNN;
		double NNmean;
		int size = (int) regionIndex[i].size();
		normalNN.resize(size,0.0);
		// create a vector of the NN values in the region
		for (int j=0;j<size;j++){
			//edgefile << "create normalNN j " << j <<endl;
		       	normalNN.push_back(this->NN[regionIndex[i][j]]);
		}
	       	NNmean = this->absmean_vector(normalNN);
		edgefile << " mean of NN in region " << i << "is " <<NNmean << endl;
                edgefile << " i iteration " << i << endl;
            for (auto j = this->regionList[i].begin(); j != this->regionList[i].end();j++) {
                //edgefile << " j iteration " << *j << endl;
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
                    tempvect1.push_back(this->NN[*itr] - NNmean);
                    count1++;
                }
                // edgefile << " after filling tempvector1 " << endl;
                for (auto itr = this->edges[edgeIndex2].begin(); itr != this->edges[edgeIndex2].end();itr++) {
                    tempvect2.push_back(this->NN[*itr] - NNmean);
                    count2++;
                }
                   std::reverse(tempvect2.begin(),tempvect2.end());
                // edgefile << " after filling tempvector2 " << endl;
                //int correlationValue = this->correlate_vector(tempvect1,tempvect2);
                //edgefile << i << " tv1 " << tempvect1.size() << " j " << *j << " tv22  " <<  tempvect2.size() << endl;
                if (tempvect1.size() == tempvect2.size() && tempvect1.size()*tempvect2.size() != 0){
                    double correlationValue = this->correlate_Eqvector(tempvect1,tempvect2);
		    result += fabs(correlationValue);
                    countResult++;
                    edgefile << " i = " << i << " c1 " << count1 << " j " << *j << " c2  " <<  count2 << " cV " << correlationValue << endl;
                    //edgefile << i << " Size1 " << edges[edgeIndex1].size() << " j " << *j << " Size2  " << edges[edgeIndex2].size() << endl;
                } //end of code if both edges are equal
                else if (tempvect1.size() > tempvect2.size() && tempvect1.size()*tempvect2.size() != 0) {
                    tempvect3 = this->equalize_vector(tempvect2,tempvect1);
                    if (tempvect3.size() == 0)
                        edgefile << " i " << i << " count1 " << tempvect1.size() << " j " << *j << " tempvect3 is zero  "   << endl;
                    double correlationValue = this->correlate_Eqvector(tempvect1,tempvect3);
		    result += fabs(correlationValue);
                    countResult++;
                    edgefile << " i " << i << " c1 " << count1 << " j " << *j << " c2  " <<  count2 << " cV " << correlationValue << endl;
                }
                else if (tempvect1.size() < tempvect2.size() && tempvect1.size()*tempvect2.size() != 0) {
                    tempvect3 = this->equalize_vector(tempvect1,tempvect2);
                    if (tempvect3.size() == 0)
                        edgefile << " i " << i << " count2 " << tempvect2.size() << " j " << *j << " tempvect3 is zero  "   << endl;
                    double correlationValue = this->correlate_Eqvector(tempvect3,tempvect2);
		    result += fabs(correlationValue);
                    countResult++;
                    edgefile << " i = " << i << " c1 " << count1 << " j " << *j << " c2  " <<  count2 << " cV " << correlationValue << endl;
                }
            } //end of single edge comparison
        } // end of loop on regions
	edgefile << " countResult "<<countResult<<endl;
	result = result / (countResult * 1.0);
	return result;
     } //end of function correlate_edges

    //function to return the correlation of two vectors
    double correlate_Eqvector(vector<double> vector1, vector<double> vector2) {
        //ofstream jfile;
        //jfile.open("correlateEqvector.out",ios::app);
        double result;
        //jfile << " In correlateEqvector " << endl;
        if (vector1.size() != vector2.size()){
            //jfile << "error: vectors must be same length" << endl;
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
           // jfile << "result = " << result << endl;
            return result;

    } //end of function correlateEqvector

    vector <double> equalize_vector (vector<double> svector, vector<double> lvector) {
        //ofstream kfile;
        //kfile.open("eqVector.out",ios::app);
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
          //  kfile << " count " << count << " lSize " << lSize << " sSize " << sSize <<  " not filled" << endl;
            result.resize(0);
            return result;
        }
        else {
            // kfile << " count " << count << " l size " << lSize << " sSize " << sSize << " filled" << endl;
            // result.push_back(svector[sSize - 1]);
            return result;
        }
    } //end of function equalize_vector


    double min_radius(int regNum) {
        point  barycentre;
        point boundHex;
        barycentre.xval = this->diff[regNum].first + centres[regNum].first;
        barycentre.yval = this->diff[regNum].second + centres[regNum].second;
        double minRadius = 100000.0;
        for (unsigned int i = 0; i < regionIndex[regNum].size();i++) {
            if (Creg[regionIndex[regNum][i]] > 0) {
                boundHex.xval = H->at(regionIndex[regNum][i]).x,
                boundHex.yval = H->at(regionIndex[regNum][i]).y;
                double boundDist = getdist(boundHex,barycentre);

                if (boundDist < minRadius)
                    minRadius = boundDist;
            }
        }
        return minRadius;
    }

     double max_radius(int regNum) {
        point  barycentre;
        point boundHex;
        barycentre.xval = this->diff[regNum].first + centres[regNum].first;
        barycentre.yval = this->diff[regNum].second + centres[regNum].second;
        double maxRadius = -100000.0;
        for (unsigned int i = 0; i < regionIndex[regNum].size();i++) {
                boundHex.xval = H->at(regionIndex[regNum][i]).x,
                boundHex.yval = H->at(regionIndex[regNum][i]).y;
                double boundDist = getdist(boundHex,barycentre);

                if (boundDist > maxRadius)
                    maxRadius = boundDist;
        }
        return maxRadius;
    }


    //sectorize over radius
    vector <double> sectorize_reg_radius (int regNum, int numSectors, int beginAngle, int endAngle) {
    ofstream dfile ("logs/sectorRadius.txt",ios::app);
    vector <double>  radiusNN;
    radiusNN.resize(numSectors,0);
    vector <int> radiusCount;
    vector <double> normalNN;
    radiusCount.resize(numSectors,0);
    double startRadius, finishRadius, radiusInc; //sector radii
    //double maxRadius = max_radius(regNum);
    double minRadius = min_radius(regNum);
    //dfile << "region " << regNum << " maxRadius used " << maxRadius << " minRadius used " << minRadius <<endl;
    radiusInc = minRadius /(1.0*numSectors);
    double startAngle, finishAngle, angleInc; //sector angles
    angleInc = 2*PI/(1.*numSectors);
    startAngle = (beginAngle%numSectors)*angleInc;
    finishAngle = (endAngle%numSectors)*angleInc;
    int size = (int) regionIndex[regNum].size();
    // to normalise the NN field
     for (int i=0;i<size;i++){
          normalNN.push_back(this->NN[regionIndex[regNum][i]]);
      }
      normalNN = meanzero_vector(normalNN);
      //for (int i=0;i<size;i++)
          // dfile << " i " << i << " normalNN[i] " << normalNN[i] << endl;

    for (int k=0;k<numSectors;k++) {
      startRadius = (k*radiusInc);
      finishRadius = (k+1)*radiusInc;

      for (int i=0; i< size;i++) {
	if (this->H->at(regionIndex[regNum][i]).phi>=startAngle && this->H->at(regionIndex[regNum][i]).phi<finishAngle) {
	  if (this->H->at(regionIndex[regNum][i]).r>=startRadius && this->H->at(regionIndex[regNum][i]).r<finishRadius) {
	      radiusCount[k]++;
              radiusNN[k] += normalNN[i];
	  } //end of if on radius
        } //end of if on angleSector
      } //end of loop over i


      if (radiusCount[k] != 0)
	radiusNN[k] = radiusNN[k];
      else
	radiusNN[k] = 0.0;
      dfile << " region " << regNum << " startRadius "<<startRadius<<"  finishRadius "<<finishRadius<< " radiusNN " << radiusNN[k] << endl;
    }//end loop on k
     return radiusNN;

  } //end of function sectorize_radius

     //sectorize over radius
    //adapted for digital output
    vector <int> sectorize_reg_Dradius (int regNum, int numSectors, int beginAngle, int endAngle) {
        ofstream dfile ( "logs/sectorRadius.txt",ios::app);
    vector <int>  radiusNN;
    radiusNN.resize(numSectors,0);
    vector <int> radiusCount;
    vector <double> radiusHold;
    vector <double> normalNN;
    radiusCount.resize(numSectors,0);
    radiusHold.resize(numSectors,0);
    double startRadius, finishRadius, radiusInc; //sector radii
    double maxRadius = max_radius(regNum);
    double minRadius = min_radius(regNum);
    dfile << "region " << regNum << " maxRadius used " << maxRadius << " minRadius used " << minRadius <<endl;
    radiusInc = minRadius /(1.0*numSectors);
    double startAngle, finishAngle, angleInc; //sector angles
    angleInc = 2*PI/(1.*numSectors);
    startAngle = (beginAngle%numSectors)*angleInc;
    finishAngle = (endAngle%numSectors)*angleInc;
    int size = (int) regionIndex[regNum].size();
    for (int i=0;i<size;i++){
          normalNN.push_back(this->NN[regionIndex[regNum][i]]);
      }
    // to normalise the NN field
    normalNN = meanzero_vector(normalNN);
    double epsilon = 0.0001*(this->maxVal(normalNN) - this->minVal(normalNN));
      //for (int i=0;i<size;i++)
          // dfile << " i " << i << " normalNN[i] " << normalNN[i] << endl;

    for (int k=0;k<numSectors;k++) {
      startRadius = (k*radiusInc);
      finishRadius = (k+1)*radiusInc;
      for (int i=0; i< size;i++) {
        if (this->H->at(regionIndex[regNum][i]).phi >= startAngle && this->H->at(regionIndex[regNum][i]).phi<finishAngle) {
          if (this->H->at(regionIndex[regNum][i]).r >= startRadius && this->H->at(regionIndex[regNum][i]).r<finishRadius) {
              radiusCount[k]++;
              //radiusCC[k] += this->CC[regionIndex[regNum][i]];
              radiusHold[k] += normalNN[i];
          } //end of if on radius
        } //end of if on angleSector
      } //end of loop over i
    }//end of loop over k
    
    for (int k=0;k<numSectors;k++){
	    startRadius = (k*radiusInc);
	    finishRadius = (k+1)*radiusInc;
	    if (radiusCount[k] == 0) {
		    radiusNN[k] = 2;
		    continue;
	    }
//	    radiusHold[k]  = radiusHold[k]  / (1.*radiusCount[k]);
        if (radiusHold[k] > epsilon)
                radiusNN[k] = 1;
        else if (radiusHold[k] < epsilon)
                radiusNN[k] = -1;
        else
                radiusNN[k] = 0;

       dfile << " region " << regNum <<" startRadius "<<startRadius<<"  finishRadius "<<finishRadius<< " radiusNN " << radiusNN[k] << endl;
    }//end loop on k
     return radiusNN;

  } //end of function sectorize_radius



 //function to count the hexes in sectors of a region via angular sectors
    vector <double> sectorize_reg_angle (int regNum, int numSectors, int beginRadius, int endRadius) {
    //std::pair<double,double> diff; //difference between seed point and CoG of region
    ofstream cfile ("logs/sectorAngle.txt",ios::app);
    // diff = this->set_polars(regNum);
    vector <double> angleNN; //average value of CC in each sector
    vector <double> normalNN;
    vector <double>  angleVal;
    vector <int> count; //number of hexes in each sector
    angleNN.resize(numSectors,0);
    count.resize(numSectors,0);
    double startAngle, endAngle, angleInc; //sector angles
    double startRadius, finishRadius,radiusInc;
    double maxRadius = max_radius(regNum);
    // double minRadius = min_radius(regNum);
    radiusInc = maxRadius/ (1.0*numSectors);
    startRadius = beginRadius*radiusInc;
    finishRadius = endRadius*radiusInc;
    angleInc = 2*PI/(1.*numSectors);
    // to normalise the NN field
    int size = (int) regionIndex[regNum].size();
     for (int i=0;i<size;i++){
          normalNN.push_back(this->NN[regionIndex[regNum][i]]);
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
          if (this->H->at(regionIndex[regNum][i]).r >= startRadius && this->H->at(regionIndex[regNum][i]).r<finishRadius) {
              if (this->H->at(regionIndex[regNum][i]).phi >= startAngle && this->H->at(regionIndex[regNum][i]).phi<endAngle) {
                  angleVal.push_back(H->at(regionIndex[regNum][i]).phi);
                  count[k]++;
          //angle[k] += this->[regionIndex[regNum][i]];
                  angleNN[k] += normalNN[i];
              }//end if on angle
        //cfile << setw(5) << angleVal[i]  <<"  ";
          }//end if on radius
      }//end loop on i
    }//end loop on k

    angleNN = meanzero_vector(angleNN);
    for (int k=0;k<numSectors;k++){
	    startAngle = k*angleInc;
	    endAngle = k*angleInc;
	    if (count[k] != 0)
		    angleNN[k] = angleNN[k]/(1.*count[k]);
	    else
		    angleNN[k] = -999.999;
       cfile << " region " << regNum <<" startRadius "<<startRadius<<"  finishRadius "<<finishRadius<< " radiusNN " << angleNN[k] << endl;
    }//end loop on k

     
     return angleNN;

  } //end of function sectorize_region

 //function to count the hexes in sectors of a region via angular sectors
    vector <int> sectorize_reg_Dangle (int regNum, int numSectors, int beginRadius, int endRadius) {
    ofstream cfile ("logs/sectorAngle.txt");
    //std::pair<double,double> diff; //difference between seed point and CoG of region
    // diff = this->set_polars(regNum);
    vector <int> angleNN; //digitized value of NN in each sector
    vector <double> angleHold;
    vector <double> normalNN;
    vector <double>  angleVal;
    vector <int> angleCount; //number of hexes in each sector
    angleNN.resize(numSectors,0);
    angleHold.resize(numSectors,0);
    angleCount.resize(numSectors);
    double startAngle, endAngle, angleInc; //sector angles
    double startRadius, finishRadius,radiusInc;
    double maxRadius = max_radius(regNum);
    // double minRadius = min_radius(regNum);
    radiusInc = maxRadius/ (1.0*numSectors);
    startRadius = beginRadius*radiusInc;
    finishRadius = endRadius*radiusInc;
    angleInc = 2*PI/(1.*numSectors);
    // to normalise the NN field
    int size = (int) regionIndex[regNum].size();
     for (int i=0;i<size;i++){
          normalNN.push_back(this->NN[regionIndex[regNum][i]]);
      }
      normalNN = meanzero_vector(normalNN);
      double epsilon = 0.0001*(this->maxVal(normalNN) - this->minVal(normalNN));
     // for (int i=0;i<size;i++)
       //    cfile << " i " << i << " normalNN[i] " << normalNN[i] << endl;

    for (int k=0;k<numSectors;k++) {
      startAngle = k*angleInc;
      endAngle = (k+1)*angleInc;
      if ((k+1) == numSectors)
         endAngle = 2*PI;

      for (int i=0; i< size;i++) {
          if (this->H->at(regionIndex[regNum][i]).r >= startRadius && this->H->at(regionIndex[regNum][i]).r<finishRadius) {
              if (this->H->at(regionIndex[regNum][i]).phi >=startAngle && this->H->at(regionIndex[regNum][i]).phi<endAngle) {
                  angleVal.push_back(H->at(regionIndex[regNum][i]).phi);
                  angleCount[k]++;
          //angle[k] += this->[regionIndex[regNum][i]];
                  angleHold[k] += normalNN[i];
              }//end if on angle
        //cfile << setw(5) << angleVal[i]  <<"  ";
          }//end if on radius
      } //end if on i
    }//end if on k

      angleHold = meanzero_vector(angleHold);

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
       cfile << " region " << regNum <<" startRadius "<<startRadius<<"  finishRadius "<<finishRadius<< " radiusNN " << angleNN[k] << endl;

    } //end loop on k
     return angleNN;

  } //end of function sectorize_region
}; // Region
