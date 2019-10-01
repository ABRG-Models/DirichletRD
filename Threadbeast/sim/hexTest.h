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
  morph::HexGrid* Hgrid;
  //class constructor
  Region (int scale, string logpath) {
//        ofstream //afile (logpath + "/debug.out" );
//        ofstream jfile (logpath + "/correlateVector.out");

	srand(time(NULL));
    this->scale = scale;


    double s = pow(2.0, scale-1);
	ds = 1.0/s;

	
    n = 0;
    Hgrid = new morph::HexGrid(this->ds, 4.0, 0.0, morph::HexDomainShape::Boundary);
    morph::ReadCurves r("./barrelAE.svg");
    Hgrid->setBoundary (r.getCorticalPath());
    cout << "before filling H " << Hgrid->num() << endl;
    n = Hgrid->num();
    cout << " max x " << Hgrid->getXmax(0.0) << " min x " << Hgrid->getXmin(0.0) << endl; 
    cout << "after  filling H " << " n = " << n <<endl;


     cout << "before hex list loop" << endl;
     // this determines the type of hex
      for (auto h : Hgrid->hexen){

          if(h.boundaryHex == true) {
		    cout << " its a boundary hex i= " << h.vi << endl;
            if (!HAS_NE(h.vi)) {
				cout << "no ne" << endl;
				cout << " nbr hex " << this->Hgrid->d_ne[h.vi] <<endl;
            } 
            else {
			   cout << "has ne" << endl;
			   cout << " nbr hex " << this->Hgrid->d_ne[h.vi] <<endl;
            }
                  
            if (!HAS_NNE(h.vi)) {
				cout << "no nne" << endl;
			   cout << " nbr hex " << this->Hgrid->d_nne[h.vi] <<endl;
            } 
            else {
			   cout << "has nne" << endl;
			   cout << " nbr hex " << this->Hgrid->d_nne[h.vi] <<endl;
            }
            if (!HAS_NNW(h.vi)) {
				cout << "no nnw" << endl;
			   cout << " nbr hex " << this->Hgrid->d_nnw[h.vi] <<endl;
            } 
            else {
			   cout << "has nnw" << endl;
			   cout << " nbr hex " << this->Hgrid->d_nnw[h.vi] <<endl;
            }
            if (!HAS_NW(h.vi)) {
				cout << "no nw" << endl;
			   cout << " nbr hex " << this->Hgrid->d_nw[h.vi] <<endl;
            } 
            else {
			   cout << "has nw" << endl;
			   cout << " nbr hex " << this->Hgrid->d_nw[h.vi] <<endl;
            }
            if (!HAS_NSW(h.vi)) {
				cout << "no nsw" << endl;
			   cout << " nbr hex " << this->Hgrid->d_nsw[h.vi] <<endl;
            } 
            else {
			   cout << "has nsw" << endl;
			   cout << " nbr hex " << this->Hgrid->d_nsw[h.vi] <<endl;
            }
            if (!HAS_NSE(h.vi)) {
				cout << "no nse" << endl;
			   cout << " nbr hex " << this->Hgrid->d_nse[h.vi] <<endl;
            } 
            else {
			   cout << "has nse" << endl;
			   cout << " nbr hex " << this->Hgrid->d_nse[h.vi] <<endl;
            }
		  }
          else {
		        cout << "its not a boundary hex" << endl;
             } //end of if dealing with outer boundary
		}
		for (auto h : Hgrid->hexen) {
		  cout << "h.vi " << h.vi << endl;
		  }
  } // end of region constructor

  void listHex(){
  cout << "in listHex" << endl;
	for (auto h : Hgrid->hexen) {
	 cout << "h.vi " << h.vi << endl;
   }
}    
}; // end of region class
