#include <morph/tools.h>
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

using std::vector;
using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::runtime_error;
using morph::Tools;
using namespace std;

// to find the maximum value of a vector
class Analysis {
//global class variables
public:
  struct extremum {
  int radialIndex;
  double radialValue;
  };
  vector<extremum> turnVal; //radial turning points

  //default constructor

  Analysis () {};

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

    // return the mean of a vector
        double mean_vector(vector<double> invector) {
        double result;
        int size = invector.size();
        //meanzero << "size " << size << endl;
        double sum = 0;
        for (int i=0; i <size; i++) {
            sum += invector[i];
        }
        sum = sum/(1.0*size);
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

	vector <double> normalise (vector <double> invector) {
	  vector <double> result;
	  unsigned int size = invector.size();
	  result.resize(size, 0);
	  double maxV = -1e7;
	  double minV = 1e7;
	  for (unsigned int i=0;i<size;i++) {
	    if (invector[i] > maxV) {maxV = invector[i];}
		if (invector[i] < minV) {minV = invector[i];}
		}
	  double scaleV = 1./(maxV - minV);
	  for (unsigned int i=0;i<size;i++) {
	    result[i] = fmin(fmax((invector[i] - minV)*scaleV,0.),1.);
	  }
	  return result;
    }




  //function find_max to find turning points both values and indices.
    int find_max(vector<double> ray, int window) {
    int size = ray.size();
    // ofstream dfile ("turn.txt",ios::app);
    //dfile <<"size = " << size <<endl;
    vector<double> smoothRay;
    smoothRay = this->smooth_vector(ray, window);
    turnVal.resize(1000);
    //cout <<" "<<iend<<iend<<flush;
    double old_slope = 0;
    double new_slope = 0;
    int count = 0;
    old_slope = smoothRay[1] - smoothRay[0];
    for (int i =2; i<=size+1;i++){

      new_slope = smoothRay[i%size]-smoothRay[(i-1)%size];
      //dfile << " " << i%size << " " << old_slope << " "<<new_slope <<endl;
      if (new_slope*old_slope < 0.) {
	turnVal[count].radialIndex = i;
	turnVal[count].radialValue = smoothRay[i]; //should really interpolate
      //  dfile << "turn index " << turnVal[count].radialIndex << " turn value " << turnVal[count].radialValue << endl;
        count++;
      }
      old_slope = new_slope;
    }
    return count;
  }
// find the zeros in a ray angular
    int find_zeroDAngle(vector<int> ray) {
    int size = ray.size();
    //ofstream zerofile ("zero.txt",ios::app);
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
      //zerofile << " radius " << i%size << " " << old_val << " "<<new_val <<endl;
      }
    return count;
  }


 // find the zeros in a ray angular
    int find_zeroDRadius(vector<int> ray) {
    int size = ray.size();
    //ofstream zerofile ("zero.txt",ios::app);
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
      //zerofile << " angle " << i%size << " " << old_val << " "<<new_val <<endl;
      }
    return count;
  }


// find the zeros in a ray angular
    int find_zeroRadius(vector<int> ray) {
    int size = ray.size();
    //ofstream zerofile ("zero.txt",ios::app);
    int old_val = 0;
    int new_val = 0;
    int count = 0;
    for (int i =0; i<size;i++){
        if (ray[i%size] != 0) {
           new_val = ray[i%size];
           if (old_val*new_val == -1) {
                   count++;
           }
	   old_val = new_val;
        }
      //zerofile << " " << i%size << " " << old_val << " "<<new_val <<endl;
      }
    return count;
  }


    // find the zeros in a ray angular
    int find_zeroAngle(vector<double> ray, int window) {
    int size = ray.size();
    // ofstream zerofile ("zero.txt",ios::app);
    vector<double> smoothRay;
    //smoothRay = this->smooth_vector(ray, window);
    smoothRay = ray;
    //zerofile <<"size = " << size <<endl;
    turnVal.resize(1000);
    double old_val = 0;
    double new_val = 0;
    int count = 0;
    old_val = smoothRay[0];
    for (int i =1; i<size+1;i++){
      new_val = smoothRay[i%size];
      //zerofile << " " << i%size << " " << old_val << " "<<new_val <<endl;
      if (new_val*old_val < 0.) {
	turnVal[count].radialIndex = i;
	turnVal[count].radialValue =smoothRay[i]; //should really interpolate
        //zerofile << "turn index " << turnVal[count].radialIndex << " turn value " << turnVal[count].radialValue << endl;
        count++;
      }
      old_val = new_val;
    }
    return count;
  }

      // find the zeros in a radial ray
    int find_zeroRadius(vector<double> ray, int window) {
    int size = ray.size();
    //ofstream zerofile ("zero.txt",ios::app);
     vector<double> smoothRay;
    // smoothRay = this->smooth_vector(ray, window);
     smoothRay = ray;
    //zerofile <<"size = " << size <<endl;

    turnVal.resize(1000);
    double old_val = 0;
    double new_val = 0;
    int count = 0;
    old_val = smoothRay[0];
    for (int i =1; i<size;i++){
      new_val = smoothRay[i%size];
      //zerofile << " " << i%size << " " << old_val << " "<<new_val <<endl;
      if (new_val*old_val < 0.) {
	turnVal[count].radialIndex = i;
	turnVal[count].radialValue = smoothRay[i]; //should really interpolate
        //zerofile << "turn index " << turnVal[count].radialIndex << " turn value " << turnVal[count].radialValue << endl;
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
}; //end of class Analysis
