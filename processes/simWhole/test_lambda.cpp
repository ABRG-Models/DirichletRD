
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <math.h>
#include <random>
#include <algorithm>
#define PI 3.14159265




//using namespace std;
   


    template <typename T>
    std::vector<size_t> sort_indexes(const std::vector<T> &v) {
     int start=0;
     int end=10;
     
     // initialize original index locations
     std::vector<size_t> idx(v.size());
     iota(idx.begin(), idx.end(), 0);

     // sort indexes based on comparing values in v
     sort(idx.begin(), idx.end(),
        [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

     return idx;
    };

  using namespace std;
  int main (int argc, char **argv) {
     int start=0;
     int end=10;
     std::default_random_engine generator;
     std::uniform_int_distribution<int> distribution(0,10);

   
    //sort_indexes<vector <int> >;
    
    
     std::vector <int> v;
     v.resize(10);

    for (int i=0;i<10;i++){
      v[i] = distribution(generator);
      cout <<v[i] << "   ";
     
    }
     cout << endl; cout << endl;

     int count=0;
     vector <int> sortedIndex;
     sortedIndex.resize(10);
     int temp;
     vector size_t myindices = sort_indexes(v);
     for (size_t & i: myindices) {
      //cout << v[i] << "  ";
      temp = i;
      sortedIndex[count] = temp;
      cout << "i= " << i << "    v[i]= " << v[i] << "  sortedIndex= " << temp <<  endl;
      count++;
     }
     cout << endl;

     for (int j=0;j<10;j++)
     cout<<"sortedIndex[j] = "<< sortedIndex[j]<<endl;
     
  }
