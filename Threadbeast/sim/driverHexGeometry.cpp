#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "newregion.h"

using namespace std;

int main()
{
  //hexGeometry hGeo;
  DRegion M(7, "./logs");
  M.hGeo->value = 1;
  hexGeometry::point a,b;
  a.first = -1.0;
  a.second = -1.0;
  b.first = 1.0;
  b.second = 1.0;
  hexGeometry::lineSegment s;
  s.start = a;
  s.end = b;
  hexGeometry::line l1,l2;
  l1 = M.hGeo->segment2line(s);
  cout << "slope "<<l1.slope<<" intercept "<<l1.intercept<<endl;
  hexGeometry::point p;
  p.first = 1.0; p.second = 0.0;
  l2 = M.hGeo->perpPoint2Line(p, l1);
  cout << "slope "<<l2.slope<<" intercept "<<l2.intercept<<endl;
  cout << "value " << M.hGeo->value<<endl;
  hexGeometry::cIntersect cI;
  cI = M.hGeo->lineIntersect(l1,l2);
  cout << " intersect " << cI.intersect << " at " << "( "<<cI.p.first<<" , "<<cI.p.second<<" )" <<endl;
  cout << " dist p to s1 " << M.hGeo->pointLineDist(p,l1) << endl;
  hexGeometry::point h1;
  h1.first = -2.0;
  h1.second = 2.0;
  hexGeometry::segDist h1dist = M.hGeo->point2Seg(h1,s);
  cout << " point h1  inside " << h1dist.inside << " distance " << h1dist.dist << endl;
  int count = 0;
  for (auto h : M.Hgrid->hexen) {
    if (M.lineSegIntersectHex(s,h)) {
	  cout << " hex at ( " << M.Hgrid->d_x[h.vi] <<" , " << M.Hgrid->d_y[h.vi] <<" ) intersects slope 1 thru origin.)" << endl;
	  count++;
	  }
	}
	cout << count << " hexes intersected" <<endl;
 
  return 0;
}
