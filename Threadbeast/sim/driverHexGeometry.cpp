#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "region.h"
#include "analysis.h"
#include "ksSolver.h"
#include <morph/display.h>
#include <morph/Config.h>
#include <cctype>
#include <locale>
#include <algorithm>
#include <string>
using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::runtime_error;
using morph::Gdisplay;

using namespace std;
using namespace morph;

int main(int argc, char **argv)
{
  //hexGeometry hGeo;
     bool Lsequence = atoi(argv[1]);
     if (argc < 2) {
          std::cout << "You need to set the bool arg to determine if the setting of the region boundaries is singular or in sequence" << endl;
          return -1;
     }
  DRegion M(5, 5, "./logs");
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
  double tol = 0.1;
  bool result;
  vector<hexGeometry::point> invect;
  vector<hexGeometry::lineSegment> triangle;
  std::pair tPoint (0.0, 0.0);
  hexGeometry::point testPoint;
  morph::BezCurvePath<float> rBoundary;
  testPoint.first = tPoint.first;
  testPoint.second = tPoint.second;
  std::pair p1 (0.866, 0.5);
  std::pair p2 (0.0, 1.0);
  std::pair p3 (-0.866, 0.5);
  std::pair p4 (-0.866 , -0.5);
  std::pair p5 (0.0 , -1.0);
  std::pair p6 (0.866 , -0.5);
  morph::BezCurve<float> bc1(p2,p1);
  rBoundary.addCurve(bc1);
  morph::BezCurve<float> bc2(p3,p2);
  rBoundary.addCurve(bc2);
  morph::BezCurve<float> bc3(p4,p3);
  rBoundary.addCurve(bc3);
  morph::BezCurve<float> bc4(p5,p4);
  rBoundary.addCurve(bc4);
  morph::BezCurve<float> bc5(p6,p5);
  rBoundary.addCurve(bc5);
  morph::BezCurve<float> bc6(p1,p6);
  rBoundary.addCurve(bc6);

  invect.push_back(M.hGeo->pair2point(p1));
  invect.push_back(M.hGeo->pair2point(p2));
  invect.push_back(M.hGeo->pair2point(p3));
  triangle.push_back(M.hGeo->createLineSegment(M.hGeo->pair2point(p2),M.hGeo->pair2point(p1)));
  triangle.push_back(M.hGeo->createLineSegment(M.hGeo->pair2point(p3),M.hGeo->pair2point(p2)));
  triangle.push_back(M.hGeo->createLineSegment(M.hGeo->pair2point(p4),M.hGeo->pair2point(p3)));
  triangle.push_back(M.hGeo->createLineSegment(M.hGeo->pair2point(p5),M.hGeo->pair2point(p4)));
  triangle.push_back(M.hGeo->createLineSegment(M.hGeo->pair2point(p6),M.hGeo->pair2point(p5)));
  triangle.push_back(M.hGeo->createLineSegment(M.hGeo->pair2point(p1),M.hGeo->pair2point(p6)));
  vector<double> angles;
  int size = 3;
  double windingNumber = 0;
  double minAngle = 2 * PI;
  int indexCount = 0;
  for (unsigned int i=0; i<size; i++) {
      hexGeometry::lineSegment firstSide = M.hGeo->createLineSegment(testPoint,invect[i]);
      hexGeometry::dLine dFirstSide = M.hGeo->segment2dLine(firstSide);
      if (dFirstSide.angle < minAngle) {
          minAngle = dFirstSide.angle;
          indexCount = i;
      }
  }
  cout << " minimum angle " << minAngle << " index " << indexCount;

   for (unsigned int i=indexCount; i < size + indexCount; i++) {
       hexGeometry::lineSegment tempSide = M.hGeo->createLineSegment(testPoint,invect[i%size]);
        hexGeometry::dLine dSide = M.hGeo->segment2dLine(tempSide);
        double correctedAngle = dSide.angle - minAngle;
        angles.push_back(correctedAngle);
    }
    cout << " before windingNumber loop" << endl;
    for (unsigned int i=0; i<size;i++){
        unsigned int lead = (i+1) % size;
        if (angles[lead] == 0.0) {
            angles[lead] += 2.0 * PI;
        }
        windingNumber += angles[lead] - angles[i];
        cout << " windinNumber " << windingNumber << " lead " << angles[lead] << " follow " << angles[i] << endl;
    }
    cout << " after windingNumber loop" << endl;
    windingNumber = windingNumber / (2.0 * PI);
    if (((1.0 - tol) < windingNumber) && (windingNumber < (1.0 + tol))) {
        result = true;
    }
    else {
        result = false;
    }
    cout << " Winding Number " << " is " << windingNumber << " in region bool " << result << endl;

    ksSolver S;
    double radius = 1.0;
    std::pair<double, double>  centroid = std::make_pair(0.0, 0.0);
    S = ksSolver(5, 5.0, "./", radius, centroid);

    vector<double> fix(3.0, 0.0);
    vector<double> eye(3.0, 0.0);
    vector<double> rot(3.0, 0.0);
    double rhoInit = 3.0;
    //double rhoInit = 5.0;
    array<float,3> cl_a = morph::Tools::getJetColorF (0.78);
    array<float,3> cl_c = morph::Tools::getJetColorF (0.28);
    array<float,3> cl_b = morph::Tools::getJetColorF (0.58);
    array<float,3> cl_d = morph::Tools::getJetColorF (0.00);
    array<float,3> offset = {{0, 0, 0}};
    if (Lsequence) {
        morph::Gdisplay isp(900, 900, 0, 0, "A boundary", rhoInit, 0.0, 0.0);
        cout << "after setting display" << endl;
        isp.resetDisplay (fix, eye, rot);
        cout << "after setting display" << endl;
        int triangleSize = triangle.size();

        for (auto& h : S.Hgrid->hexen) {
            for (unsigned int i=0; i<triangleSize; i++) {
                hexGeometry::point p;
                //p.first = M.Hgrid->d_x[h.vi];
                p.first = h.position()[0];
                p.second = h.position()[1];
                // p.second = M.Hgrid->d_y[h.vi];
                // if (M.hGeo->hexIntersectLineSegment(regionLineSegments[i], p, hexWidth)) {
                if (M.hGeo->hexIntersectLineSegment(triangle[i],h)){
                    h.setFlags(HEX_IS_REGION_BOUNDARY);
                    cout << " hex.x " << p.first << " hex.y " << p.second <<  " hex id " << h.vi << " side " << i << endl;
                }
            }
        }
        int boundaryCount = 0;
        for (auto h : S.Hgrid->hexen) {
            if (h.testFlags(HEX_IS_BOUNDARY) == true) {
                isp.drawHex (h.position(), (h.d/2.0f), cl_a);
                boundaryCount++;
            }
        }

        boundaryCount = 0;
        for (auto h : S.Hgrid->hexen) {
            if (h.testFlags(HEX_IS_REGION_BOUNDARY) == true) {
                isp.drawHex (h.position(), (h.d/2.0f), cl_c);
                boundaryCount++;
            }
        }

        cout << "boundaryCount "<<boundaryCount<<endl;
        usleep (10000000); // ten seconds
        isp.redrawDisplay();
        usleep (10000000); // ten seconds
        isp.saveImage("./Tesselation.png");
        isp.closeDisplay();
    }

    for (auto& h : S.Hgrid->hexen) {
        if (h.testFlags(HEX_IS_BOUNDARY) == true) {
            h.resetUserFlags();
        }
    }
    morph::Gdisplay misp(900, 900, 0, 0, "A boundary", rhoInit, 0.0, 0.0);
    cout << "after setting display" << endl;
    misp.resetDisplay (fix, eye, rot);
    cout << "after setting display" << endl;

    std::vector<std::list<morph::Hex>::iterator> regionpHex;
    std::pair<float,float> regionCentroid;
    regionpHex = S.Hgrid->getRegion(rBoundary, regionCentroid, true);
    //plot stuff here.
    int boundaryCount = 0;
    for (auto h : S.Hgrid->hexen) {
        if (h.testFlags(HEX_IS_BOUNDARY) == true) {
            misp.drawHex (h.position(), (h.d/2.0f), cl_a);
            boundaryCount++;
        }
    }
    boundaryCount = 0;
    for (auto h : regionpHex) {
     //if (M.Creg[h.vi] ==  1 ) {
        if (h->testFlags(HEX_IS_REGION_BOUNDARY) == true) {
            misp.drawHex (h->position(), (h->d/2.0f), cl_c);
            boundaryCount++;
        }
    }
    cout << "boundaryCount "<<boundaryCount<<endl;
usleep (10000000); // ten seconds
    misp.redrawDisplay();
    usleep (10000000); // ten seconds
    misp.saveImage("./Tesselation1.png");
    misp.closeDisplay();
    return 0;
}
