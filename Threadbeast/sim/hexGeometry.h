/* 
 * hexGeometry class
 * Author: John Brooke
 *
 * Date 2019/10
 *
 * Helper class for HexGrid and Dregion 
 * to be included from Dregion 
 * (maybe change to HexGrid later)
 */
using morph::HexGrid;
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
 
class hexGeometry 
{
public: 
  struct point {
    double first;
	double second;

 };


  struct  lineSegment {
    point start;
	point end;
  };

  struct line {
    double slope;
	double intercept;
  };

  struct cIntersect {
  bool intersect;
  point p;
  };

  struct segDist {
    bool inside;
	double dist;
	};

  int value;

  hexGeometry() {
   int value = 0;
   point start;
   start.first = 0.0;
   start.second = 0.0;
   }
  
  
  //returns line from lineSegment
  line  segment2line(lineSegment s) {
    line result;
	result.slope = ((s.start).second - (s.end).second)/((s.start).first - (s.end).first);
	result.intercept = (s.start).second - result.slope*(s.start).first;
	return result;
	}

//returns line through point p perpendicular to line l
   line perpPoint2Line (point p, line l) {
     line result;
     if (l.slope == 0) {
	   result.slope = 1;
	   result.intercept = p.second;
	   }
	 else {
	   result.slope = - 1.0 / l.slope;
	   result.intercept = p.second - result.slope * p.first;
	  }
	  return result;
	}

//returns a structure giving the result of the intersection of two lines.
	cIntersect lineIntersect (line l1, line l2) {
	  cIntersect result;
	  if (l1.slope == l2.slope) {
	    result.intersect = 0;
		result.p.first = -100000.0;
		result.p.second = -100000.0;
		}
	   else {
	     result.intersect = 1;
		 result.p.first = (l2.intercept - l1.intercept) / (l1.slope - l2.slope);
		 result.p.second = (l1.slope*l2.intercept - l1.intercept*l2.slope) / (l1.slope - l2.slope);
        }
		return result;
     }

//gives the distance of a point from a line. 
	 double pointLineDist (point p1, line l1) {
	   double result;
	   line l2 = perpPoint2Line(p1,l1);
	   cIntersect cI = lineIntersect(l1, l2);
	   point p2 = cI.p;
	   result = getdist(p1,p2);
	   return result;
	   }

// returns a structure that determines distance of the point from the line derived from a 
// segment and if the point is within the rectangle defined by the segment
	 segDist point2Seg (point p, lineSegment s) {
	   segDist result;
	   line l1 = segment2line(s);
	   double segX = fabs(s.start.first - s.end.first);
	   double segY = fabs(s.start.second - s.end.second);
	   double pX = fabs(p.first - s.end.first);
	   double pY = fabs(p.second - s.end.second);
	   if ((pX <= segX) && (pY <= segY)) {
	     result.inside = 1;
		 }
	   else {
	     result.inside = 0;
		 }
	   result.dist = pointLineDist(p,l1);
	   return result;
	  }
/* now moved to DRegion
	 bool lineIntersectHex(lineSegment s, Hex h, HexGrid hGrid) {
	   bool result;
	   double d = h.getD();
	   point hexCentre;
	   hexCentre.first = hGrid->d_x(h.vi);
	   hexCentre.second = hGrid->d_y(h.vi);
	   result = point2Seg(hexCentre,s).inside;
	   return result;
	   }
*/

private:
 
//distance between two points
double  getdist(point p1, point p2) {
   double result;
   result = sqrt((p1.first - p2.first)*(p1.first - p2.first) + (p1.second-p2.second)*(p1.second-p2.second));
   return result;
   }
	 

	 

 
	 
};


