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
  double tol = 0.00001;

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

  struct horzLine {
    double slope = 0;
	double intercept;
	};

  struct vertLine {
    double slope = 999.999;
	double xintercept;
	};

  struct cIntersect {
  bool intersect;
  point p;
  };

  struct segDist {
    bool inside = false;
	double dist = 999.999;
	};

  int value;

  hexGeometry() {
   //int value = 0;
   //point start;
   //start.first = 0.0;
   //start.second = 0.0;
   }
  
  //returns line from lineSegment
  unsigned int  segmentLineType(lineSegment s) {
    unsigned int  result;
	double denominator = s.end.first - s.start.first;
	double numerator = s.end.first - s.start.first;
	if (denominator == 0)
	{
	  return result == 0;
	}
	else if ( numerator == 0)  
	{
	  return result = 1;
	}
	else 
	{
	  return result = 2;
	}

  }	

  // returns line from lineSegment, assumes have tested for verticality
    line segment2line(lineSegment s)
	{
      line result;
	  double denominator = s.end.first - s.start.first;
	  double numerator = s.end.first - s.start.first;
	  result.slope = numerator / denominator;
	  result.intercept = (s.end.first * s.start.second - s.end.second * s.start.first)/denominator;
	  return result;
	}

// returns a vertLine from a line segment
 vertLine segment2vert(lineSegment s) 
 {
   vertLine result;
   result.xintercept = s.end.first;
   return result;
 }

 horzLine segment2horz(lineSegment s) {
   horzLine result;
   result.intercept = s.end.second;
   return result;
 }

	//returns lineSegment from two points
	lineSegment createLineSegment(point a, point b)
	{ 
	  lineSegment result;
	  result.start = a;
	  result.end = b;
	  return result;
	}

	bool pointOnLine(point p, line l) {
	  bool result;
	  double residual = p.second - p.first * l.slope - l.intercept;
	  if (residual < tol)
	  {
	    result = 1;
	  }
	  else 
	  {
	    result = 0;
	  }
	  return result;
	}

	bool pointOnhorzLine(point p, horzLine h) 
	{
	  bool result;
	  if ( (p.second - h.intercept) < tol)
	  {
	    result = 1;
	  }
	  else 
	  {
	    result = 0;
	  }
	  return result;
	}  
	
	bool pointOnvertLine(point p, vertLine h) 
	{
	  bool result;
	  if ( (p.first - h.xintercept) < tol)
	  {
	    result = 1;
	  }
	  else 
	  {
	    result = 0;
	  }
	  return result;
	}  



//returns line through point p perpendicular to line l
   line perpPoint2Line (point p, line l) {
     line result;
	   result.slope = - 1.0 / l.slope;
	   result.intercept = p.second - result.slope * p.first;
	  return result;
	}

// returns horizontal line through point perpendicular to vertLine
  horzLine perpPoint2vertLine( point p, vertLine h) {
    horzLine result;
	result.intercept = p.second;
	return result;
  }

// returns a vertLine through a point perpendicular to horzLine
  vertLine perpPoint2horzLine(point p, horzLine h) 
  {
    vertLine result;
	result.xintercept = p.first;
	return result;
  }


//returns a structure giving the result of the intersection of two lines.
	cIntersect lineIntersect (line l1, line l2) {
	  cIntersect result;
	  if (l1.slope == l2.slope) {
	    result.intersect = 0;
		result.p.first = -999.999;
		result.p.second = -999.999;
		}
	   else {
	     result.intersect = 1;
		 result.p.first = (l2.intercept - l1.intercept) / (l1.slope - l2.slope);
		 result.p.second = (l1.slope*l2.intercept - l1.intercept*l2.slope) / (l1.slope - l2.slope);
        }
		return result;
     }

// returns a structure giving the result of the intersection of a line and a horzLine
  cIntersect horzIntersect (horzLine h, line l) {
    cIntersect result;
	result.intersect = 1;
	result.p.first = (h.intercept - l.intercept) / l.slope;
	result.p.second = h.intercept;
	return result;
 }

// returns a structure giving the result of the intersection of a line and a vertLine
  cIntersect vertIntersect (vertLine h, line l) {
    cIntersect result;
	result.intersect = 1;
	result.p.first = h.xintercept;
	result.p.second = h.xintercept * l.slope + l.intercept;
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

//gives the distance of a point from a vertLine
  double pointVertDist(point p1, vertLine h1) {
    double result;
	result = fabs(p1.first - h1.xintercept);
	return result;
  }

 //gives the distance of a point from a horizontal line
   double pointHorzDist(point p1, horzLine h1) {
     double result;
	 result = fabs(p1.second - h1.intercept);
	 return result;
   }

// returns a structure that determines distance of the point from the line derived from a 
// segment and if the point is within the rectangle defined by the segment
	 segDist point2Seg (point p, lineSegment s) {
	   segDist result;
	   double segX = fabs(s.start.first - s.end.first);
	   double segY = fabs(s.start.second - s.end.second);
	   //double pX = fmax(fabs(p.first - s.end.first),fabs(p.first - s.start.first));
	   //double pY = fmax(fabs(p.second - s.end.second),fabs(p.second - s.start.second));
	   double pX = (p.first - s.end.first)*(p.first - s.start.first);
	   double pY = (p.second - s.end.second)*(p.second - s.start.second);
	     
	   line l1 = segment2line(s);
	   //if ((pX <= segX) || (pY <= segY)) {
	   if (segX < this->tol) {
		 if (pY < 0.0) 
		   result.inside = 1;
		 else 
		   result.inside = 0;
	     vertLine h1 = segment2vert(s);
		 result.dist = pointVertDist(p,h1);
	   }
	   else if (segY < this->tol) {
		 if (pX < 0.0) 
		   result.inside = 1;
		 else 
		   result.inside = 0;
	     horzLine h1 = segment2horz(s);
		 result.dist = pointHorzDist(p,h1);
	   }
	   else {
	   if ((pX <= 0.0 ) || (pY <= 0.0)) {
	     result.inside = 1;
		 }
	   else {
	     result.inside = 0;
		 }
	   result.dist = pointLineDist(p,l1);
	   }
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

//private:
 
//distance between two points
double  getdist(point p1, point p2) {
   double result;
   result = sqrt((p1.first - p2.first)*(p1.first - p2.first) + (p1.second-p2.second)*(p1.second-p2.second));
   return result;
   }
	 

	 

 
	 
};


