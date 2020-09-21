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
#define PI 3.1415926535897932

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
    double tol = 0.0001;

    enum lineType  {VERTICAL, HORIZONTAL, SLANTED};

    const double INFINITE = 9999999.9999999;

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

    /*
     * a directed line, angle is angle to x axis
     * intercept is cut on y axis if slope not infinite
     * intercept is cut on x axis if angle = +- pi/2
     */
    struct dLine {
        double angle;
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


    //constructor
    hexGeometry() {
    }


    point pair2point(std::pair<double,double> inpair) {
        point result;
        result.first = inpair.first;
        result.second = inpair.second;
        return result;
    }

    std::pair<double, double> point2pair(point inPoint) {
        std::pair<double, double> result;
        result.first = inPoint.first;
        result.second = inPoint.second;
        return result;
    }


    /*!
     * Euclidean distance between two points
     */
    double distance(point a, point b) {
        double result;
        result = sqrt((a.first-b.first)*(a.first-b.first) + (a.second-b.second)*(a.second-b.second));
        return result;
    }

  //returns linetype from lineSegment
    lineType segmentLineType(lineSegment s) {
        lineType  result;
        double denominator = s.end.first - s.start.first;
        double numerator = s.end.first - s.start.first;
        if (denominator == 0)
        {
            return result = VERTICAL;
        }
        else if ( numerator == 0)
        {
            return result = HORIZONTAL;
        }
        else
        {
            return result = SLANTED;
        }

    }

    lineType lineLineType (line l) {
        lineType result;
        if (l.slope == INFINITE) {
            result = VERTICAL;
        }
        else if (l.slope == 0) {
            result = HORIZONTAL;
        }
        else {
            result = SLANTED;
        }
        return result;

    }
  // returns line from lineSegment, tests for verticality and horizontality
    line segment2line(lineSegment s)
	{
        line result;
        double denominator = s.end.first - s.start.first;
        double numerator = s.end.second - s.start.second;
        if (fabs(denominator) <= this->tol) {
           result.slope = INFINITE;
           result.intercept = s.end.first;
           return result;
        }
        else if (fabs(numerator) < this->tol) {
            result.slope = 0.0;
            result.intercept = s.end.second;
            return result;
        }
        else {
            result.slope = numerator / denominator;
            result.intercept =  (s.end.first * s.start.second - s.end.second * s.start.first)/denominator;
            return result;
        }
    }

    dLine segment2dLine (lineSegment s) {
        dLine result;
        double angle, intercept;
        double denominator = s.end.first - s.start.first;
        double numerator = s.end.second - s.start.second;
        if (numerator >= 0) {
            angle = atan2(numerator, denominator);
        }
        else {
            angle = atan2(numerator, denominator) + 2.0 * PI;
        }

        if (denominator == 0.0) {
            intercept = s.end.first;
        }
        else if (numerator >= 0.0) {
            intercept = s.end.second - s.end.first * tan(angle);
        }
        else {
            intercept = s.start.first + s.start.first * tan(angle);
        }
        result.angle = angle;
        result.intercept = intercept;
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

    //returns the point that is the start of a lineSegment
    point startSegment(lineSegment s) {
        return s.start;
    }

    //returns the point that is the end of a lineSegment
    point endSegment(lineSegment s) {
        return s.end;
    }

    //returns a signed angle between 2 lines
    double subtendLines( dLine start , dLine end) {
        return start.angle - end.angle;
    }

    //returns signed angle between two segments
    double subtendSegments(lineSegment s1, lineSegment s2) {
        dLine dL1 = segment2dLine(s1);
        dLine dL2 = segment2dLine(s2);
        return subtendLines(dL1,dL2);
    }

    //true if point is left of a lineSegment
    bool pointLeftSegment(point p, lineSegment s) {
        double phi, psi;
        lineSegment sp = createLineSegment(startSegment(s), p);
        psi = segment2dLine(s).angle;
        phi = segment2dLine(sp).angle;
        if (phi <= PI) {
            if ((phi <= psi) && (psi <= phi + PI)) {
                return true;
            }
            else {
                return false;
            }
        }
        else {
            if ((phi <= psi && psi <= 2*PI) || (psi <= phi - PI)) {
                return true;
            }
            else {
                return false;
            }
        }
    }


//returns a bool if point is on line or not
    bool pointOnLine(point p, line l) {
        bool result = false;
        double residual = p.second - p.first * l.slope - l.intercept;
        if (lineLineType(l) == VERTICAL && (fabs(p.first - l.intercept) < this->tol)) {
                result = true;
                return result;
        }
        if (residual < this->tol)
        {
           result = true;
        }
       return result;
    } //end of pointOnLine

    bool pointOnhorzLine(point p, horzLine h)
	{
        bool result;
        if ( (p.second - h.intercept) < this->tol)
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
        if ( fabs((p.first - h.xintercept)) < this->tol)
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
        lineType lt1 = lineLineType(l1);
        lineType lt2 = lineLineType(l2);
        if ((l1.slope == l2.slope) && (l1.intercept == l2.intercept)) {
            if (lt2 == VERTICAL && lt1 == VERTICAL) {
            //    std::cout << " lines vertical same intercept " << endl;
            }
            result.intersect = true;
            result.p.first = INFINITE;
            result.p.second = INFINITE;
            return result;
		}
        else if ((l1.slope == l2.slope) && (l1.intercept != l2.intercept)) {
            if (lt2 == VERTICAL && lt1 == VERTICAL) {
            //    std::cout << " lines vertical different intercept " << endl;
            }
            result.intersect = false;
            result.p.first = INFINITE;
            result.p.second = INFINITE;
            return result;
        }
        else if (lt1 == VERTICAL && lt2 != VERTICAL) {
            result.intersect = true;
            result.p.first = l1.intercept;
            // std::cout << " hex side is vertical " << l1.intercept << endl;
            result.p.second = l2.slope * l1.intercept + l2.intercept;
            return result;
        }
        else if (lt1 != VERTICAL && lt2 == VERTICAL) {
            result.intersect = true;
            result.p.first = l2.intercept;
            // std::cout << " hexagon side is vertical " << l2.intercept << endl;
            result.p.second = l1.slope * l2.intercept + l1.intercept;
            return result;
        }
        else {
            result.intersect = true;
            result.p.first = (l2.intercept - l1.intercept) / (l1.slope - l2.slope);
            result.p.second = (l1.slope*l2.intercept - l1.intercept*l2.slope) / (l1.slope - l2.slope);
            //std::cout << " neither side  is vertical " << l2.intercept << endl;
            return result;
        }
    }

    /*!
     * returns true if point is within the box defined by a line segment
     */
    bool pointInLineSegmentBox(point p, lineSegment s) {
        bool result = false;
        double xsign = (p.first - s.start.first) * (p.first - s.end.first);
        double ysign = (p.second - s.start.second) * (p.second - s.end.second);
        if (s.start.first > s.end.first - this->tol && s.start.first < s.end.first + this->tol) {
            xsign = 0.0;
            ysign = (p.second - s.start.second) * (p.second - s.end.second);
        }
        else if (s.start.second > s.end.second - this->tol && s.start.second < s.end.second + this->tol) {
            xsign = (p.first - s.start.first) * (p.first - s.end.first);
            ysign = 0.0;
        }
        else {
            xsign = (p.first - s.start.first) * (p.first - s.end.first);
            ysign = (p.second - s.start.second) * (p.second - s.end.second);
        }
        if ((xsign <= 0.0) && (ysign <= 0.0)) {
            result = true;
        }
        return result;
    }

    /*!
     * returns with bool member of cIntersect true if line segments intersect with their
     * point of intersection inside both bounding boxes. Point member given by the
     * line intersection.
     */
    cIntersect lineSegmentIntersect (lineSegment s1, lineSegment s2) {
        cIntersect result;
        line l1 =  segment2line(s1);
        line l2 = segment2line(s2);
        result = lineIntersect(l1, l2);
        //return result;
        if (result.intersect == true && pointOnLine(result.p,l1) && pointOnLine(result.p,l2)) {
            return result;
        }
        else {
            result.intersect = false;
            return result;
        }
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

     std::vector<lineSegment> hexSides(morph::Hex h) {
        vector<lineSegment> result;
        result.resize(6);
        double lr = static_cast<double> (h.getLR());
        double sr = static_cast<double> (h.getSR());
        point point1, point2;
        point1.first = h.x + sr; point1.second = h.y + lr/2.0;
        point2.first = h.x; point2.second = h.y + lr;
        result[0] = createLineSegment(point1, point2);
        point1 = point2;
        point2.first = h.x - sr; point2.second = h.y + lr/2.0;
        result[1] = createLineSegment(point1, point2);
        point1 = point2;
        point2.first = h.x - sr; point2.second = h.y - lr/2.0;
        result[2] = createLineSegment(point1, point2);
        point1 = point2;
        point2.first = h.x; point2.second = h.y - lr;
        result[3] = createLineSegment(point1, point2);
        point1 = point2;
        point2.first = h.x + sr; point2.second = h.y - lr/2.0;
        result[4] = createLineSegment(point1, point2);
        point1 = point2;
        point2.first = h.x + sr; point2.second = h.y + lr/2.0;
        result[5] = createLineSegment(point1, point2);
        //now return the vector of line segments
        return result;
    }


     std::vector<lineSegment> hexSides(point p, double d) {
        vector<lineSegment> result;
        result.resize(6);
        double lr = 2.0 * d / sqrt(3.0);
        double sr = d;
        point point1, point2;
        point1.first = p.first + sr; point1.second = p.second + lr/2.0;
        point2.first = p.first; point2.second = p.second + lr;
        result[0] = createLineSegment(point1, point2);
        point1 = point2;
        point2.first = p.first - sr; point2.second = p.second + lr/2.0;
        result[1] = createLineSegment(point1, point2);
        point1 = point2;
        point2.first = p.first - sr; point2.second = p.second - lr/2.0;
        result[2] = createLineSegment(point1, point2);
        point1 = point2;
        point2.first = p.first; point2.second = p.second - lr;
        result[3] = createLineSegment(point1, point2);
        point1 = point2;
        point2.first = p.first + sr; point2.second = p.second - lr/2.0;
        result[4] = createLineSegment(point1, point2);
        point1 = point2;
        point2.first = p.first + sr; point2.second = p.second + lr/2.0;
        result[5] = createLineSegment(point1, point2);
        //now return the vector of line segments
        return result;
    }

    bool pointInHexBox (point p, morph::Hex h) {
        bool result;
        double upperx = h.x + h.getSR() + this->tol;
        double lowerx = h.x - h.getSR() - this->tol;
        double uppery = h.y + h.getLR() + this->tol;
        double lowery = h.y - h.getLR() - this->tol;
        result = (p.first >= lowerx && p.first <= upperx && p.second >= lowery && p.second <= uppery);
        return result;
    }



    bool hexIntersectLineSegment(lineSegment s, point p, double d) {
        bool result = false;
        vector<lineSegment> hexSides = this->hexSides(p, d);
        for (int i=0; i<6; i++) {
            if (lineSegmentIntersect(hexSides[i], s).intersect) {
                result = true;
                break;
            }
        }
        return result;
    }

    bool hexIntersectLineSegment(lineSegment s, morph::Hex h) {
        bool result = false;
        point pt;
        pt.first = h.x;
        pt.second = h.y;
        vector<line> hexLines;
        for (int i=0; i<6; i++) {
            hexLines.push_back(this->segment2line(this->hexSides(h)[i]));
        }
        line l = this->segment2line(s);
        int intersected = 0;
        for (int i=0; i<6; i++) {
            pt = lineIntersect(hexLines[i], l).p;
            if (lineIntersect(hexLines[i],l).intersect && pointInHexBox(pt, h) && pointInLineSegmentBox(pt, s)) {
                intersected++;
                //cout << " line intersectd side i " << i <<  endl;
            }
        }
        if (intersected > 0) {
           //cout << "no of sides intersected " << intersected << endl;
           result = true;
        }
        return result;
    }





//private:

//distance between two points
double  getdist(point p1, point p2) {
   double result;
   result = sqrt((p1.first - p2.first)*(p1.first - p2.first) + (p1.second-p2.second)*(p1.second-p2.second));
   return result;
   }


};


