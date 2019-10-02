pair<float,float> v1 = make_pair (-0.28f, 0.0f); 
pair<float,float> v2 = make_pair (0.28f, 0.0f);
pair<float,float> v3 = make_pair (0.28f, 0.45f); 
pair<float,float> v4 = make_pair (-0.28f, 0.45f); 
BezCurve c1 = BezCurve(v1,v2);
BezCurve c2 = BezCurve(v2,v3);
BezCurve c3 = BezCurve(v3,v4);
BezCurve c4 = BezCurve(v4,v1);
BezCurvePath bound;
bound.addCurve(c1);
bound.addCurve(c2);
bound.addCurve(c3);
bound.addCurve(c4);

