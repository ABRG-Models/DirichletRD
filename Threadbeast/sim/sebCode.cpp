pair<float,float> v1 = make_pair (-0.28f, 0.0f); 
pair<float,float> v2 = make_pair (0.28f, 0.0f);
pair<float,float> v3 = make_pair (0.28f, 0.45f); 
pair<float,float> v4 = make_pair (-0.28f, 0.45f); 
cout << "ater making pairs" << endl;;
morph::BezCurve c1(v1,v2);
morph::BezCurve c2(v2,v3);
morph::BezCurve c3(v3,v4);
morph::BezCurve c4(v4,v1);
cout << "after making BezCurves" << endl;
morph::BezCurvePath bound;
cout << "after making BezCurvePath" << endl;
bound.addCurve(c1);
cout << "after adding curve" << endl;
bound.addCurve(c2);
cout << "after adding curve" << endl;
bound.addCurve(c3);
cout << "after adding curve" << endl;
bound.addCurve(c4);
cout << "after adding curve" << endl;
  Hgrid = new HexGrid(this->ds, 4.0, 0.0, morph::HexDomainShape::Boundary);
  cout << "after creating HexGrid"<<endl;
  //morph::ReadCurves r("./barrelAE.svg");
  //Hgrid->setBoundary (r.getCorticalPath());
  Hgrid->setBoundary (bound);
  cout << "before filling H " << Hgrid->num() << endl;
